import os
import re
import time
import json
import zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry


class UniProtAPIClient:
    POLLING_INTERVAL = 3
    API_URL = "https://rest.uniprot.org"

    def __init__(self):
        retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(max_retries=retries))

    def check_response(self, response):
        """
        Ensure the response is valid, otherwise raise an error.
        """
        try:
            response.raise_for_status()
        except requests.HTTPError:
            print(response.json())
            raise

    def submit_id_mapping(self, from_db, to_db, ids):
        """
        Submit an ID mapping job to UniProt.
        :param from_db: Database to map from (e.g., UniProtKB_AC-ID)
        :param to_db: Database to map to (e.g., ChEMBL)
        :param ids: List of IDs to map
        :return: jobId for tracking the status of the job
        """
        request = self.session.post(
            f"{self.API_URL}/idmapping/run",
            data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
        )
        self.check_response(request)
        return request.json()["jobId"]

    def get_next_link(self, headers):
        """
        Retrieve the 'next' link from the headers for pagination.
        :param headers: Response headers
        :return: URL of the next page, or None if not present
        """
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)
        return None

    def check_id_mapping_results_ready(self, job_id):
        """
        Poll the API to check if the ID mapping job is completed.
        :param job_id: The job ID for tracking
        :return: True if the results are ready, False otherwise
        """
        while True:
            request = self.session.get(f"{self.API_URL}/idmapping/status/{job_id}")
            self.check_response(request)
            job_status = request.json()
            if "jobStatus" in job_status:
                if job_status["jobStatus"] in ("NEW", "RUNNING"):
                    print(f"Job {job_id} is running. Retrying in {self.POLLING_INTERVAL}s")
                    time.sleep(self.POLLING_INTERVAL)
                else:
                    raise Exception(f"Job failed with status: {job_status['jobStatus']}")
            else:
                return bool(job_status.get("results") or job_status.get("failedIds"))

    def get_id_mapping_results_link(self, job_id):
        """
        Retrieve the link to the ID mapping job results.
        :param job_id: The job ID
        :return: The URL to access the results
        """
        url = f"{self.API_URL}/idmapping/details/{job_id}"
        request = self.session.get(url)
        self.check_response(request)
        return request.json()["redirectURL"]

    def decode_results(self, response, file_format, compressed):
        """
        Decode the results based on the format and compression.
        :param response: The API response
        :param file_format: The format of the results (json, tsv, etc.)
        :param compressed: Whether the results are compressed
        :return: Decoded results
        """
        if compressed:
            decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
            if file_format == "json":
                return json.loads(decompressed.decode("utf-8"))
            elif file_format == "tsv":
                return [line for line in decompressed.decode("utf-8").split("\n") if line]
        else:
            if file_format == "json":
                return response.json()
            elif file_format == "tsv":
                return [line for line in response.text.split("\n") if line]
        return response.text

    def get_batch(self, batch_response, file_format, compressed):
        """
        Fetch the next batch of results if pagination is needed.
        :param batch_response: Initial response
        :param file_format: The format of the results
        :param compressed: Whether the results are compressed
        :return: Generator of batch results
        """
        batch_url = self.get_next_link(batch_response.headers)
        while batch_url:
            batch_response = self.session.get(batch_url)
            batch_response.raise_for_status()
            yield self.decode_results(batch_response, file_format, compressed)
            batch_url = self.get_next_link(batch_response.headers)

    def combine_batches(self, all_results, batch_results, file_format):
        """
        Combine multiple batches of results into a single dataset.
        :param all_results: All previous results
        :param batch_results: New batch results
        :param file_format: Format of the results
        :return: Combined results
        """
        if file_format == "json":
            for key in ("results", "failedIds"):
                if key in batch_results and batch_results[key]:
                    all_results[key] += batch_results[key]
        elif file_format == "tsv":
            return all_results + batch_results[1:]
        return all_results

    def get_id_mapping_results_search(self, url):
        """
        Fetch the ID mapping results in a paginated manner.
        :param url: URL for retrieving the results
        :return: Final combined results
        """
        parsed = urlparse(url)
        query = parse_qs(parsed.query)
        file_format = query["format"][0] if "format" in query else "json"
        query["size"] = [500]  # Fetch the maximum number of entries per page
        compressed = query.get("compressed", ["false"])[0].lower() == "true"
        
        parsed = parsed._replace(query=urlencode(query, doseq=True))
        url = parsed.geturl()

        request = self.session.get(url)
        self.check_response(request)
        results = self.decode_results(request, file_format, compressed)

        total = int(request.headers.get("x-total-results", 500))
        fetched = len(results["results"]) if "results" in results else 0

        while fetched < total:
            next_link = self.get_next_link(request.headers)
            if not next_link:
                break
            request = self.session.get(next_link)
            self.check_response(request)
            batch_results = self.decode_results(request, file_format, compressed)
            results = self.combine_batches(results, batch_results, file_format)
            fetched = len(results["results"])

        return results
    
    def map_entries(self, from_db, to_db, ids):
        """
        Map IDs between two databases and fetch the results.
        :param from_db: Source database
        :param to_db: Target database
        :param ids: List of IDs to map
        :return: Mapped results
        """
        job_id = self.submit_id_mapping(from_db, to_db, ids)
        if self.check_id_mapping_results_ready(job_id):
            link = self.get_id_mapping_results_link(job_id)
            return self.get_id_mapping_results_search(link)

    def write_results(self, results, outdir, batch_idx):
        """
        Write the results to a file in JSON format.
        :param results: Results to write
        :param outdir: Output directory
        :param batch_idx: Index of the batch
        """
        with open(f"{outdir}/results_{batch_idx}.json", "w") as f:
            json.dump(results, f, indent=4)

    def batch_search(self, from_db, to_db, ids, outdir):
        """
        Perform a batch search for multiple protein IDs.
        :param from_db: Source database
        :param to_db: Target database
        :param ids: List of IDs to search
        :param outdir: Output directory
        :return: Dictionary with counts of results and failed IDs
        """
        os.makedirs(outdir, exist_ok=True)
        # Batches of 500 proteins
        batches = [ids[i:i+500] for i in range(0, len(ids), 500)]
        failed = []
        for batch_idx, batch in enumerate(batches):
            try:
                results = self.map_entries(from_db, to_db, batch)
                if results:
                    self.write_results(results, outdir, batch_idx)
                else:
                    failed += batch
            except Exception as e:
                failed += batch
                print(e)

        # Write failed ids
        with open(f"{outdir}/failed_ids.txt", "w") as f:
            f.write("\n".join(failed))

        return {"results": len(os.listdir(outdir)), "failedIds": len(failed)}

    def single_entity_search(self, from_db, to_db, protein_id):
        """
        Perform a search for a single protein ID.
        :param from_db: Source database
        :param to_db: Target database
        :param protein_id: Single protein ID to search
        :return: Mapped result
        """
        try:
            result = self.map_entries(from_db, to_db, [protein_id])
            if result and result.get("results"):
                return result["results"]
            else:
                print(f"No results found for protein ID: {protein_id}")
        except Exception as e:
            print(f"Error occurred during search for protein ID {protein_id}: {str(e)}")
