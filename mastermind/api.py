"""
Mastermind wrapper class for api calls
"""
import os
import requests
import time
from io import BufferedReader
from typing import List, Optional, Union
from .schemas import JobResponse
from time import perf_counter
from contextlib import contextmanager

@contextmanager
def measure_time() -> float:
    start = perf_counter()
    yield lambda: perf_counter() - start

WAIT_INTERVAL_SECONDS = 10
DEFAULT_API_URL = "https://mastermind.genomenon.com/api/v2"

class MasterMind:

    __slots__ = ["mm_api_url", "token", "assembly", "curr_job"]

    __endpoints__ = [
        "suggestions",
        "counts",
        "articles",
        "genes",
        "variants",
        "diseases",
        "phenotypes",
        "therapies",
        "article_info",
        "file_annotations/counts",
    ]

    def __init__(
        self,
        token: str,
        mm_api_url: str = DEFAULT_API_URL,
        assembly: str = "grch38",
    ):
        self.mm_api_url = mm_api_url
        self.token = token
        self.assembly = assembly
        self.curr_job: Optional[JobResponse] = None

    def api_request(
        self,
        url: Optional[str] = None,
        endpoint: Optional[str] = None,
        params: Optional[dict] = None,
        request_type: str = "GET",
        extra_headers: dict = None,
        data: Optional[Union[dict, BufferedReader]] = None,
    ) -> requests.Response:
        if endpoint and not url:
            if not any([x in endpoint for x in self.__endpoints__]):
                raise ValueError("Invalid endpoint given %s", endpoint)
            url = append_url(self.mm_api_url, endpoint)
        elif url and not endpoint:
            url = url
        else:
            raise ValueError("Cannot supply both endpoint and url")
        
        header = {"X-API-TOKEN": self.token}
        if extra_headers:
            header.update(extra_headers)
        if not data:
            response = requests.request(
                request_type, url=url, params=params, headers=header
            )
        else:
            response = requests.request(
                request_type, data=data, url=url, params=params, headers=header
            )
        return response

    def vcf_annotate(self, vcf_path: str, out_vcf_path: str, assembly: str = None) -> None:
        """Given a VCF input use the mastermind apis to annotate the vcf.

        There are Multiple steps for annotating a Mastermind VCF
        1. Submit a job to genomenon
        """
        # check exists
        if not os.path.exists(vcf_path):
            raise ValueError("cannot find %s", vcf_path)
        file_name = os.path.basename(vcf_path)
        if not assembly:
            assembly = self.assembly

        if not self.curr_job:
            self.curr_job = self.create_job(file_name, assembly)

        if self.curr_job.state == "created":
            self.upload_vcf(self.curr_job.upload_url, vcf_path)

        # We need to keep pinging the job until it's succeeded
        # Let's also measure how long this takes
        with measure_time() as t:
            # Keep checking until state is succeeded
            while self.curr_job.state != "succeeded":
                self.check_job()
                time.sleep(WAIT_INTERVAL_SECONDS)
                print(f"\nTime elapsed: {t():.4f} secs\n")
                print(f"Job state: {self.curr_job.state}")
                if self.curr_job.state == "failed":
                    raise ValueError("VCF Annotate failed")
            
        print(f"\nFinal Time elapsed: {t():.4f} secs\n")
        self.download_vcf(out_vcf_path=out_vcf_path)


    def upload_vcf(self, upload_url: str, vcf_path: str) -> None:
        with open(vcf_path, "rb") as open_f:
            response = self.api_request(
                url=upload_url,
                request_type="PUT",
                extra_headers={"Content-Type": "application/octet-stream"},
                data=open_f,
            )
            if response.status_code != 200:
                raise ValueError("Could not upload vcf. Response code %d", response.status_code)
    
    def download_vcf(self, out_vcf_path: str) -> None:
        if not self.curr_job:
            raise ValueError("No job set. Cannot download vcf")
        if out_vcf_path[-3:] != ".gz":
            out_vcf_path = out_vcf_path + ".gz"
        response = self.api_request(endpoint=f"file_annotations/counts/{self.curr_job.job_id}/download")
        if response.status_code != 200:
            raise ValueError("Could not download annotated vcf. Status code %d", response.status_code)
        with open(out_vcf_path, "wb") as out_f:
            out_f.write(response.content)
        self.curr_job = None 


    def create_job(self, filename: str, assembly: str = None) -> JobResponse:
        if not assembly:
            assembly = self.assembly
        params = {"assembly": assembly, "filename": filename}
        response = self.api_request(
            endpoint="file_annotations/counts", params=params, request_type="POST"
        )
        return JobResponse(**response.json())
    
    def check_job(self) -> None:
        if not self.curr_job:
            raise ValueError("No job to check")
        response = self.api_request(endpoint=f"file_annotations/counts/{self.curr_job.job_id}")
        # Update the curr job
        self.curr_job = JobResponse(**response.json())


def fix_url(part):
    if part[0] == "/":
        # Remove initial / if exists
        part = part[1:]
    if part[-1] == "/":
        return part
    else:
        return part + "/"


def append_url(baseurl: str, *args: List[str]) -> str:
    built_url: str = fix_url(baseurl)
    for arg in args:
        built_url += fix_url(arg)
    # Remove last /?
    built_url = built_url[0:-1]
    return built_url
