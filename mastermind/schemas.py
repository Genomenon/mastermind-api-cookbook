# pylint: disable=E0611, R0903
from pydantic import BaseModel


class JobResponse(BaseModel):
    """JobResponse Schema
    The mastermind: file_annotations/counts api for job creation and polling
    returns the following output via json. This is used for extracting
    and validating the data returned"""

    job_id: str
    state: str
    assembly: str
    upload_url: str
    download_url: str
    input_filename: str
    created_at: str
    job_url: str
