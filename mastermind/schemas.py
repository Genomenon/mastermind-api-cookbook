from pydantic import BaseModel


class JobResponse(BaseModel):
    job_id: str
    state: str
    assembly: str
    upload_url: str
    download_url: str
    input_filename: str
    created_at: str
    job_url: str
