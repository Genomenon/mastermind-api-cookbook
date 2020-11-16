import pytest
from mastermind.api import MasterMind
from mastermind.cred import get_creds


@pytest.fixture(scope="session")
def mmapi():
    creds = get_creds()
    mm_api = MasterMind(token=creds.token)
    return mm_api
