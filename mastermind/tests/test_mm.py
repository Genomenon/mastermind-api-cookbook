import os
import pytest
from mastermind.api import MasterMind
from mastermind.cred import get_creds, CredSchema


cd = os.path.dirname(os.path.realpath(__file__))
yaml_file = os.path.join(cd, "data", "mastermind.yaml")
vcf = os.path.join(cd, "data", "some_vcf.vcf")
out_vcf = os.path.join(cd, "data", "out", "out_vcf.vcf")


def test_get_creds():
    cred: CredSchema = get_creds(yaml_file)
    assert cred
    assert cred.token == "HEREISMYTOKEN"


def test_init_mm_slot_test():
    cred: CredSchema = get_creds(yaml_file)
    mm = MasterMind(token=cred.token)

    assert mm.token
    assert mm.mm_api_url == "https://mastermind.genomenon.com/api/v2"
    assert mm.assembly == "grch38"
    with pytest.raises(AttributeError):
        mm.dog = False
    mm.assembly = "grch37"


@pytest.mark.functional
def test_create_job(mmapi: MasterMind):
    resp = mmapi.create_job("somefile")
    assert resp.job_id
    assert resp.state == "created"
    assert resp.assembly
    assert resp.upload_url
    assert resp.input_filename
    assert resp.job_url


@pytest.mark.functional
def test_vcf_annotate(mmapi: MasterMind):
    mmapi.vcf_annotate(vcf_path=vcf, out_vcf_path=out_vcf)
    assert mmapi.curr_job == None
