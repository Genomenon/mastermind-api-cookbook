"""
Get credentials for mastermind API
"""

import os
from ruamel.yaml import YAML
from pydantic import BaseModel


class CredSchema(BaseModel):
    token: str


def get_creds(
    yaml_file: str = os.path.expanduser("~/mastermind.yaml"),
) -> CredSchema:
    if not os.path.exists(yaml_file):
        raise ValueError(
            "Path to yaml_file doesn't exist: %s. Please supply one with the token inside",
            yaml_file,
        )

    parser = YAML(typ="safe")
    data = None
    with open(yaml_file) as of:
        data = of.read()
    data = parser.load(data)
    cred = CredSchema(**data)
    return cred
