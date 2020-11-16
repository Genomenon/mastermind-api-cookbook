import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--functional", action="store_true", default=False, help="run functional tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "functional: mark test as functional to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--functional"):
        # --functional given in cli: do not skip slow tests
        return
    skip_functional = pytest.mark.skip(reason="need --functional option to run")
    for item in items:
        if "functional" in item.keywords:
            item.add_marker(skip_functional)
