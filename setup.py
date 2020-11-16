from setuptools import setup

setup(
    name="mastermind",
    version="0.0.1",
    description="Testing Mastermind",
    url="",
    author="",
    author_email="chris.ieng@vcgs.org.au",
    license="",
    packages=["mastermind"],
    # install_requires=requirements,
    entry_points={"console_scripts": ["mastermind = mastermind.cli:cli"]},
    zip_safe=False,
)
