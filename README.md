# Mastermind API Cookbook

This repository consists of sample scripts for interacting with the
Mastermind RESTful JSON API.

Each Python file in this repository contains a description of what the
file does and how to run it at the top of the file.

For more information about the API for the Mastermind Genomic Search
Engine, see [the live API demo and
documentation](https://mastermind.genomenon.com/api/).

[Contact Us](https://www.genomenon.com/contact/) to get your API key to
get started.

# Testing

```bash
pip install -r requirements.txt
pytest # Unit tests
pytest --functional # Functional tests
```

# Api 

Simple API wrapper 

```
mmapi = MasterMind(token=token, mm_api_url=mm_api_url, assembly=assembly)
mmapi.vcf_annotate(vcf_path=vcf, out_vcf_path=out_vcf)
```


# Cli Tool

```bash
pip install -e . 
mastermind   
Usage: mastermind [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  annotate-with-evidence  Annotate a VCF with Evidence
```



