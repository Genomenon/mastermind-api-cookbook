import click
from mastermind.api import MasterMind, DEFAULT_API_URL


def print_help() -> None:
    ctx = click.get_current_context()
    click.echo(ctx.get_help())
    ctx.exit()


@click.command(help="Annotate a VCF with Evidence")
@click.option("--token", help="Mastermind API Token")
@click.option("--assembly", default="grch38", help="Genome Build")
@click.option("--mm_api_url", default=DEFAULT_API_URL, help="Genome Build")
@click.option("--vcf", help="VCF to annotate")
@click.option("--out_vcf", help="Path to output")
def annotate_with_evidence(token, assembly, mm_api_url, vcf, out_vcf):
    if not all([token, vcf, out_vcf]):
        print_help()
    mmapi = MasterMind(token=token, mm_api_url=mm_api_url, assembly=assembly)
    mmapi.vcf_annotate(vcf_path=vcf, out_vcf_path=out_vcf, assembly=assembly)


@click.group()
def cli():
    pass


cli.add_command(annotate_with_evidence)


if __name__ == "__main__":
    cli()
