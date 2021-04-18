"""
CLI.
"""

import typer


cli = typer.Typer()


@cli.command()
def start():
    """
    Starts a server that accepts arbitrary Python code.
    """
    pass


if __name__ == "__main__":
    cli()
