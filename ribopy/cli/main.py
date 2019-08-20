# -*- coding: utf-8 -*-
import click

@click.group()
@click.option('--debug/--no-debug', default = False)
def cli(debug):
    pass
    
from .create   import *
from .merge    import *
from .info     import *
from .dump     import *
from .plot     import *
from .metadata import *
from .rnaseq    import *



if __name__ == "__main__":

    cli()


