"""
Metadata for chemgrid.
"""
import logging
from pathlib import Path
# importlib.metadata is compat with Python 3.8 only
from importlib_metadata import PackageNotFoundError, metadata as __load

from chemgrid.graphics import ChemGraphicsKit
from chemgrid.payload import Payload, ConcretePayload
from chemgrid.service import Client, Server, ChemGrid
from chemgrid.models import BaseChem, Concrete as Chem

logger = logging.getLogger('chemgrid')

try:
    metadata = __load(Path(__file__).parent.name)
    __status__ = "Development"
    __copyright__ = "Copyright 2020"
    __date__ = "2020-06-03"
    __uri__ = metadata["home-page"]
    __title__ = metadata["name"]
    __summary__ = metadata["summary"]
    __license__ = metadata["license"]
    __version__ = metadata["version"]
    __author__ = metadata["author"]
    __maintainer__ = metadata["maintainer"]
    __contact__ = metadata["maintainer"]
except PackageNotFoundError:
    logger.error("Could not import package metdata for chemgrid. Is it not installed?")
