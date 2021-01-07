"""
Metadata for chemserve.
"""
import logging
from pathlib import Path
from importlib.metadata import PackageNotFoundError
from importlib.metadata import metadata as __load

from chemserve.graphics import ChemGraphicsKit
from chemserve.payload import Payload, ConcretePayload
from chemserve.service import Client, Server, ChemServe
from chemserve.models import BaseChem, Concrete as Chem

pkg = Path(__file__).absolute().parent.name
logger = logging.getLogger(pkg)
_metadata = None
try:
    _metadata = __load(Path(__file__).absolute().parent.name)
    __status__ = "Development"
    __copyright__ = "Copyright 2020–2021"
    __date__ = "2020-06-03"
    __uri__ = _metadata["home-page"]
    __title__ = _metadata["name"]
    __summary__ = _metadata["summary"]
    __license__ = _metadata["license"]
    __version__ = _metadata["version"]
    __author__ = _metadata["author"]
    __maintainer__ = _metadata["maintainer"]
    __contact__ = _metadata["maintainer"]
except PackageNotFoundError:  # pragma: no cover
    logger.error(f"Could not load package metadata for {pkg}. Is it installed?")
