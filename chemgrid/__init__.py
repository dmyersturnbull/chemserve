"""
Metadata for chemgrid-service.
"""

from pathlib import Path, PurePath
from typing import Union, Iterable, Optional

# importlib.metadata is compat with Python 3.8 only
from importlib_metadata import metadata as __load

metadata = __load(Path(__file__).parent.name)
__status__ = "Development"
__copyright__ = "Copyright 2020"
__date__ = "2020-05-24"
__uri__ = metadata["home-page"]
__title__ = metadata["name"]
__summary__ = metadata["summary"]
__license__ = metadata["license"]
__version__ = metadata["version"]
__author__ = metadata["author"]
__maintainer__ = metadata["maintainer"]
__contact__ = metadata["maintainer"]


def resource(*nodes: Union[PurePath, str]) -> Path:
    """Gets a path of a resource file under resources/ directory."""
    return Path(Path(__file__).parent, "resources", *nodes)


class QueryBuilder:
    def __init__(self):
        self._rows = 8
        self._cols = 12
        self._compounds = []
        self._labels = []
        self._desalt = True
        self._deduplicate = True
        self._path = None
        
    def grid(self, rows: int, cols: int):
        self._rows = rows
        self._cols = cols
        return self
    
    def add(self, inchi_inchikey_smiles: str, label: Optional[str] = None):
        self._compounds.append(inchi_inchikey_smiles)
        self._labels.append(label)
        return self
    
    def keep_duplicates(self):
        self._deduplicate = False
        return self
        
    def keep_salts(self):
        self._desalt = False
        return self
        
    def to(self, path: Union[Path, str]):
        self._path = path
        return self
        
    def payload():
        pass
    
    def send():
        pass


class Chemgrid:

    @classmethod
    def is_alive() -> bool:
        pass

    @classmethod
    def new(cls, rows: int = 8, cols: int = 12) -> QueryBuilder:
        return QueryBuilder().rows(rows).cols(cols)

    @classmethod
    def send(
        cls,
        json_data: dict
    ):
        pass
