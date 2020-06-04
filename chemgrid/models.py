
from __future__ import annotations
from typing import Union, Optional, Iterable
from pathlib import Path

from chemgrid._rdkit_imports import Chem, MoleculeConversionError
from chemgrid._concrete_fns import ChemWrap
from chemgrid.base_chem import BaseChem as _BaseChem


class Concrete(ChemWrap):

    @classmethod
    def from_mol_file(
        cls, path: Union[Path, str], name: Optional[str], key: Optional[str]
    ) -> Concrete:
        mol = Chem.MolFromMolFile(str(path))
        return Concrete(mol, name, key)

    @classmethod
    def join(cls, sources: Iterable[Concrete], label_sep: Optional[str] = "; "):
        smiles = ".".join([source.smiles for source in sources])
        if label_sep is None:
            name, key = None, None
        else:
            name = (
                None
                if any([s.name is None for s in sources])
                else "; ".join([s.name for s in sources])
            )
            key = (
                None
                if any([s.key is None for s in sources])
                else "; ".join([s.key for s in sources])
            )
        return cls.of(smiles, name, key)

    @classmethod
    def of(
        cls, source: Union[Chem.Mol, str, Concrete], name: Optional[str], key: Optional[str]
    ) -> Concrete:
        if isinstance(source, Concrete):
            return source.copy(name=name, key=key)
        if isinstance(source, str) and source.startswith("InChI="):
            mol = Chem.MolFromInchi(source)
        elif isinstance(source, str):
            mol = Chem.MolFromSmiles(source)
        elif isinstance(source, Chem.Mol):
            mol = source
        else:
            raise TypeError(
                "Type {} (for instance {}, key {}, name {}) is not valid for mol".format(
                    type(source), source, key, name
                )
            )
        if mol is None:
            raise MoleculeConversionError(
                "{} (key {}, name {}) was not converted correctly by SMILES or InChI".format(
                    source, key, name
                )
            )
        return Concrete(mol, name, key)


class BaseChem(_BaseChem):

    def __init__(self, inchi_or_smiles: str, name: Optional[str], key: Optional[str]):
        self._inchi_or_smiles = inchi_or_smiles
        self._name = name
        self._key = key

    @property
    def inchi_or_smiles(self) -> str:
        return self._inchi_or_smiles

    @property
    def name(self) -> str:
        return self._name

    @property
    def key(self) -> str:
        return self._key


__all__ = ['Concrete', 'BaseChem']
