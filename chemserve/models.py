from __future__ import annotations
from typing import Union, Optional, Iterable
from pathlib import Path

from chemserve._rdkit_imports import Chem, MoleculeConversionError
from chemserve._concrete_fns import ChemWrap
from chemserve.base_chem import BaseChem as _BaseChem


class Concrete(ChemWrap):
    @classmethod
    def from_mol_file(
        cls, path: Union[Path, str], name: Optional[str] = None, key: Optional[str] = None
    ) -> Concrete:
        mol = Chem.MolFromMolFile(str(path))
        return Concrete(mol, name, key)

    @classmethod
    def join(cls, sources: Iterable[Concrete], label_sep: Optional[str] = "; "):
        smiles = ".".join([source.smiles for source in sources])
        if label_sep is None:
            name, key = None, None
        else:
            name = f"[{[i.name for i in sources]}"
            key = f"[{[i.key for i in sources]}"
        return cls.of(smiles, name, key)

    @classmethod
    def of(
        cls,
        source: Union[Chem.Mol, str, ChemWrap],
        name: Optional[str] = None,
        key: Optional[str] = None,
    ) -> Concrete:
        if isinstance(source, ChemWrap):
            return source._copy_as(cls, replace=True, name=name, key=key)
        if isinstance(source, str) and source.startswith("InChI="):
            mol = Chem.MolFromInchi(source)
        elif isinstance(source, str):
            mol = Chem.MolFromSmiles(source)
        elif isinstance(source, Chem.Mol):
            mol = source
        else:
            raise TypeError(
                f"Type {type(source)} (for instance {source}, key {key}, name {name}) is not valid for mol"
            )
        if mol is None:
            raise MoleculeConversionError(
                f"{source} (key {key}, name {name}) was not converted correctly by SMILES or InChI"
            )
        return cls(mol, name, key)


class BaseChem(_BaseChem):
    def __init__(self, inchi_or_smiles: str, name: Optional[str] = None, key: Optional[str] = None):
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


__all__ = ["Concrete", "BaseChem"]
