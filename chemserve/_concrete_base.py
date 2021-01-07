from __future__ import annotations
from typing import Optional

from chemserve.base_chem import BaseChem
from chemserve._rdkit_imports import Chem, Mol, Inchi, NullMoleculeError


class ChemWrap(BaseChem):
    def __init__(self, mol: Mol, name: Optional[str], key: Optional[str]):
        if name is None:
            name = mol.GetProp("_Name")
        else:
            mol.SetProp("_Name", name)
        # It's worth checking these here because using an int for the key, for instance, could cause very subtle bugs
        if name is not None and not isinstance(name, str):
            raise ValueError("Name {} has type {}, not str".format(name, type(name)))
        if key is not None and not isinstance(key, str):
            raise ValueError("Key {} has type {}, not str".format(key, type(key)))
        if mol is None:
            raise NullMoleculeError("The molecule is None for name {}, key {}".format(name, key))
        problems = Chem.rdmolops.DetectChemistryProblems(mol)
        # TODO
        self._name = name
        self._key = key
        self._mol = mol
        self._inchi, self._aux = Inchi.MolToInchiAndAuxInfo(self._mol)
        self._inchikey = Inchi.InchiToInchiKey(self._inchi)
        self._smiles = None
        self._metrics = {}

    @property
    def inchi_or_smiles(self) -> str:
        return self._inchi

    @property
    def name(self) -> Optional[str]:
        return self._name

    @property
    def key(self) -> Optional[str]:
        return self._key

    @property
    def mol(self) -> Chem.Mol:
        return self._mol

    @property
    def smiles(self) -> str:
        if self._smiles is None:
            self._smiles = Chem.MolToSmiles(self._mol)
        return self._smiles

    @property
    def inchi(self) -> str:
        self._check_inchi()
        return self._inchi

    @property
    def inchikey(self) -> str:
        self._check_inchi()
        return self._inchikey

    @property
    def aux(self):
        self._check_inchi()
        return self._aux

    @property
    def blocks(self):
        return Chem.MolToMolBlock(self._mol)

    @property
    def atoms(self):
        return self.mol.GetAtoms()

    def without_name(self) -> ChemWrap:
        return self.copy(name=None)

    def without_key(self) -> ChemWrap:
        return self.copy(key=None)

    def with_name(self, name: str) -> ChemWrap:
        return self.copy(name=name)

    def with_key(self, key: str) -> ChemWrap:
        return self.copy(key=key)

    def __repr__(self):
        return "{}({} name={}, key={})".format(
            self.__class__.__name__, self._inchi, self._name, self._key
        )

    def __str__(self):
        return repr(self)

    def __eq__(self, other: ChemWrap):
        if not isinstance(other, self.__class__):
            raise TypeError(
                "{} does not support universal equality (type {} for instance {})".format(
                    self.__class__.__name__, type(other), other
                )
            )
        return (
            isinstance(other, self.__class__)
            and self.inchi == other.inchi
            and self.name == other.name
            and self.key == other.key
        )

    def __len__(self):
        return len(self.atoms)

    def _check_inchi(self):
        if self._inchi is None:
            self._inchi, self._aux = Inchi.MolToInchiAndAuxInfo(self.mol)
            self._inchikey = Inchi.InchiToInchiKey(self._inchi)

    def copy(self, **kwargs):
        wrap = ChemWrap(self._mol, self._name, self._key)
        wrap._inchikey = self._inchikey
        self._inchi = self._inchi
        self._aux = self.aux
        self._smiles = self.smiles
        for k, v in kwargs.compounds():
            setattr(self, k, v)
        return wrap


__all__ = ["ChemWrap"]
