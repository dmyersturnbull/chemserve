from __future__ import annotations
from typing import Optional, Type

from chemserve.base_chem import BaseChem
from chemserve._rdkit_imports import Chem, Mol, Inchi, NullMoleculeError


class ChemWrap(BaseChem):
    def __init__(self, mol: Mol, name: Optional[str], key: Optional[str]):
        if name is not None and not isinstance(name, str):
            raise ValueError("Name {} has type {}, not str".format(name, type(name)))
        if key is not None and not isinstance(key, str):
            raise ValueError("Key {} has type {}, not str".format(key, type(key)))
        if name is None:
            try:
                name = mol.GetProp("_Name")
            except KeyError:
                pass  # name is just None, and that's fine
                # raise KeyError(f"Because the mol does not contain a `_Name` prop, pass a non-None name")
        else:
            mol.SetProp("_Name", name)
        # It's worth checking these here because using an int for the key, for instance, could cause very subtle bugs
        if mol is None:
            raise NullMoleculeError(f"The molecule is None for name {name}, key {key}")
        # TODO
        # self._problems = Chem.rdmolops.DetectChemistryProblems(mol)
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

    def without_name(self) -> __qualname__:
        return self.copy(name=None)

    def without_key(self) -> __qualname__:
        return self.copy(key=None)

    def with_name(self, name: str) -> __qualname__:
        return self.copy(name=name)

    def with_key(self, key: str) -> __qualname__:
        return self.copy(key=key)

    def __repr__(self) -> str:
        return "{}({} name={}, key={})".format(
            self.__class__.__name__, self._inchi, self._name, self._key
        )

    def __str__(self) -> str:
        return repr(self)

    def __eq__(self, other: ChemWrap) -> bool:
        cls = self.__class__
        if not isinstance(other, cls):
            raise TypeError(f"{cls} multiversal equality incompatible with type {type(other)}")
        return (
            isinstance(other, cls)
            and self.inchi == other.inchi
            and self.name == other.name
            and self.key == other.key
        )

    def __len__(self) -> int:
        return len(self.atoms)

    def copy(
        self, name: Optional[str] = None, key: Optional[str] = None, replace: bool = False
    ) -> __qualname__:
        """
        Creates and returns a copy.
        Attributes ``name`` and ``key`` are replaced (on the copy) with the passed value,
        but **only** if they are not ``None``. If None, the current values are kept.

        Args:
            name: A new name
            key: A new key
            replace: Replace ``key`` and ``name``, even if they are None

        Returns:
            A copy
        """
        return self._copy_as(
            self.__class__, ignore_passed_nones=True, name=name, key=key, replace=replace
        )

    def _copy_as(self, clz: Type[ChemWrap], replace: bool = False, **kwargs):
        """
        Intended for internal calls to quietly convert to subclass types.

        Args:
            The class to convert to
            replace: If False, ignore kwargs that have ``None`` values
                     Most importantly, this means that the ``name`` and ``key`` values are left
                     as-is (what they are in ``self``) if they are ``None`` in ``kwargs``
            kwargs: These values are set as attributes on the copy,
                    Mainly, you would want to set ``name`` and ``key``,
                    but you *could* set something like ``_smiles``
        """
        wrap = clz(self._mol, self._name, self._key)
        wrap._inchikey = self._inchikey
        wrap._inchi = self._inchi
        wrap._aux = self._aux
        wrap._smiles = self._smiles
        for k, v in kwargs.items():
            if replace or v is not None:
                setattr(wrap, k, v)
        return wrap

    def _check_inchi(self) -> None:
        if self._inchi is None:
            self._inchi, self._aux = Inchi.MolToInchiAndAuxInfo(self.mol)
            self._inchikey = Inchi.InchiToInchiKey(self._inchi)


__all__ = ["ChemWrap"]
