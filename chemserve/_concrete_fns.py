from __future__ import annotations
from typing import Iterable, Union
from pathlib import Path
from copy import copy

from chemserve._rdkit_imports import Chem, SaltRemover, logger
from chemserve._concrete_base import ChemWrap as _Wrap
from chemserve.fingerprints import Fingerprint


class ChemWrap(_Wrap):
    """"""

    @property
    def formal_charge(self) -> int:
        return Chem.rdmolops.GetFormalCharge(self._mol)

    @property
    def formula(self):
        return self._descriptor("CalcMolFormula")

    @property
    def n_stereocenters(self) -> int:
        return self._descriptor("CalcNumAtomStereoCenters")

    @property
    def n_rings(self) -> int:
        return self._descriptor("CalcNumRings")

    @property
    def n_aromatic_rings(self) -> int:
        return self._descriptor("CalcNumAromaticRings")

    @property
    def n_unspecified_stereocenters(self) -> int:
        return self._descriptor("CalcNumUnspecifiedAtomStereoCenters")

    def _descriptor(self, fn_name):
        if fn_name not in self._metrics:
            fn = getattr(Chem.rdMolDescriptors, fn_name)
            self._metrics[fn_name] = fn(self._mol)
        return self._metrics[fn_name]

    def write_sdf(self, path: Union[str, Path]) -> None:
        w = Chem.SDWriter(path)
        w.write(self._mol)

    def ecfp(self, radius: int, bits: int = 2048) -> Fingerprint:
        """
        Calculates an ECFP (Morgan) fingerprint.

        Args:
            radius: Recommended values are 2, 3, and 4
            bits: The number of bits to encode the fingerprint with
        """
        return self._morgan(radius=radius, nBits=bits, useFeatures=False)

    def ecfp_with_fcfp(self, radius: int, bits: int = 2048) -> Fingerprint:
        """
        Calculates a fingerprint that combines ECFP with "functional-class fingerprint" features.
        See `the rdkit FPCP features <https://rdkit.readthedocs.io/en/latest/GettingStartedInPython.html#feature-definitions-used-in-the-morgan-fingerprints>`_.

        Args:
            radius: Recommended values are 2, 3, and 4
            bits: The number of bits to encode the fingerprint with
        """
        return self._morgan(radius=radius, nBits=bits, useFeatures=True)

    def _morgan(self, **kwargs) -> Fingerprint:
        from rdkit.Chem import AllChem

        fp1 = AllChem.GetMorganFingerprintAsBitVect(self._mol, **kwargs)
        return Fingerprint(fp1)

    def is_substruct_of(self, f: ChemWrap) -> bool:
        return self.mol.HasSubstructMatch(f.mol)

    def without_stereochem(self) -> __qualname__:
        mol = copy(self._mol)
        Chem.rdmolops.RemoveStereochemistry(mol)
        return ChemWrap(mol, self._name, self._key)

    def with_2d_coords(self, hydrogens: bool = True) -> __qualname__:
        return self._with_coords(False, hydrogens=hydrogens)

    def with_3d_coords(self, hydrogens: bool = True) -> __qualname__:
        return self._with_coords(True, hydrogens=hydrogens)

    def _with_coords(self, three_d: bool, hydrogens: bool) -> __qualname__:
        if hydrogens:
            mol = Chem.AddHs(self._mol)
        else:
            mol = copy(self._mol)
        # noinspection PyUnresolvedReferences
        from rdkit.Chem import AllChem

        if three_d:
            AllChem.EmbedMolecule(mol)
        else:
            AllChem.Compute2DCoords(mol)
        return ChemWrap(mol, name=self._name, key=self._key)

    def desalt(self) -> __qualname__:
        # don't make a copy because the inchi might be different
        desalted = SaltRemover.SaltRemover().StripMol(self._mol)
        return self.__class__(desalted, self._name, self._key)

    def deduplicate(self) -> __qualname__:
        parts = set(self.smiles.split("."))
        joined = ".".join(parts)
        mol = Chem.MolFromSmiles(joined)
        return self.__class__(mol, self.name, self.key)

    def longest_component(self) -> __qualname__:
        parts = set(self.smiles.split("."))
        if len(parts) > 1:
            logger.debug(
                f"SMILES {self.smiles} (name={self.name}, key={self.key}) has {len(parts)} unique parts"
            )
        longest = ""
        for p in parts:
            if len(p) > len(longest):
                longest = p
        mol = Chem.MolFromSmiles(longest)
        return ChemWrap(mol, self.name, self.key)

    @classmethod
    def mcs(cls, items: Iterable[ChemWrap], **kwargs) -> __qualname__:
        """
        Finds a maximum common substructure using rdkit ``Chem.rdFMCS``.

        Args:
            items: Molecules
            kwargs: Passed directly to ``rdFMCS``
        """
        return cls(
            Chem.rdFMCS([i.mol for i in items], **kwargs),
            name=f"[{[i.name for i in items]}",
            key=f"[{[i.key for i in items]}",
        )


__all__ = ["ChemWrap"]
