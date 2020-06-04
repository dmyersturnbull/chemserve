from __future__ import annotations
from typing import Union
from pathlib import Path
from copy import copy

from chemgrid._rdkit_imports import Chem, SaltRemover, logger
from chemgrid._concrete_base import ChemWrap as _Wrap


class ChemWrap(_Wrap):

    @property
    def formal_charge(self) -> int:
        return Chem.rdmolops.GetFormalCharge(self._mol)

    @property
    def formula(self):
        return self._descriptor('CalcMolFormula')

    @property
    def n_stereocenters(self) -> int:
        return self._descriptor('CalcNumAtomStereoCenters')

    @property
    def n_rings(self) -> int:
        return self._descriptor('CalcNumRings')

    @property
    def n_aromatic_rings(self) -> int:
        return self._descriptor('CalcNumAromaticRings')

    @property
    def n_unspecified_sterocenters(self) -> int:
        return self._descriptor('CalcNumUnspecifiedAtomStereoCenters')

    def _descriptor(self, fn_name):
        if fn_name not in self._metrics:
            fn = getattr(Chem.rdMolDescriptors, fn_name)
            self._metrics[fn_name] = fn(self._mol)
        return self._metrics[fn_name]

    def write_sdf(self, path: Union[str, Path]) -> None:
        w = Chem.SDWriter(path)
        w.write(self._mol)

    def fingerprint(self, radius: int, features: bool) -> str:
        # noinspection PyUnresolvedReferences
        from rdkit.Chem import AllChem
        fp1 = AllChem.GetMorganFingerprint(self._mol, radius, useFeatures=features)
        return fp1

    def without_sterochem(self) -> ChemWrap:
        mol = copy(self._mol)
        Chem.rdmolops.RemoveStereochemistry(mol)
        return ChemWrap(mol, self._name, self._key)

    def with_2d_coords(self, with_hydrogens: bool = True) -> ChemWrap:
        return self._with_coords(False, with_hydrogens=with_hydrogens)

    def with_3d_coords(self, with_hydrogens: bool = True) -> ChemWrap:
        return self._with_coords(True, with_hydrogens=with_hydrogens)

    def kekulize(self) -> ChemWrap:
        mol = Chem.Kekulize(self._mol)
        # don't make a copy because the inchi might be different
        new = ChemWrap(mol, self._name, self._key)
        # I don't know whether this is necessary
        new._smiles = Chem.MolToSmiles(self._mol, kekuleSmiles=True)
        return new

    def _with_coords(self, three_d: bool, with_hydrogens: bool) -> ChemWrap:
        if with_hydrogens:
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

    def desalt(self) -> ChemWrap:
        # don't make a copy because the inchi might be different
        desalted = SaltRemover.SaltRemover().StripMol(self._mol)
        return ChemWrap.of(desalted, self._name, self._key)

    def deduplicate(self) -> ChemWrap:
        parts = set(self._smiles.split("."))
        joined = ".".join(parts)
        mol = Chem.MolFromSmiles(joined)
        return ChemWrap(mol, self.name, self.key)

    def longest_component(self) -> ChemWrap:
        parts = set(self._smiles.split("."))
        if len(parts) > 1:
            logger.debug(
                "SMILES {} (name={}, key={}) contains {} nonidentical independent components".format(
                    self.smiles, self.name, self.key, len(parts)
                )
            )
        longest = ""
        for p in parts:
            if len(p) > len(longest):
                longest = p
        mol = Chem.MolFromSmiles(longest)
        return ChemWrap(mol, self.name, self.key)

    def equiv_to(self, other: ChemWrap) -> bool:
        return self.kekulize().inchikey == other.kekulize().inchikey


__all__ = ['ChemWrap']
