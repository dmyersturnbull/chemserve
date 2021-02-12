from __future__ import annotations
from typing import Iterator, List


class Fingerprint:
    """
    Just a simple wrapper for rdkit fingerprints.
    A bit unnecessary, but convenient when you're using them a lot.
    """

    def __init__(self, fp):
        self._fp = fp

    ###
    # NOTE: Often there is both a magic function (like __str__)
    # and a property (for a more fluent style; like .string)
    # These have duplicate code because delegating causes a small performance hit
    # This class is small, so that's easy to maintain -- just copy and paste
    ###

    @property
    def bytes(self) -> bytes:
        return self._fp.ToBinary()

    @property
    def numpy(self):
        import numpy as np

        # NOTE: frombuffer will NOT work correctly for bool arrays
        # also, fromiter is much slower than creating a list first
        # This is appears to be the fastest way to create an array here
        return np.array(list(self._fp), dtype=bool)

    @property
    def list(self) -> List[bool]:
        return list(map(bool, self._fp))

    @property
    def string(self) -> str:
        return self._fp.ToBitString()

    @property
    def base64(self) -> str:
        return self._fp.ToBase64()

    @property
    def n_bits(self) -> int:
        return self._fp.GetNumBits()

    @property
    def n_on(self) -> int:
        return self._fp.GetNumOnBits()

    @property
    def n_off(self) -> int:
        return self._fp.GetNumOffBits()

    # TODO: Consider changing to hold bytes or numpy array, and implement | and &
    # def __ror__(self, other: Fingerprint) -> Fingerprint:
    # https://bugs.python.org/issue19251

    def __len__(self) -> int:
        return self._fp.GetNumBits()

    def __str__(self) -> str:
        return self._fp.ToBitString()

    def __repr__(self) -> str:
        return self._fp.ToBitString()

    def __bytes__(self) -> bytes:
        return self._fp.ToBinary()

    def __iter__(self) -> Iterator[bool]:
        return iter(map(bool, self._fp))


__all__ = ["Fingerprint"]
