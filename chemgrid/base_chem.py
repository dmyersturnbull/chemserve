import abc
from functools import total_ordering

@total_ordering
class BaseChem(metaclass=abc.ABCMeta):

    @property
    def inchi_or_smiles(self) -> str:
        """
        Returns:
            Either InChI or SMILES
        """
        raise NotImplementedError()

    @property
    def name(self) -> str:
        raise NotImplementedError()

    @property
    def key(self) -> str:
        raise NotImplementedError()

    def __eq__(self, other):
        # TODO too harsh
        if not isinstance(other, self.__class__):
            raise TypeError("Universal equality unsupported")
        return (self.inchi_or_smiles, self.name, self.key) == (other.inchi_or_smiles, other.name, other.key)

    def __hash__(self):
        return hash((self.inchi_or_smiles, self.name, self.key))

    def __lt__(self, other):
        # TODO too harsh
        if not isinstance(other, self.__class__):
            raise TypeError("Universal equality unsupported")
        # the most useful for display, presumably
        return (self.name, self.key, self.inchi_or_smiles) < (other.name, other.key, other.inchi_or_smiles)


__all__ = ['BaseChem']
