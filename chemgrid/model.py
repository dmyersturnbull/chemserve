from rdkit.Chem import SaltRemover, MolToSmiles
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.inchi import MolToInchiAndAuxInfo, InchiToInchiKey

MolLike = Union[Chem.Mol, str]


class MoleculeError(Exception):
    pass


class MoleculeConversionError(MoleculeError):
    pass


class NullMoleculeError(MoleculeConversionError):
    pass


class ChemTools:

    # TODO this is badly named
    @classmethod
    def convert(cls, compound: MolLike) -> Chem.Mol:
        mol = cls.try_convert(compound)
        if mol is None:
            raise MoleculeConversionError(
                "{} was not converted correctly by SMILES or InChI".format(compound)
            )
        return mol

    @classmethod
    def try_convert(cls, compound: MolLike) -> Optional[Chem.Mol]:
        if isinstance(compound, Chem.Mol):
            return compound
        if isinstance(compound, str) and compound.startswith("InChI="):
            return cls._try_convert_single(compound, Chem.MolFromInchi)
        elif isinstance(compound, str):
            return Chem.MolFromSmiles(compound)
        else:
            raise TypeError(str(type(compound)))

    @classmethod
    def _try_convert_single(
        cls, smiles_or_inchi: str, convert_fn: Callable[[str], Chem.Mol], log_fn=logger.error
    ):
        try:
            mol = convert_fn(smiles_or_inchi)
            if mol is None:
                raise NullMoleculeError("Created a null molecule from {}".format(mol))
            else:
                return mol
        except:
            log_fn("Failed to convert {} to mol".format(smiles_or_inchi), exc_info=True)


class WrappedInchi:
    def __init__(self, compound: MolLike):
        # make rdkit_inchi etc. first, then set inchi etc., then correct inchi etc. if passed and different than rdkit's result
        self.mol = ChemTools.convert(compound)
        self.rdkit_inchi, self.aux = MolToInchiAndAuxInfo(self.mol)
        self.rdkit_inchikey = InchiToInchiKey(self.rdkit_inchi)
        self.inchi, self.inchikey = self.rdkit_inchi, self.rdkit_inchikey
        self.rdkit_smiles = MolToSmiles(self.mol)
        self.smiles = self.rdkit_smiles
        if (
            isinstance(compound, str)
            and compound.startswith("InChI=")
            and self.rdkit_inchi != compound
        ):
            self.inchi, self.inchikey = compound, InchiToInchiKey(compound)
        elif isinstance(compound, str) and self.rdkit_smiles != compound:
            self.smiles = compound

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__, self.inchikey)

    def __str__(self):
        return repr(self)

    def __eq__(self, other):
        return (
            isinstance(other, self.__class__)
            and self.inchi == other.inchi
        )


class ChemSimplifer:
    """
    De-salts, removes duplicate substructures, and potentially does more.
    """

    def __init__(self, aggressive: bool = False):
        """
        :param aggressive: If there are > 1 nonidentical components, takes the longest and logs with level 'caution'; otherwise, only deduplciates
        """
        self.aggressive = aggressive

    def desalt(self, mol: MolLike) -> Chem.Mol:
        mol = ChemTools.convert(mol)
        return SaltRemover.SaltRemover().StripMol(mol)

    def simplify(self, compound: MolLike) -> Chem.Mol:
        desalted = self.desalt(compound)
        parts = set(Chem.MolToSmiles(desalted).split("."))
        if self.aggressive:
            if len(parts) > 1:
                logger.caution(
                    "SMILES {} contains {} nonidentical independent components".format(
                        desalted, len(parts)
                    )
                )
            longest = Tools.longest_str(parts)
        else:
            # just keep the unique parts
            longest = ".".join(parts)
        return Chem.MolFromSmiles(longest)


__all__ = [
    "ChemTools",
    "ChemSimplifer",
    "ChemblApi",
    "WrappedInchi",
    "MolLike",
    "MoleculeError",
    "MoleculeConversionError",
    "NullMoleculeError",
]
