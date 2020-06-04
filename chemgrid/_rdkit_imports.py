import logging

logger = logging.getLogger("chemgrid")


try:
    import cairosvg
except (ImportError, OSError):
    # if the backend is not installed, throws an OSError
    logger.debug("Failed importing cairosvg.", exc_info=True)
    logger.info("Did not find cairosvg. Assuming 'client' mode.")
    cairosvg = None

# if it's raw PIL instead of Pillow, it won't import in Python 3
try:
    from PIL.Image import Image
except ImportError:
    logger.debug("Failed importing pillow.", exc_info=True)
    logger.info("Did not find pillow. Assuming 'client' mode.")
    Image = None

try:
    from rdkit import Chem
    from rdkit.Chem import SaltRemover
    from rdkit.Chem import Mol
    from rdkit.Chem import Draw
    import rdkit.Chem.inchi as Inchi
except ImportError:
    logger.debug("Failed importing rdkit.", exc_info=True)
    logger.info("Did not find rdkit. Assuming 'client' mode.")
    Chem = None
    Draw = None
    Mol = None
    Inchi = None
    SaltRemover = None


class MoleculeError(Exception):
    pass


class MoleculeConversionError(MoleculeError):
    pass


class NullMoleculeError(MoleculeConversionError):
    pass
