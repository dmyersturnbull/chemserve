import logging
from pathlib import Path
from typing import Optional, Union, Iterable, Sequence, Callable
from chemserve.models import Concrete

from chemserve._rdkit_imports import Draw, cairosvg, Image

logger = logging.getLogger("chemserve")
PathLike = Union[Path, str]


class DrawingOptionSets:
    @classmethod
    def nature_springer(cls):
        opts = Draw.DrawingOptions()
        opts.bgColor = None
        opts.bondLineWidth = 0.6
        opts.dblBondLengthFrac = 0.18
        opts.atomLabelFontSize = 8
        return opts

    def __repr__(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__


class ChemGraphicsKit:
    def __init__(self, opts=None, scale: float = 1.0):
        self._opts = opts
        self._scale = scale

    def save_svgs(self, chems: Iterable[Concrete], directory: PathLike) -> None:
        directory = Path(directory)
        directory.mkdir(exist_ok=True, parents=True)
        for chem in chems:
            label = self._label(chem)
            self.save_pretty(chem, directory / (label + ".svg"))

    def save_pretty(self, chem: Concrete, path: PathLike) -> None:
        path = Path(path)
        accepted = {".svg": None, ".png": cairosvg.svg2png, ".pdf": cairosvg.svg2pdf}
        mol = chem.mol
        save_fn = accepted.get(path.suffix)
        if save_fn is None:
            Draw.MolToFile(mol, str(path), options=self._opts)
        else:
            tmppath = path.with_suffix(path.suffix + ".tmp")
            Draw.MolToFile(mol, str(tmppath), options=self._opts)
            save_fn(url=str(tmppath), write_to=str(path), scale=self._scale)

    def draw_grid(
        self,
        chems: Iterable[Union[Concrete, Iterable[Concrete]]],
        n_rows: int,
        n_cols: int,
        opts=None,
    ) -> Image:
        wrapped = self._wrap_up(chems)
        mols = [chem.mol for chem in wrapped]
        labels = [self._label(chem) for chem in wrapped]
        return Draw.MolsToGridImage(mols, legends=labels, molsPerRow=n_cols, opts=opts)

    def draw_single(self, chem: Concrete) -> Image:
        label = self._label(chem)
        return Draw.MolsToGridImage([chem.mol], legends=[label], molsPerRow=1)

    def _wrap_up(self, chems: Iterable[Union[Concrete, Iterable[Concrete]]]) -> Sequence[Concrete]:
        """
        Wraps up potentially nested list by ``ChemWrap.join``.
        """
        wrapped = []
        for chem in chems:
            chem = Concrete.join([chem] if isinstance(chem, Concrete) else chem)
            wrapped.append(chem)
        return wrapped

    def _label(self, chem: Concrete) -> str:
        if chem.name is not None:
            return chem.name
        elif chem.key is not None:
            return chem.key
        else:
            return chem.inchikey

    def __repr__(self):
        return "{}@{}".format(self.__class__.__name__, hex(id(self)))

    def __str__(self):
        return "{}@{}".format(self.__class__.__name__, hex(id(self)))


__all__ = ["ChemGraphicsKit", "DrawingOptionSets"]
