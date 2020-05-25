import cairosvg
from PIL.Image import Image
from rdkit import Chem
from rdkit.Chem import Draw
from chemgrid.chem import *

# TODO this is unused
opts = Draw.DrawingOptions()
opts.bgColor = None
opts.bondLineWidth = 0.6
opts.dblBondLengthFrac = 0.18
opts.atomLabelFontSize = 8


class ChemGraphicsKit:
    def __init__(self, simplifier: Optional[Callable[[str], str]] = None):
        if simplifier is None:
            simplifier = lambda s: s
        self.simplifier = simplifier

    def save_svgs(
        self,
        compounds: Iterable[CompoundLike],
        labels: Union[None, Sequence[str], CompoundNamer] = None,
        directory: PLike = "compound-svgs",
    ) -> None:
        compounds = InternalTools.fetch_all_ids_unchecked(Compounds, compounds)
        directory = Tools.prepped_dir(directory)
        if labels is None:
            labels = [str(c) for c in compounds]
        if isinstance(labels, CompoundNamer):
            dct = labels.fetch(compounds)
            labels = [dct[c] for c in compounds]
        for compound, label in Tools.zip_strict(compounds, labels):
            try:
                self.save_pretty(compound, directory / (label + ".svg"))
            except MoleculeConversionError:
                logger.exception("Failed to convert c{} / {}".format(compound, label))

    def save_pretty(self, compound: CompoundLike, path: PLike, scale: float = 1.0) -> None:
        mol = self.simplifier(ChemTools.convert(compound))
        path = Path(path)
        if path.suffix == ".svg":
            Draw.MolToFile(mol, str(path), options=opts)
        else:
            save_fn = Tools.or_raise(
                {".png": cairosvg.svg2png, ".pdf": cairosvg.svg2pdf}.get(path.suffix),
                IllegalFilenameError(str(path)),
            )
            with Tools.tmppath() as tmp_path:
                Draw.MolToFile(mol, str(tmp_path), options=opts)
                save_fn(url=str(tmp_path), write_to=str(path), scale=scale)

    def draw_single(self, compound: MolLike, label: str):
        mol = ChemTools.convert(compound)
        return Draw.MolsToGridImage([mol], legends=[label], molsPerRow=1)

    def draw_grid(
        self,
        compounds: Iterable[Union[MolLike, Iterable[MolLike]]],
        legends: Union[None, Iterable[str], Iterable[Iterable[str]]] = None,
        n_cols: int = 12,
    ):
        compounds = list(compounds)
        if legends is None:
            legends = [str(i) for i in range(len(compounds))]
        if len(compounds) != len(legends):
            raise LengthMismatchError(
                "{} SMILES but {} legends".format(len(compounds), len(legends))
            )
        mols = []
        the_legends = []
        for compound, legend in Tools.zip_list(compounds, legends):
            smiles = []
            for c in InternalTools.flatten_smart(compound):
                smiles.append(Chem.MolToSmiles(self.simplifier(ChemTools.convert(c))))
            mols.append(Chem.MolFromSmiles(".".join(smiles)))
            the_legends.append(", ".join(InternalTools.flatten_smart(legend)))
        return Draw.MolsToGridImage(mols, legends=[str(k) for k in the_legends], molsPerRow=n_cols)


__all__ = ["ChemSimplifer", "ChemGraphicsKit"]
