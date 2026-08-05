"""
Shared colour scheme for every figure in this pipeline.

Built on the LaCroix palettes, mainly Pamplemousse, so figures match the rest of
the bird paper.  Import from here rather than defining colours per script: the
point is that a locus, a gene segment or a "neutral background" looks the same
in every panel of every figure.

Deliberate choices
------------------
* No green-vs-red oppositions.  Red/green is the least colourblind-safe pairing
  there is, and it was previously used for the possible/impossible contrast.
  That contrast is now teal vs rose, which stays legible under deuteranopia and
  protanopia and is easier to look at.
* Neutral things -- background distributions, unused genes, expected-by-chance
  values -- are light grey, never a saturated colour.
* Sequential scales run from a visible tint of the category colour to the full
  colour, never from white, so low values do not vanish.
"""
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

# ─── LaCroix palettes ────────────────────────────────────────────────────────
PAMPLEMOUSSE = ["#EA7580", "#F6A1A5", "#F8CD9C", "#1BB6AF", "#088BBE", "#172869"]
CERISE_LIMON = ["#EE4244", "#F8D961", "#B6D944", "#638E6E", "#3C5541", "#132157"]
PURE = ["#AFDFEF", "#54BCD1", "#1BB6AF", "#0099D5", "#007BC3", "#172869"]
BERRY = ["#B25D91", "#CB87B4", "#EFC7E6", "#1BB6AF", "#088BBE", "#172869"]

# ─── semantic assignments used across the paper ──────────────────────────────
LOCUS = {
    "IGH": "#87b4dc",
    "IGL": "#638E6E",
    "IGK": "#F6A1A5",
    "TRA": "#EA7580",
    "TRB": "#F8CD9C",
    "TRG": "#088BBE",
}

SEGMENT = {"V": "#172869", "D": "#088BBE", "J": "#1BB6AF"}

# Possible vs impossible / matched vs mismatched: blue against rose.
# Not teal, because teal is reserved for the J segment and the two would collide
# in the arc plots where a J marker sits alongside the arcs.
YES = "#088BBE"          # donor available, matched reference, "as expected"
NO = "#EA7580"           # donor deleted, mismatched reference, "impossible"
BOTH = "#B8BCC4"         # shared between two sets

# neutral / background
GREY = "#D9D9D9"         # unused gene, empty category
GREY_DARK = "#9BA1A9"    # expected-by-chance, null distributions
INK = "#172869"          # darkest accent, axis emphasis

# two-class comparisons that are not yes/no
CLASS_A = "#088BBE"      # e.g. differences outside tracts (SHM)
CLASS_B = "#EA7580"      # e.g. differences inside tracts (conversion)


# Heatmaps use viridis REVERSED, so yellow marks low counts and dark purple
# marks high ones. Perceptually uniform and colourblind-safe in both directions.
HEATMAP = "viridis_r"


def ramp(colour, light=0.30):
    """
    Sequential colormap from a visible tint of `colour` up to `colour`.

    Starting at white makes small counts invisible, which was a real problem on
    the locus maps; `light` is how far toward white the pale end sits.
    """
    base = mcolors.to_rgb(colour)
    pale = tuple(1 - light * (1 - c) for c in base)
    return LinearSegmentedColormap.from_list("ramp", [pale, base])


def locus_ramp(locus):
    return ramp(LOCUS.get(locus, PAMPLEMOUSSE[0]))


def cycle(palette=None):
    """Categorical colours for arbitrary series."""
    return list(palette or PAMPLEMOUSSE)


def save_figure(fig, path, dpi=200, formats=("pdf", "png", "svg")):
    """
    Write one figure in every format we need, from a single output path.

    Callers pass whatever filename the workflow asked for (usually .pdf) and get
    siblings in the other formats next to it.  Keeping this in one place means a
    figure can never end up existing as a PDF but not an SVG: PDF for LaTeX, PNG
    for quick viewing and slides, SVG for editing panels in Illustrator or
    Inkscape without re-running the analysis.
    """
    import os

    stem, ext = os.path.splitext(path)
    written = []
    for fmt in formats:
        target = path if ext.lower() == "." + fmt else f"{stem}.{fmt}"
        # dpi only affects the raster output; the vector formats ignore it.
        fig.savefig(target, bbox_inches="tight", dpi=dpi, format=fmt)
        written.append(target)
    return written


class MultiPageFigures:
    """
    Drop-in replacement for PdfPages that also writes each page as PNG and SVG.

    SVG has no notion of multiple pages, so a multi-page report cannot simply be
    re-saved in that format: each page becomes its own file, <stem>_p01.svg and
    so on, while the PDF still holds the whole document.  Used for the alignment
    and rank-comparison reports, which run to dozens of pages.
    """

    def __init__(self, path, dpi=200, formats=("png", "svg")):
        import os
        from matplotlib.backends.backend_pdf import PdfPages

        self._pdf = PdfPages(path)
        self._stem = os.path.splitext(path)[0]
        self._dpi = dpi
        self._formats = formats
        self._n = 0

    def savefig(self, fig, **kw):
        self._n += 1
        self._pdf.savefig(fig, **kw)
        bbox = kw.get("bbox_inches", "tight")
        for fmt in self._formats:
            fig.savefig(f"{self._stem}_p{self._n:02d}.{fmt}",
                        bbox_inches=bbox, dpi=self._dpi, format=fmt)

    def __enter__(self):
        self._pdf.__enter__()
        return self

    def __exit__(self, *exc):
        return self._pdf.__exit__(*exc)
