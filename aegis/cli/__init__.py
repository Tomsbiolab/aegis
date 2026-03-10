import typer

from . import extract
from . import list as list_cmd
from . import merge
from . import motif_search
from . import orthology
from . import overlap
from . import prune
from . import reformat
from . import rename
from . import subset
from . import summary
from . import symbols
from . import tidy
from . import tidy_genome

app = typer.Typer(
    name="aegis",
    help="Annotation Extraction Genomic Integration Suite.",
    add_completion=False,
    no_args_is_help=True,
)

app.add_typer(extract.app, name="extract", help="Extract sequences from a genome based on an annotation file.")
app.add_typer(list_cmd.app, name="list", help="List genes or transcripts from an annotation file.")
app.add_typer(merge.app, name="merge", help="Merge two annotation files.")
app.add_typer(motif_search.app, name="motifs", help="Search for DNA motifs in promoter regions.")
app.add_typer(orthology.app, name="orthology", help="Pairwise orthology analysis between two annotations from different genome assembly.")
app.add_typer(overlap.app, name="overlap", help="Detect overlaps between annotations associated to the same genome assembly.")
app.add_typer(prune.app, name="prune", help="Remove chosen gene or transcript models from an annotation.")
app.add_typer(reformat.app, name="reformat", help="Reformat an annotation file between GFF and GTF formats.")
app.add_typer(rename.app, name="rename", help="Rename gene/transcript/subfeature IDs in an annotation.")
app.add_typer(subset.app, name="subset", help="Subset an annotation, and optionally its corresponding genome assembly, in various different ways. May your testing be lite.")
app.add_typer(summary.app, name="summary", help="Output summary statistics for an annotation.")
app.add_typer(symbols.app, name="symbols", help="Add gene symbols to an annotation.")
app.add_typer(tidy.app, name="tidy", help="Clean and reformat an annotation file, allowing customisable format flavours.")
app.add_typer(tidy_genome.app, name="tidy-genome", help="Clean and reformat a genome FASTA file.")

if __name__ == "__main__":
    app()
