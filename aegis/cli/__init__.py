import typer

from .extract import main as extract_main
from .merge import main as merge_main
from .motif_search import main as motifs_main
from .orthology import main as orthology_main
from .overlap import main as overlap_main
from .prune import main as prune_main
from .reformat import main as reformat_main
from .rename import main as rename_main
from .split import main as split_main
from .subset import main as subset_main
from .summary import main as summary_main
from .symbols import main as symbols_main
from .tidy import main as tidy_main
from .tidy_genome import main as tidy_genome_main
from . import list as list_cmd

app = typer.Typer(
    name="aegis",
    help="Annotation Extraction Genomic Integration Suite.",
    add_completion=False,
    no_args_is_help=True,
)

# Single-command modules: register their main function directly as a command
app.command(name="extract", help="Extract sequences from a genome based on an annotation file.")(extract_main)
app.command(name="merge", help="Merge two annotation files.")(merge_main)
app.command(name="motifs", help="Search for DNA motifs in promoter regions.")(motifs_main)
app.command(name="orthology", help="Pairwise orthology analysis between two annotations from different genome assembly.")(orthology_main)
app.command(name="overlap", help="Detect overlaps between annotations associated to the same genome assembly.")(overlap_main)
app.command(name="prune", help="Remove chosen gene or transcript models from an annotation.")(prune_main)
app.command(name="reformat", help="Reformat an annotation file between GFF and GTF formats.")(reformat_main)
app.command(name="rename", help="Rename gene/transcript/subfeature IDs in an annotation.")(rename_main)
app.command(name="split", short_help="Split an annotation and/or genome assembly based on genomic feature names, haplotypes, or patterns.")(split_main)
app.command(name="subset", help="Subset an annotation, and optionally its corresponding genome assembly, in various different ways. May your testing be lite.")(subset_main)
app.command(name="summary", help="Output summary statistics for an annotation.")(summary_main)
app.command(name="symbols", help="Add gene symbols to an annotation.")(symbols_main)
app.command(name="tidy", help="Clean and reformat an annotation file, allowing customisable format flavours.")(tidy_main)
app.command(name="tidy-genome", help="Clean and reformat a genome FASTA file.")(tidy_genome_main)

# Multi-command module: keep as a sub-group
app.add_typer(list_cmd.app, name="list", help="List genes or transcripts from an annotation file.")

if __name__ == "__main__":
    app()
