"""
Utility functions and callbacks for the AEGIS CLI suite.
"""
from typing import List, Sequence, Union, Optional

import typer
from typing_extensions import Annotated

from ..utils.genefunctions import NCBI_GENETIC_CODES, TAXONOMY_ORGANELLE_CODES


def split_callback(value: Union[str, Sequence[str], None]) -> List[str]:
    """
    Callback for Typer/Click options that accept comma-separated strings
    (e.g., '--features gene,CDS') or repeated flags (e.g., '-f gene -f CDS').

    Always returns a clean list of trimmed strings.
    """
    if not value:
        return []
    if isinstance(value, (list, tuple, set)):
        result: List[str] = []
        for item in value:
            if isinstance(item, str):
                result.extend([x.strip() for x in item.split(",") if x.strip()])
            elif item is not None:
                result.append(str(item).strip())
        return result
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    return [str(value).strip()]


# ---------------------------------------------------------------------------
# Genetic code options shared by every command that translates CDSs
# ---------------------------------------------------------------------------
GENETIC_CODES_PANEL = "Genetic Codes"


def taxonomy_callback(value: str) -> str:
    if value.lower() not in TAXONOMY_ORGANELLE_CODES:
        raise typer.BadParameter(f"Unknown taxonomy '{value}'. Choose from: {', '.join(TAXONOMY_ORGANELLE_CODES)}.")
    return value.lower()


def genetic_code_callback(value: Optional[int]) -> Optional[int]:
    if value is not None and value not in NCBI_GENETIC_CODES:
        raise typer.BadParameter(f"Unsupported NCBI genetic code {value}. Supported tables: {', '.join(str(t) for t in NCBI_GENETIC_CODES)}.")
    return value


TaxonomyOption = Annotated[str, typer.Option(
    "-tax", "--taxonomy", callback=taxonomy_callback, rich_help_panel=GENETIC_CODES_PANEL,
    help="Taxonomic group preset for organelle genetic codes: 'plant' (default: nuclear=1, mito=1, plastid=11), 'vertebrate' (nuclear=1, mito=2), 'invertebrate' (nuclear=1, mito=5), or 'yeast' (nuclear=1, mito=3). Specific codes can be individually customized with --genetic-code, --mito-code, or --plastid-code."
)]
GeneticCodeOption = Annotated[int, typer.Option(
    "-gc", "--genetic-code", callback=genetic_code_callback, rich_help_panel=GENETIC_CODES_PANEL,
    help="NCBI genetic code table number for nuclear translation (default: 1)."
)]
MitoCodeOption = Annotated[Optional[int], typer.Option(
    "-mc", "--mito-code", callback=genetic_code_callback, rich_help_panel=GENETIC_CODES_PANEL,
    help="NCBI genetic code table number for mitochondrial contigs (overrides --taxonomy default)."
)]
PlastidCodeOption = Annotated[Optional[int], typer.Option(
    "-pc", "--plastid-code", callback=genetic_code_callback, rich_help_panel=GENETIC_CODES_PANEL,
    help="NCBI genetic code table number for plastid/chloroplast contigs (overrides --taxonomy default)."
)]
AutoOrganelleCodesOption = Annotated[bool, typer.Option(
    "--auto-organelle-codes/--no-auto-organelle-codes", rich_help_panel=GENETIC_CODES_PANEL,
    help="Automatically use mitochondrial and plastid genetic codes for contigs detected as organelles from their name or FASTA description (default: True)."
)]
MitochondriaChromsOption = Annotated[List[str], typer.Option(
    "--mitochondria-chroms", callback=split_callback, rich_help_panel=GENETIC_CODES_PANEL,
    help="Explicit list or comma-separated names of mitochondrial chromosomes/scaffolds to translate with --mito-code."
)]
ChloroplastChromsOption = Annotated[List[str], typer.Option(
    "--chloroplast-chroms", callback=split_callback, rich_help_panel=GENETIC_CODES_PANEL,
    help="Explicit list or comma-separated names of chloroplast/plastid chromosomes/scaffolds to translate with --plastid-code."
)]


def initiator_methionine_callback(value: str) -> str:
    if value not in ("canonical", "all", "none"):
        raise typer.BadParameter(f"Invalid value '{value}'. Choose from: canonical, all, none.")
    return value


InitiatorMethionineOption = Annotated[str, typer.Option(
    "--initiator-methionine", callback=initiator_methionine_callback, rich_help_panel=GENETIC_CODES_PANEL,
    help="Which start codons are translated as M when they initiate a 5'-complete protein: 'canonical' (default: the curated start codons of the table, e.g. GTG/TTG in plastids), 'all' (also NCBI's alternative initiators, e.g. CTG/TTG in table 1), or 'none' (literal translation of the first codon)."
)]
