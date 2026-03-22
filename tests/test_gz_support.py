import gzip
import shutil
from pathlib import Path

from aegis.annotation import Annotation
from aegis.genome import Genome
from .conftest import TEST_DATA_DIR

def test_gff3_gz_support(tmp_path):
    gff3_path = TEST_DATA_DIR / "input/annotation/minimal.gff3"
    
    gz_path = tmp_path / "minimal.gff3.gz"
    with open(gff3_path, 'rb') as f_in:
        with gzip.open(gz_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out) # type: ignore
            
    # Load Annotation from the .gz file
    annot = Annotation(str(gz_path), quiet=True)
    
    assert "gene" in annot.features
    assert annot.features["gene"] >= 1
    assert "mRNA" in annot.features
    assert annot.features["mRNA"] >= 1

def test_gtf_gz_support(tmp_path):
    gtf_path = TEST_DATA_DIR / "input/annotation/convert_basic.gtf"
    
    gz_path = tmp_path / "convert_basic.gtf.gz"
    with open(gtf_path, 'rb') as f_in:
        with gzip.open(gz_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out) # type: ignore
            
    # Load Annotation from the .gtf.gz file
    annot = Annotation(str(gz_path), quiet=True)
    
    assert "gene" in annot.features
    assert annot.features["gene"] >= 1
    assert "mRNA" in annot.features or "transcript" in annot.features

def test_fasta_gz_support(tmp_path):
    fasta_path = TEST_DATA_DIR / "input/fasta/minimal.fasta"
    
    gz_path = tmp_path / "minimal.fasta.gz"
    with open(fasta_path, 'rb') as f_in:
        with gzip.open(gz_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out) # type: ignore
            
    # Load Genome from the .fasta.gz file
    genome = Genome("minimal", str(gz_path), quiet=True)
    
    assert len(genome.scaffolds) > 0
