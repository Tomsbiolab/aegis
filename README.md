# Aegis: Annotation Extraction Genomic Integration Suite

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python Version](https://img.shields.io/badge/python-3.7%2B-brightgreen.svg)](https://www.python.org/downloads/)

**Aegis** is a powerful and flexible Python-based suite for the manipulation, analysis, and integration of genomic annotations. It provides a robust, object-oriented framework for working with genomic data, enabling complex analyses and data transformations with intuitive, high-level commands.

## Key Features

- **Object-Oriented Design**: Aegis represents genomic features (genes, transcripts, exons, etc.) as a hierarchical system of Python classes, providing a clean and intuitive API for data manipulation.
- **Comprehensive Annotation Handling**: Seamlessly parse, process, and export genomic annotations in GFF3 format.
- **Sequence Extraction**: Easily extract sequences for any genomic feature, including genes, transcripts, CDS, proteins, and promoters.
- **Comparative Genomics**: Perform comparative analyses, such as identifying orthologs and syntenic regions between different species.
- **Extensible and Modular**: The modular design of Aegis allows for easy extension and integration with other bioinformatics tools and pipelines.

## The Aegis Class System

The core of Aegis is its sophisticated class system, which models the hierarchical nature of genomic annotations. This object-oriented approach provides several key advantages over traditional, line-by-line processing of annotation files:

- **Intuitive Data Representation**: Genomic features are not just lines in a file; they are objects with properties and relationships. A `Gene` object contains `Transcript` objects, which in turn contain `Exon` and `CDS` objects. This makes the code more readable, maintainable, and less error-prone.
- **Data Integrity**: The class system enforces data consistency. For example, when a `Gene` object is updated, all its associated `Transcript` and sub-feature objects are updated accordingly, ensuring that the annotation remains coherent.
- **Complex Queries and Manipulations**: The object-oriented structure allows for complex queries and manipulations that would be difficult to perform with traditional text-based tools. For example, you can easily retrieve all coding transcripts for a specific gene, or calculate the total length of all exons in a given transcript.
- **Code Reusability**: The class-based design promotes code reusability. Once you have defined a class for a specific genomic feature, you can reuse it in different parts of your analysis pipeline.

### Core Classes

- **`Genome`**: Represents a genome, containing a collection of `Scaffold` objects.
- **`Scaffold`**: Represents a chromosome or scaffold, containing the sequence and a collection of `Gene` objects.
- **`Annotation`**: The main container for genomic annotations, holding a collection of `Gene` objects.
- **`Gene`**: Represents a gene, containing one or more `Transcript` objects.
- **`Transcript`**: Represents a transcript, containing `Exon`, `CDS`, and `UTR` objects.
- **`Exon`**, **`CDS`**, **`UTR`**, **`Intron`**: Represent the sub-features of a transcript.
- **`Protein`**, **`Promoter`**: Represent other biological features of interest.

## Installation

To install Aegis, clone the repository and install the required dependencies:

```bash
git clone https://github.com/your-username/aegis.git
cd aegis
pip install -r requirements.txt
```

## Usage

Aegis is designed to be used as a library in your Python scripts. Here is a simple example of how to load an annotation and extract the sequences of all genes:

```python
from aegis.annotation import Annotation
from aegis.genome import Genome

# Load the genome and annotation
genome = Genome("my_genome", "path/to/genome.fasta")
annotation = Annotation("my_annotation", "path/to/annotation.gff3", genome)

# Generate and export gene sequences
annotation.generate_sequences(genome)
annotation.export_genes()
```

## Contributing

Contributions to Aegis are welcome! Please feel free to submit a pull request or open an issue on the GitHub repository.

## License

Aegis is licensed under the GNU General Public License v3.0. See the `LICENSE.md` file for more details.
