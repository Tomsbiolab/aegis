
from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..annotation import Annotation

from pathlib import Path
import warnings

# 1. Create a Base Class to hold shared state and methods
class AnnotationComponent:
    _annot: Annotation
    def __init__(self, annotation:Annotation):
        self._annot = annotation

    def _resolve_output_path(self, filepath: str | None, output_dir: str | None, filename: str | None, 
        feature_type: str = "", suffix:str = "", extra_suffixes: list[str] | None = None, 
        extension: str = ".fasta", use_annot_dir: bool = False, 
        subfolder_name: str = "features", subfolder: bool = False
    ) -> Path:
        
        # Exact same logic as your original code
        if filepath:
            p = Path(filepath)
            p.parent.mkdir(parents=True, exist_ok=True)
            if (output_dir is not None) or subfolder or filename or use_annot_dir:
                warnings.warn(f"Exact output file path ({filepath}) was specified. Ignoring output_dir, use_annot_dir, subfolder, and filename arguments.")
            return p
        else:
            if use_annot_dir:
                export_folder = Path(self._annot.path)
                if output_dir is not None:
                    warnings.warn(f"Both 'use_annot_dir={use_annot_dir}' and 'output_dir={output_dir}' were provided. Defaulting to the annotation directory ({self._annot.path}).")
            else:
                export_folder = Path(output_dir or ".")

            if subfolder:
                export_folder = export_folder / subfolder_name.strip("/")

            export_folder.mkdir(parents=True, exist_ok=True)

            if not filename:
                if feature_type:
                    feature_type = f"_{feature_type}"
                tag_str = "".join([f"_{t}" for t in (extra_suffixes or []) if t])
                filename = f"{self._annot.id}{suffix}{feature_type}{tag_str}{extension}"
            
            if not Path(filename).suffix:
                filename = f"{filename}{extension}"

            return export_folder / filename