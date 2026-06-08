
from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..annotation import Annotation

# 1. Create a Base Class to hold shared state and methods
class AnnotationComponent:
    _annot: Annotation
    def __init__(self, annotation:Annotation):
        self._annot = annotation

    def _resolve_output_path(self, *args, **kwargs):
        return self._annot._resolve_output_path(*args, **kwargs)

    def _resolve_output_dir(self, *args, **kwargs):
        return self._annot._resolve_output_dir(*args, **kwargs)