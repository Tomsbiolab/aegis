"""
Utility functions and callbacks for the AEGIS CLI suite.
"""
from typing import List, Sequence, Union, Optional


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
