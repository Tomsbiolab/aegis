"""
Created on Thu Jan 19 15:26:49 2023

@authors: David Navarro, Antonio Santiago
"""
from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any, TextIO

import pickle
import re
import subprocess
import sys
import gzip

import pandas as pd

from pathlib import Path
from collections import Counter
from tqdm import tqdm

def pickle_load(file):
    f = open(file, "rb")
    item = pickle.load(f)
    f.close()
        
    return item

def pickle_save(file, item):
    f = open(file, "wb")
    pickle.dump(item, f)
    f.close()

def covered_length(intervals) -> int:
    """
    Number of bases covered by 1-based inclusive (start, end) intervals, each base counted once however
    many intervals overlap it.
    """
    total = 0
    cur_start = cur_end = None
    for start, end in sorted(intervals):
        if cur_end is None or start > cur_end:
            if cur_end is not None:
                total += cur_end - cur_start + 1
            cur_start, cur_end = start, end
        else:
            cur_end = max(cur_end, end)
    if cur_end is not None:
        total += cur_end - cur_start + 1
    return total

def shared_bases(features1, features2, same_frame:bool=False) -> int:
    """
    Bases covered by both lists of features, each counted once, so it matches the covered_length
    sizes of exons and CDSs even when rows of one list overlap (e.g. a -1 ribosomal slippage).
    same_frame only counts the bases of pairs of rows read in the same frame.
    """
    shared = []
    for f1 in features1:
        for f2 in features2:
            if same_frame and f1.frame != f2.frame:
                continue
            start = max(f1.start, f2.start)
            end = min(f1.end, f2.end)
            if start <= end:
                shared.append((start, end))
    return covered_length(shared)

def unambiguous_numbers(numbered) -> bool:
    """
    Whether (key, number) pairs can keep their numbers when renaming: every key has a number
    and no number is shared by different keys (features with the same key, e.g. an exon shared
    by transcripts, may share it).
    """
    owners = {}
    for key, number in numbered:
        if number is None or owners.setdefault(number, key) != key:
            return False
    return True

def count_occurrences(string, char):
    return Counter(string)[char]

def find_all_occurrences(pattern, text, flags: int = 0):
    matches = []
    for match in re.finditer(pattern, text, flags):
        matches.append((match.start(), match.end(), match.group()))

    return matches

def start_progress_bar(total: int, description: str, quiet: bool = False, colour: str = "92"):
    disable = quiet or not sys.stderr.isatty()

    return tqdm(total=total, disable=disable, bar_format=(
            f'\033[1;{colour}m{description}:\033[0m '
            '{percentage:3.0f}%|'
            f'\033[1;{colour}m{{bar}}\033[0m| '
            '{n}/{total} [{elapsed}<{remaining}]'))

def read_headerless_tsv(file, columns:list, **kwargs) -> pd.DataFrame:
    """
    Reads a tab separated table without a header (e.g. DIAMOND or MCscan output) and names its
    columns. A file without rows (e.g. a search with no hits) gives an empty table with those columns.
    """
    try:
        df = pd.read_csv(file, sep="\t", header=None, **kwargs)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=columns)
    df.columns = columns
    return df

def run_command(working_directory: Path, command: list):
    """
    Executes an external command (e.g. Liftoff, DIAMOND) in working_directory.

    Args:
        working_directory (Path): The working directory for the command.
        command (list): The command and its arguments as a list of strings.

    Raises:
        subprocess.CalledProcessError: If the command fails.
    """
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True, cwd=working_directory)
        return result
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {' '.join(command)}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        raise


def open_file(file_path:Any, mode:str='r', encoding:str|None=None) -> TextIO:
    """
    Transparently opens a text file or a gzipped text file, depending on the extension.
    """
    if str(file_path).endswith('.gz'):
        if 'r' in mode and 'b' not in mode:
            mode += 't'
        return gzip.open(file_path, mode, encoding=encoding) # type: ignore
    return open(file_path, mode, encoding=encoding) # type: ignore


def read_file_with_fallback(file_path, encodings=('utf-8', 'cp1252', 'latin-1'), chunk_size=1 << 20):
    """
    Returns the first encoding that decodes the whole file (read in chunks), as a single byte that is not
    valid in an encoding can be anywhere (e.g. a cp1252 en dash, 0x96, deep into an otherwise utf-8 GFF).
    cp1252 (Windows text) is tried before latin-1, which decodes any byte.
    """
    for enc in encodings:
        try:
            with open_file(file_path, 'r', encoding=enc) as f:
                while f.read(chunk_size):
                    pass
                return enc
        except UnicodeDecodeError:
            continue

    raise ValueError(f"Not able to decodify '{file_path}'")