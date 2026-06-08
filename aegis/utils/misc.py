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

def count_occurrences(string, char):
    return Counter(string)[char]

def find_all_occurrences(pattern, text):
    matches = []
    for match in re.finditer(pattern, text):
        matches.append((match.start(), match.end(), match.group()))

    return matches

def start_progress_bar(total: int, description: str, quiet: bool = False, colour: str = "92"):
    disable = quiet or not sys.stderr.isatty()

    return tqdm(total=total, disable=disable, bar_format=(
            f'\033[1;{colour}m{description}:\033[0m '
            '{percentage:3.0f}%|'
            f'\033[1;{colour}m{{bar}}\033[0m| '
            '{n}/{total} [{elapsed}<{remaining}]'))

def run_command(working_directory: Path, command: list):
    """
    Executes a generic command inside a Docker container.

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


def read_file_with_fallback(file_path, encodings=['utf-8', 'ascii', 'latin-1'], sample_size=100000):
    """
    Tries several encodings to find the suitable one.
    """
    for enc in encodings:
        try:
            with open_file(file_path, 'r', encoding=enc) as f:
                f.read(sample_size)
                return enc
        except UnicodeDecodeError:
            continue

    raise ValueError(f"Not able to decodify '{file_path}'")