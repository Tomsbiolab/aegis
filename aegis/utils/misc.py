"""
Created on Thu Jan 19 15:26:49 2023

@authors: David Navarro, Antonio Santiago
"""

import pickle
import re
import subprocess

from pathlib import Path
from collections import Counter

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


def read_file_with_fallback(file_path, encodings=['utf-8', 'latin-1', 'ascii']):
    """
    Tries several encodings to find the suitable one.
    """
    for enc in encodings:
        try:
            with open(file_path, 'r', encoding=enc) as f:
                f.readlines()
                return enc
        except UnicodeDecodeError:
            continue

    raise ValueError(f"Not able to decodify '{file_path}'")