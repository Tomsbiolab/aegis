#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 15:26:49 2023

@authors: David Navarro, Antonio Santiago
"""

import subprocess
import shlex
import pickle
import os

def pickle_load(file):
    f = open(file, "rb")
    item = pickle.load(f)
    f.close()
        
    return item


def pickle_save(file, item):
    f = open(file, "wb")
    pickle.dump(item, f)
    f.close()

    
def bash_run(command, tag="", folder_out="", extension:str=".txt", internal_output=False, standard_output=True, error_output=True):
    proc = subprocess.run(shlex.split(command), capture_output=True, encoding="utf-8")

    if standard_output:
        if proc.stdout != "":
            f = open(folder_out + tag + extension, "w", encoding="utf-8")
            f.write(proc.stdout)
            f.close()
            
    if error_output:
        error_folder = folder_out + "error_files/"
        if proc.stderr != "":
            os.system("mkdir -p " + error_folder)
            f = open(error_folder + tag + "_stderr.txt", "w", encoding="utf-8")
            f.write(proc.stderr)
            f.close()
    
    if internal_output:
        return proc
    