#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import subprocess, yaml

def run(command:str, shell:bool = True, executable:str = '/bin/bash', Popen:bool = False):
    """This function is just a useful wrapper around subprocess.Popen, subprocess.run

    Args:
        command (str): Any command to execute on Linux
        shell (bool, optional): keyword of Popen and Run. Defaults to True.
        executable (str, optional): keyword of Popen and Run. Defaults to '/bin/bash'.
        Popen (bool, optional): If True it will launch popen if not Run. Defaults to False.

    Raises:
        RuntimeError: In case of non-zero exit status.

    Returns:
        object: the processes returned by Popen or Run.
    """
    if Popen:
        #In this case you could acces the pid as: run.pid
        process = subprocess.Popen(command, shell = shell, executable = executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
    else:
        process = subprocess.run(command, shell = shell, executable = executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
        returncode = process.returncode
        if returncode != 0:
            print(f'Command {command} returned non-zero exit status {returncode}')
            raise RuntimeError(process.stderr)
    return process

