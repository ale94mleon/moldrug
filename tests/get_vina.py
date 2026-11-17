"""Download the AutoDock-Vina executable used by the test suite.

AutoDock-Vina is the only non-pip dependency of moldrug. Upstream publishes a
binary per operating system *and* architecture, so both have to be detected:
macOS runners on CI are arm64, and an x86_64 binary does not run there.
"""
import os
import platform
import stat
import subprocess  # nosec
import sys
import urllib.request

VINA_VERSION = '1.2.5'
BASE_URL = 'https://github.com/ccsb-scripps/AutoDock-Vina/releases/download'


def run(command):
    subprocess.run(command, shell=True, executable='/bin/bash',  # nosec
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def get_asset_name() -> str:
    """Name of the AutoDock-Vina asset for the running platform.

    Returns
    -------
    str
        E.g. ``vina_1.2.5_linux_x86_64``.

    Raises
    ------
    NotImplementedError
        For an operating system or architecture upstream does not build for.
    """
    system = platform.system().lower()
    if system == 'darwin':
        operating_system = 'mac'
    elif system == 'linux':
        operating_system = 'linux'
    else:
        raise NotImplementedError(
            f"Only MacOS and Linux are supported at the moment, got {platform.system()!r}. "
            "Install AutoDock-Vina yourself and point the tests at it.")

    machine = platform.machine().lower()
    if machine in ('x86_64', 'amd64'):
        architecture = 'x86_64'
    elif machine in ('aarch64', 'arm64'):
        architecture = 'aarch64'
    else:
        raise NotImplementedError(
            f"No AutoDock-Vina build for the {machine!r} architecture. "
            "Install AutoDock-Vina yourself and point the tests at it.")

    return f"vina_{VINA_VERSION}_{operating_system}_{architecture}"


def download(path: str = 'vina') -> str:
    """Download the AutoDock-Vina executable for the running platform.

    Parameters
    ----------
    path : str, optional
        Where to write the executable, by default 'vina' in the current
        directory.

    Returns
    -------
    str
        ``path``, for convenience.

    Raises
    ------
    RuntimeError
        If the download failed or produced an empty file.
    """
    asset = get_asset_name()
    url = f"{BASE_URL}/v{VINA_VERSION}/{asset}"
    try:
        urllib.request.urlretrieve(url, path)  # nosec
    except Exception as error:
        raise RuntimeError(f"Could not download AutoDock-Vina from {url}: {error}")

    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise RuntimeError(f"Downloaded AutoDock-Vina from {url} but {path} is empty.")

    # a+x, on top of whatever the file already has
    mode = os.stat(path).st_mode
    os.chmod(path, mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


if __name__ == "__main__":
    print(f"Downloaded {download()} ({get_asset_name()})", file=sys.stderr)
