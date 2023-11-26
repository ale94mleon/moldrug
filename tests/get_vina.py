import subprocess
import sysconfig


def run(command):
    subprocess.run(command, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def download():
    operating_system = sysconfig.get_platform().lower()

    if operating_system.startswith('macos'):
        run("""
            curl -L https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_mac_x86_64 > vina
            chmod a+x vina
        """)
    elif operating_system.startswith('linux'):
        run("""
        wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
        mv vina_1.2.5_linux_x86_64 vina
        chmod a+x vina
        """)
    else:
        raise NotImplementedError('Only MacOS and Linux support of the moment')


if __name__ == "__main__":
    download()
