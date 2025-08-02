# 💿 Installation

## Via pip (standard)

## Getting the last "stable" version

Strongly recommended; although this project is still in beta state.

```bash
pip install moldrug
```

````{note}
To enable calculation of a cost function on [a cluster](https://jobqueue.dask.org/en/latest/clusters-api.html) install optional set of dependencies

```sh
pip install moldrug[cluster]
```
````

````{admonition} Linux 🐧
:class: Tip

Recent Linux distributions require creation of virtual environment before calling `pip`.  
Additionally the default Python might be too new for Moldrug. Python 3.12 is not yet supported (see pull requests [#9](https://github.com/ale94mleon/moldrug/pull/9) and [#10](https://github.com/ale94mleon/moldrug/pull/10)).

For example on Ubuntu 24.04 LTS, one needs to perform following steps

1. Get older Python
   ```sh
   sudo add-apt-repository ppa:deadsnakes/ppa
   sudo apt update
   sudo apt install python3.11-venv
   ```
2. Create a virtual environment
   ```sh
   python3.11 -m venv .venv
   source .venv/bin/activate
   ```
3. Only then call `pip`
   ```sh
   pip install moldrug
   ```

The other Linux distributions might require installation of older Python and creation of virtual environment as well.
````

## Getting the development version

```bash
pip install -U git+https://github.com/ale94mleon/moldrug.git@main
```

## Via conda

moldrug is also available through conda. However, the pip installation is the recommended one.

```bash
conda create -n moldrug -c ale94mleon -c conda-forge moldrug
```

```{note}
MacOS users may face some problems trying to install because of the AutoDock-Vina dependency. If that is the case, just download the executable as it is explained next.
```

If some dependencies are missing, please install them through pip. Some of them might be:

```bash
pip install meeko crem pyyaml scipy tqdm
```

## AutoDock-Vina dependency and related software

[AutoDock-Vina](https://vina.scripps.edu/) is the only non-pip dependency required for `moldrug`. [AutoDock-Vina](https://vina.scripps.edu/) is only relevant when the built-in cost functions of moldrug from {py:mod}`moldrug.fitness` are used. [AutoDock-Vina](https://vina.scripps.edu/) is an ongoing project, and it is advisable to stay up-to-date by regularly checking for the latest [release](https://github.com/ccsb-scripps/AutoDock-Vina/releases/).

As of the creation of this documentation, the most recent version is [v1.2.5](https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.5). We will be using this version as a demonstration, but we strongly recommend using the most recent release for optimal performance and features. For detailed information, please visit the
[AutoDock-Vina installation instruction](https://autodock-vina.readthedocs.io/en/latest/installation.html).

````{tab} Linux 🐧
```bash
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
chmod a+x vina_1.2.5_linux_x86_64
./vina_1.2.5_linux_x86_64 --help
```
````

````{tab} MacOS 🍏
```bash
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_mac_x86_64
chmod a+x vina_1.2.5_mac_x86_64
./vina_1.2.5_mac_x86_64 --help
```

```{note}
MacOS users might need to allow the execution of the application on ``Privacy & Security`` depending on the MacOS version.
```
````

````{tab} Windows 🪟
Please, download from [release](https://github.com/ccsb-scripps/AutoDock-Vina/releases/). Conda installation may not work.
````

## Converting pdb to pdbqt for the receptor

This step can be achieved through [OpenBabel](https://github.com/openbabel/openbabel) or [ADFR](https://ccsb.scripps.edu/adfr/downloads/). We recommend ADFR. Depending on the platform, you should be able to access the program `prepare_receptor` in different ways. In my case, it lies on `/Users/$USER/ADFRsuite-1.0/bin/prepare_receptor`. Then you can convert your ``pdb`` with:

```bash
/Users/$USER/ADFRsuite-1.0/bin/prepare_receptor -r your_protein.pdb -o your_protein.pdbqt
```

Check [this tutorial](https://ccsb.scripps.edu/adfr/how-to-create-a-pdbqt-for-my-receptor/) for more information. MacOS users may need to use ADFRsuite inside a docker container to avoid Gatekeeper's restrictions (see [this](https://alejandro.netlify.app/post/adfr/) article).

## Getting docking box information

To perform docking you must provide `boxcenter` (`center` for AutoDock-Vina) and `boxsize` (`size` for AutoDock-Vina) to the cost functions defined in {py:mod}`moldrug.fitness`. For that, two PyMOl plugins are useful: [GetBox](https://raw.githubusercontent.com/ale94mleon/GetBox-PyMOL-Plugin/refs/heads/master/GetBox.py) and/or [autodock](https://github.com/ADplugin/ADplugin/blob/master/autodock.py). Details of their installation and use are not discussed here, please visit their corresponding repositories for more information.

## Working with a docker container

- Use the [Docker configuration file on GitHub](https://github.com/ale94mleon/moldrug/blob/main/Dockerfile).
- Visit the [moldrug's Docker container](https://hub.docker.com/r/ale94mleon/4moldrug).

Finally, `pip install moldrug` inside the container.

## Running tests

```sh
pip install -e .[dev]
pytest tests
```

To test execution on a cluster, change the `execution_mode` variable in `test_moldrug.py`. The tests you can run on the cluster afterward are called `test_single_receptor_command_line` and `test_multi_receptor`. Each cluster is different, and the current configuration is tailored to our two internal clusters. To test on your own cluster, you’ll also need to modify the `SLURMCluster` configuration in `test_moldrug.py`. During the execution of each test, multiple entries should appear in the Slurm queue.
