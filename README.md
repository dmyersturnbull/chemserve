# Chemgrid-service

[![Version status](https://img.shields.io/pypi/status/chemgrid-service)](https://pypi.org/project/chemgrid-service/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/dmyersturnbull/chemgrid-service?include_prereleases&label=GitHub)](https://github.com/dmyersturnbull/chemgrid-service/releases)
[![Conda-Forge](https://img.shields.io/conda/vn/conda-forge/chemgrid-service?label=Conda-Forge)](https://anaconda.org/conda-forge/chemgrid-service)
[![Documentation status](https://readthedocs.org/projects/chemgrid-service/badge/?version=latest&style=flat-square)](https://chemgrid-service.readthedocs.io/en/stable/)

Tiny daemon that listens on a socket or port and desalts, deduplicates, and draws compound structures in a grid.
Can output to SVG, PNG, or PDF.

This was written to isolate needed functionality from rdkit
so that the rest of a [package of mine](https://github.com/dmyersturnbull/chemfish) could be installed purely through Poetry.
Otherwise, I had some dependencies only on PyPi and others only on Conda, which caused problems.

Depends on rdkit, so it can be installed through conda but not pip.
You can build it with Poetry (or just pip install), but it won't work without rdkit.
The conda recipe specifies rdkit as a dependency.


[New issues](https://github.com/dmyersturnbull/chemgrid-service/issues) and pull requests are welcome.
Please refer to the [contributing guide](https://github.com/dmyersturnbull/chemgrid-service/blob/master/CONTRIBUTING.md).
Generated with [Tyrannosaurus](https://github.com/dmyersturnbull/tyrannosaurus) (but stripped down for conda).
