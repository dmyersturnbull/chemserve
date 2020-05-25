# Chemgrid

[![Version status](https://img.shields.io/pypi/status/chemgrid)](https://pypi.org/project/chemgrid/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/dmyersturnbull/chemgrid?include_prereleases&label=GitHub)](https://github.com/dmyersturnbull/chemgrid/releases)
[![Conda-Forge](https://img.shields.io/conda/vn/conda-forge/chemgrid?label=Conda-Forge)](https://anaconda.org/conda-forge/chemgrid)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/decorateme)](https://pypi.org/project/decorateme/)
[![Documentation status](https://readthedocs.org/projects/chemgrid/badge/?version=latest&style=flat-square)](https://chemgrid.readthedocs.io/en/stable/)

Tiny daemon that listens on a socket or port and desalts, deduplicates, and draws compound structures in a grid.

This was written to isolate needed functionality from rdkit
so that the rest of a [package of mine](https://github.com/dmyersturnbull/chemfish) could be installed through Poetry.
Otherwise, I had some dependencies only on PyPi and others only on Conda, causing problems.

There are two components:
- The server, installable through conda. It just sits and waits.
- The client, which is just a wrapper to build and send JSON payloads.

```python
Chemgrid.new(8, 12).add('InChI=1S/H2O/h1H2', 'water').send("water.svg")
```

This will send a payload like this:

```json
{
	"rows": "8",
	"columns": "12",
	"compounds": [
		{"inchi": "InChI=1S/H2O/h1H2", "smiles": "", "label": "water"}
	],
	"output": "water.svg"
}
```

The client is installable through pip/poetry. However, the daemon won't work this way!
The server depends on rdkit, so it cannot be installed through pip.
The conda recipe specifies rdkit as a dependency.


[New issues](https://github.com/dmyersturnbull/chemgrid/issues) and pull requests are welcome.
Please refer to the [contributing guide](https://github.com/dmyersturnbull/chemgrid/blob/master/CONTRIBUTING.md).
Generated with [Tyrannosaurus](https://github.com/dmyersturnbull/tyrannosaurus).
