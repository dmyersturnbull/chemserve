# Chemgrid

[![Version status](https://img.shields.io/pypi/status/chemgrid)](https://pypi.org/project/chemgrid/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/chemgrid)](https://pypi.org/project/chemgrid/)
[![Docker](https://img.shields.io/docker/v/dmyersturnbull/chemgrid?color=green&label=DockerHub)](https://hub.docker.com/repository/docker/dmyersturnbull/chemgrid)
[![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/dmyersturnbull/service-it?include_prereleases&label=GitHub)](https://github.com/dmyersturnbull/service-it/releases)
[[![Conda-Forge](https://img.shields.io/conda/vn/conda-forge/chemgrid?label=Conda-Forge)](https://anaconda.org/conda-forge/chemgrid)
![Latest version on PyPi](https://badge.fury.io/py/chemgrid.svg)](https://pypi.org/project/chemgrid/)
[![Documentation status](https://readthedocs.org/projects/service-it/badge/?version=latest&style=flat-square)](https://service-it.readthedocs.io/en/stable/)
[![Build & test](https://github.com/dmyersturnbull/service-it/workflows/Build%20&%20test/badge.svg)](https://github.com/dmyersturnbull/service-it/actions)
[![Maintainability](https://api.codeclimate.com/v1/badges/cb9bc2733ece01a8800b/maintainability)](https://codeclimate.com/github/dmyersturnbull/service-it/maintainability)
[![Coverage](https://coveralls.io/repos/github/dmyersturnbull/service-it/badge.svg?branch=master)](https://coveralls.io/github/dmyersturnbull/service-it?branch=master)

Turn any rdkit cheminformatics function into a daemon that listens on a port for JSON queries.
You can use it to make a REST API that returns converts SMILES to InChI, desalts molecules, calculates fingerprints,
or draws structures on a grid.

Your code receives a `Payload` with parameters (`.params`) and a list of compounds (`.items`),
each wrapped up into a class that lazily stores the SMILES, inchi, inchikey, name, mol object, and optional identifier.
Each has convenience functions for things like getting fingerprints and deduplicating structures.

Of course, rdkit only needs to be installed on the server.
You can safely import chemgrid on the client without rdkit.

**⚠ Status:** Not functional yet.

#### Examples

Here's an example that draws compounds on a grid:

```python

from chemgrid import ChemGrid, ChemGraphicsKit, BaseChem

def draw(payload):
    chems = payload.items
    rows, cols = payload.bool_param('rows'), payload.bool_param('cols')
    path = payload.str_param('path')
    grid = ChemGraphicsKit().draw_grid(chems, rows, cols)
    grid.save(path)

server = ChemGrid.server(draw, port=1633)

client = server.client()  # just shorthand

chems = [
    BaseChem('InChI=1S/H2O/h1H2', 'water', 'U0001')
]
client.send_multiple(chems, params=None)
```

This particular example is actually built-in:
Use `ChemGrid.drawing_service(1633)`.
Anyway, here's an example that desalts & deduplicates molecules and calculates ECFP fingerprints:

```python
from chemgrid import ChemGrid, BaseChem, Chem

def fingerprint(chem: Chem, params):
    return chem.desalt().deduplicate().fingerprint(radius=2, features=False)
# the 'simple' in simple_server just means it takes one compound at a time
server = ChemGrid.simple_server(fingerprint, port=1633)

# make a client (ordinarily this would happen on another process or remotely)
client = ChemGrid.client(1633)

# tell the client to handle replies
def callback(query, reply):
    print(query, reply)
client.add_callback(callback)

# send something
water = BaseChem('InChI=1S/H2O/h1H2', name='water', key='U0001')
client.send_one(water)
```


#### Original purpose: isolating rdkit

This was written to isolate needed functionality from rdkit
so that the rest of an [analysis package](https://github.com/dmyersturnbull/chemfish) could be installed through Poetry.
Otherwise, I had some dependencies only on PyPi and others only on Conda,
so each would overwrite the other's dependencies.

There are two components:
- The server, installable through conda. It just sits and waits.
- The client, which is just a wrapper to build and send JSON payloads.

The client is installable through pip/poetry. However, the daemon won't work this way!
The server depends on rdkit, so it cannot be installed through pip.
The conda recipe specifies rdkit as a dependency.


#### Contributing

[New issues](https://github.com/dmyersturnbull/chemgrid/issues) and pull requests are welcome.
Please refer to the [contributing guide](https://github.com/dmyersturnbull/chemgrid/blob/master/CONTRIBUTING.md).
Generated with [Tyrannosaurus](https://github.com/dmyersturnbull/tyrannosaurus).
