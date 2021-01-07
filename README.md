# chemserve

[![Version status](https://img.shields.io/pypi/status/chemserve)](https://pypi.org/project/chemserve/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/chemserve)](https://pypi.org/project/chemserve/)
[![Docker](https://img.shields.io/docker/v/dmyersturnbull/chemserve?color=green&label=DockerHub)](https://hub.docker.com/repository/docker/dmyersturnbull/chemserve)
[![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/dmyersturnbull/service-it?include_prereleases&label=GitHub)](https://github.com/dmyersturnbull/service-it/releases)
[![Conda-Forge](https://img.shields.io/conda/vn/conda-forge/chemserve?label=Conda-Forge)](https://anaconda.org/conda-forge/chemserve)
[![Latest version on PyPi](https://badge.fury.io/py/chemserve.svg)](https://pypi.org/project/chemserve/)
[![Documentation status](https://readthedocs.org/projects/service-it/badge/?version=latest&style=flat-square)](https://service-it.readthedocs.io/en/stable/)
[![Build & test](https://github.com/dmyersturnbull/service-it/workflows/Build%20&%20test/badge.svg)](https://github.com/dmyersturnbull/service-it/actions)
[![Maintainability](https://api.codeclimate.com/v1/badges/cb9bc2733ece01a8800b/maintainability)](https://codeclimate.com/github/dmyersturnbull/service-it/maintainability)
[![Coverage](https://coveralls.io/repos/github/dmyersturnbull/service-it/badge.svg?branch=master)](https://coveralls.io/github/dmyersturnbull/service-it?branch=master)

Easily craft an rdkit [facade](https://en.wikipedia.org/wiki/Facade_pattern), resulting in conda-free code.

You make two uncoupled packages:
- A facade package with a conda dependency on rdkit. It can run as a daemon or server and respond to JSON requests.
- A pip-installable main codebase that can communicate with the server.


#### Advantages / rationale:

This has a few extremely useful properties:

- **Modularity:** Wrapping complex systems (like rdkit) in a [facade](https://en.wikipedia.org/wiki/Facade_pattern)
  is a great idea. Tests for your client code in fact *can’t* depend on rdkit, leaving them uncoupled.
- **Conda independence:** Rdkit is not on PyPi†, but lots of packages are only on PyPi, which has approximately
  10× the number of packages.  
  Unfortunately, mixing them can [silently break](https://dmyersturnbull.github.io/#-the-python-build-landscape)
  Conda’s solver. Your facade only needs rdkit, and your client can use any PyPi package.
- **Performance:** Conda is heavy, and its dependency resolution is slow (and less reliable) compared to
  [Poetry’s](https://python-poetry.org).
  Your facade’s public functions likely only delegate to a few rdkit calls, so updates can probably be infrequent.
  Meanwhile, you can make frequent changes to the client with speedy testing.
- **Tooling:** Without Conda, you can use modern build, isolation, and CI/CD tools like Poetry.  
  (See [Tyrannosaurus](https://github.com/dmyersturnbull/tyrannosaurus) for an example of a modern CI/CD pipeline
  built with Poetry and [Github Actions](https://github.com/features/actions).
  It also shows how to [generate a Conda-Forge recipe](https://tyrannosaurus.readthedocs.io/en/stable/anaconda.html)
  for your facade package.)
- **Asynchronous code:** Asynchronous calls can be made to numerically heavy rdkit code:
  Your client sends a request and independently handles payloads whenever they’re shipped back.

† Rdkit has an [open issue](https://github.com/rdkit/rdkit/issues/1812) about pip installation.
Unfortunately, the Pip package would likely still be massive and slow to install.


#### Baked-in functionality:

This package provides some [convenient](https://github.com/dmyersturnbull/chemserve/blob/master/chemserve/models.py)
[wrapper](https://github.com/dmyersturnbull/chemserve/blob/master/chemserve/_concrete_base.py)
[classes](https://github.com/dmyersturnbull/chemserve/blob/master/chemserve/graphics.py) to help build your facade.
You can pretty easily turn your daemon into a REST server and handle web requests.
There’s a baked-in server that generates structure SVGs and saves them to disk.

It’s all simpler than it sounds: See the examples below.


#### Example 1: draw compounds in a grid

Your server receives a `Payload` with parameters (`.params`) and a list of compounds (`.items`),
each wrapped up into a class that lazily stores the SMILES, inchi, inchikey, name, mol object, and optional identifier.
Each has convenience functions for things like getting fingerprints and deduplicating structures.

Here’s an example that draws compounds on an m×n grid.
(This particular example is actually built-in as `chemserve.drawing_service(port)`.)

```python

from chemserve import chemserve, ChemGraphicsKit, BaseChem, Payload


def draw(payload: Payload) -> None:
  """
  Draws the compounds in a grid and saves SVG or PNG to disk.
  The payload should have params ``rows``, ``cols``, and ``path``.
  """
  rows, cols = payload.bool_param("rows"), payload.bool_param("cols")
  path = payload.str_param("path")
  grid = ChemGraphicsKit().draw_grid(payload.compounds, rows, cols)
  grid.save(path)

server = chemserve.server(draw, port=1633)
client = server.client()  # just shorthand

client.send_multiple([BaseChem("InChI=1S/H2O/h1H2", "water", "U0001")], params=None)
```

#### Example 2: Calculate fingerprints and send them back

Here’s an example that desalts & deduplicates molecules and calculates ECFP fingerprints:

```python
from chemserve import chemserve, BaseChem, Chem

def fingerprint(chem: Chem, params):
    return chem.desalt().deduplicate().fingerprint(radius=2, features=False)
# the "simple" in simple_server just means it takes one compound at a time
server = chemserve.simple_server(fingerprint, port=1633)

# make a client (ordinarily this would happen on another process or remotely)
client = chemserve.client(1633)

# tell the client to handle replies
def callback(query, reply):
    print(query, reply)
client.add_callback(callback)

# send something
water = BaseChem("InChI=1S/H2O/h1H2", name="water", key="U0001")
client.send_one(water)
```


#### The situation that sparked this package

This was written to isolate needed functionality from rdkit
so that the rest of an [analysis package](https://github.com/dmyersturnbull/sauronlab)
could be installed through Poetry.
Otherwise, I had some dependencies only on PyPi and others only on Conda,
so each would overwrite the other’s dependencies.

There are two components:
- The server, installable through conda. It just sits and waits.
- The client, which is just a wrapper to build and send JSON payloads.

The client is installable through pip/poetry. However, the daemon won’t work this way!  
The server depends on rdkit, so it cannot be installed through pip.
The conda recipe specifies rdkit as a dependency.


#### Contributing

[New issues](https://github.com/dmyersturnbull/chemserve/issues) and pull requests are welcome.
Please refer to the [contributing guide](https://github.com/dmyersturnbull/chemserve/blob/master/CONTRIBUTING.md)
and [security policy](https://github.com/dmyersturnbull/chemserve/blob/master/SECURITY.md).  
Generated with [Tyrannosaurus](https://github.com/dmyersturnbull/tyrannosaurus).
