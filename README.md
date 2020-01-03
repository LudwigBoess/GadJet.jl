[![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](LICENSE.md)
[![Documentation Status](https://readthedocs.org/projects/gadjetjl/badge/?version=latest)](https://gadjetjl.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/LudwigBoess/GadJet.jl.svg?branch=master)](https://travis-ci.org/LudwigBoess/GadJet.jl)
[![codecov.io](https://codecov.io/gh/LudwigBoess/GadJet.jl/coverage.svg?branch=master)](https://codecov.io/gh/LudwigBoess/GadJet.jl?branch=master)

# GadJet.jl


Quickstart
==========

Reading Data
------------

If you want to read a simulation snapshot into memory with GadJet.jl, it's as easy as this:

```julia
    data = read_snap(filename)
```

This will return a dictionary with the header information in `data["Header"]` and the blocks sorted by particle type.

As an example, this is how you would access the positions of the gas particles:

```julia
    data["Parttype0"]["POS"]
```

If you only want to read a specific block for a single particle type (e.g. positions of gas particles) you can use the function with a specified blockname and particle type like so:

```julia
    pos = read_snap(filename, "POS", parttype=0)
```

This will return an array of the datatype of your simulation, usually Float32.

Quick Visualisation
-------------------

For a quick glimpse at your data you can use the glimpse function (yes, I thought hard about this one...)

```julia
    image = glimpse(filename)
```

This will return a 200x200 pixel image of the whole box, centered on the center of mass.

If you want to look at a specific range you can provide an array with the center coordinates as `center_pos = [x, y, z]` and the extent in x, y and z direction with `dx, dy, dz`.

```julia
    image = glimpse(filename, center_pos, dx, dy, dz)
```


Contributing
============

If you want to contribute to this project I would greatly appreciate your help!

Please feel free to contact me: lboess@usm.lmu.de

Stuff I want to improve in the near future:

- Speedup of the sph mapping.
- Add unit conversion.
