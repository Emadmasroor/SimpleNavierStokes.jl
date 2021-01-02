# Package SimpleNavierStokes

| Status | Coverage |
| :----: | :----: |
| [![Build Status](https://travis-ci.org/emadmasroor/SimpleNavierStokes.jl.svg?branch=master)](https://travis-ci.org/emadmasroor/SimpleNavierStokes.jl) | [![codecov.io](http://codecov.io/github/emadmasroor/SimpleNavierStokes.jl/coverage.svg?branch=master)](http://codecov.io/github/emadmasroor/SimpleNavierStokes.jl?branch=master) |

Hello! Welcome to this simple package. It is currently under active development.

## Installation Instructions

This package is not (yet) part of the official Julia package repository, `Pkg`. To add the package to your Julia environment, you'll need to enter the full package URL like this:

`] add https://github.com/emadmasroor/SimpleNavierStokes.jl`

in the Julia REPL. Then, type 

`using SimpleNavierStokes`.

That's it!

## Running your first simulation

Try the following code to get started:

```
sol1 = LidDrivenCavity()
``` 

This will solve the Navier-Stokes equations on a square 32x32 grid and store the results in a custom type, `SimpleNavierStokes.Results`. The first time you run this (inside a particular instance of the Julia REPL), it will be rather slow, while the functions compile. The next few times, it should be much faster.

To visualize the results, call the following function directly on `sol1`:

```
ShowStreamlines(sol1)
```

`ShowStreamlines()` is a custom function, also exported by `SimpleNavierStokes` (and therefore available for you to use), which operates on objects of type `Results`. This will visualize the steady-state flow in a square lid-driven cavity.

![Lid-driven cavity result](../assets/sample.png?raw=true)
