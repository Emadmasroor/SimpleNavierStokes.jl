# Package SimpleNavierStokes

| Status | Coverage |
| :----: | :----: |
| [![Build Status](https://travis-ci.org/emadmasroor/SimpleNavierStokes.jl.svg?branch=master)](https://travis-ci.org/emadmasroor/SimpleNavierStokes.jl) | [![codecov.io](http://codecov.io/github/emadmasroor/SimpleNavierStokes.jl/coverage.svg?branch=master)](http://codecov.io/github/emadmasroor/SimpleNavierStokes.jl?branch=master) |

Hello! Welcome to this simple package. It will export a function `LidDrivenCavity`, and in the future it will support other canonical flows.

## Instructions

Add the package to your Julia environment by executing `] add https://github.com/emadmasroor/SimpleNavierStokes.jl.git` in the Julia REPL. Then, type `using SimpleNavierStokes`.

That's it! You can now call the function `LidDrivenCavity()` to execute this classic benchmark problem on a 32 x 32 grid. There are optional keyword arguments, but you'll have to dig into the documentation to see the details.

