My OpenFOAM solvers and problems
========

Currently, the latest version of [OpenFOAM+](https://www.openfoam.com) is supported.

Installation
-------
The repository can be cloned by
```
    git clone git@github.com:olegrog/openfoam.git
```
To build an augmented Docker image, run
```
    cd openfoam
    docker build -f Dockerfile.ubuntu -t my_openfoam-plus .
```
The new image can be used within an intermediate container as
```
    docker run -it --rm -v="$(pwd)":/home/ofuser my_openfoam-plus
```
