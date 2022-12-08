My OpenFOAM solvers and problems
========

Currently, the latest version of [OpenFOAM+](https://www.openfoam.com) is supported.

Installation
-------
The repository can be cloned by
```
    git clone https://github.com/olegrog/my_openfoam.git
```
To build an augmented Docker image, run
```
    cd openfoam
    docker build -f Dockerfile.ubuntu -t my_openfoam-plus .
```
The new image can be used within an intermediate container as
```
    docker run -it --rm -u="$(id -u):$(id -g)" -v="$(pwd):/home/openfoam" my_openfoam-plus
```

Mac OS pre-installation
-------
By default, Mac OS uses case-insensitive file system, which can cause name conflicts in OpenFOAM.
To create a case-sensitive volume, run
```
    sudo curl -o /usr/local/bin/openfoam-macos-file-system http://dl.openfoam.org/docker/openfoam-macos-file-system
    sudo chmod 755 /usr/local/bin/openfoam-macos-file-system
    sudo openfoam-macos-file-system -v my_openfoam create
    mkdir my_openfoam
    sudo openfoam-macos-file-system -v my_openfoam mount
```
