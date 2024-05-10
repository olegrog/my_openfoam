My OpenFOAM solvers and problems
========

Currently, the latest version of [OpenFOAM+](https://www.openfoam.com) is supported.

Installation
-------
The repository can be cloned by
```shell
    git clone https://github.com/olegrog/my_openfoam.git
```
To build an augmented Docker image, run
```shell
    cd openfoam
    docker build -f Dockerfile.ubuntu -t my_openfoam-plus .
```
The new image can be used within an intermediate container as
```shell
    docker run -it --rm -u="$(id -u):$(id -g)" -v="$(pwd):/home/openfoam" my_openfoam-plus
```

MacOS pre-installation
-------
By default, MacOS uses case-insensitive file system, which can cause name conflicts in OpenFOAM.
To create a case-sensitive volume, run
```shell
    sudo curl -o /usr/local/bin/openfoam-macos-file-system http://dl.openfoam.org/docker/openfoam-macos-file-system
    sudo chmod 755 /usr/local/bin/openfoam-macos-file-system
    sudo openfoam-macos-file-system -v my_openfoam create
    mkdir my_openfoam
    sudo openfoam-macos-file-system -v my_openfoam mount
```

ParaView can be installed by `brew install paraview` or manually from [the official site](https://www.paraview.org/download/).
In the latter case, it is worth to add these lines to `.bash_profile`:
```shell
    paraview_path=$(find /Applications -name 'ParaView-*' -maxdepth 1 | sort -V | tail -1)
    [ -x "$paraview_path/Contents/MacOS/paraview" ] || echo "ParaView is not found."
    export PATH="$PATH:$paraview_path/Contents/MacOS:$paraview_path/Contents/bin"
```

Shell commands
-------
It is convenient to add some functions to `.bash_profile` to run the Docker container (`ofp`) and open ParaView (`pf`):
```shell
ofp() {
    local dir="my_openfoam"
    local image_name="my_openfoam-plus"
    local user="openfoam"
    local timezone
    timezone="$(readlink /etc/localtime | sed 's_.*zoneinfo/__')"
    [ -d "$HOME/$dir/run" ] || sudo openfoam-macos-file-system -v $dir mount
    docker run -it --rm                   \
        --user="$(id -u):$(id -g)"        \
        --env USER=$USER                  \
        --env TZ="$timezone"              \
        --volume="$HOME/$dir:/home/$user" \
        $image_name
}

pf() {
    local file="system/controlDict" dir
    [ -f "$file" ] || { echo "There is no $file"; return 1; }
    dir=$(find . -maxdepth 1 -name "*.foam")
    [[ $dir ]] || { dir="$(pwd | xargs basename).foam"; touch "$dir"; }
    paraview "$dir" > /dev/null 2>&1 &
}
```

