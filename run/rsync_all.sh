#!/bin/bash

srv=hil

dir=$(pwd)
dir=${dir#$HOME/}
dest=_fine
[[ $# > 0 ]] && dest=$1

mkdir -p $HOME/$dir/$dest
cases=$(ssh $srv "cd $dir; ls | grep ^[0-9]")
cases=.

for c in $cases; do
    echo "Copying $c to $dest..."
    rsync -auvz $srv:$dir/$c/{system,constant} $HOME/$dir/$dest/$c/
    rsync -auvz $srv:$dir/$c/0.0* $HOME/$dir/$dest/$c/
    rsync -auvz $srv:$dir/$c/log.* $HOME/$dir/$dest/$c/
done
