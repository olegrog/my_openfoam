#!/bin/bash

file="constant/polyMesh/boundary"
old="$file.gmsh"
mv $file $old

sed '/physicalType/d' $old | awk '
BEGIN {
    type=""
}
{
    if ($1 == "frontAndBack")  type="empty";
    if ($1 == "bottom" || $1 == "side") type="wall";
    if ($1 == "type" && type) {
        print "        type            "type";";
        type="";
    } else print;
}' | sed '/\ttype/s/ //g' > $file
