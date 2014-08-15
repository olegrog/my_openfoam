#!/bin/bash

file="constant/polyMesh/boundary"
old="$file.gmsh"
mv $file $old

sed '/physicalType/d' $old | awk '
BEGIN {
    type=""
} 
{
    if ($1 == "back" || $1 == "front")  type="empty";
    if ($1 == "symmetry")               type="symmetryPlane";
    if ($1 == "outer" || $1 == "inner") type="wall";
    if ($1 == "type" && type) {
        print "        type            "type";";
        type="";
    } else print;
}' | sed '/\ttype/s/ //g' > $file
