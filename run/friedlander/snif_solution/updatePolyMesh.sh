#!/bin/bash

file="constant/polyMesh/boundary"
old="$file.gmsh"
mv $file $old

sed '/physicalType/d' $old | awk '
BEGIN {
    type=""
} 
{
    if ($1 == "frontAndBack")   type="empty";
    if ($1 == "right")          type="wall";
    if ($1 == "left")           type="wall";
    if ($1 == "symmetry")       type="symmetryPlane";
    if ($1 == "cold")           type="wall";
    if ($1 == "hot")            type="wall";
    if ($1 == "short")          type="wall";
    if ($1 == "long")           type="wall";
    if ($1 == "type" && type) {
        print "        type            "type";";
        type="";
    } else print;
}' | sed '/\ttype/s/ //g' > $file
