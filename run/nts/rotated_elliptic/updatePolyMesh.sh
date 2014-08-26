#!/bin/bash

file="constant/polyMesh/boundary"
old="$file.gmsh"
mv $file $old

sed '/physicalType/d' $old | awk '
BEGIN {
    type=""
    name=""
} 
{
    if ($1 == "back" || $1 == "front")  type="empty";
    if ($1 == "symmetry")               type="symmetryPlane";
    if ($1 == "outer" || $1 == "inner") type="wall";
    if ($1 == "type" && type) {
        print "        type            "type";";
        if (type == "cyclic") {
            if (name == "top") {
                neighbour = "bottom";
            } else {
                neighbour = "top";
            }
            print "        neighbourPatch  "neighbour";";
            print "        matchTolerance  0;";
            print "        transform       rotational;";
            print "        rotationAxis    (0 0 1);";
            print "        rotationCentre  (0 0 0);";
            print "        rotationAngle   180;";
        }
        type="";
    } else print;
}' | sed '/\ttype/s/ //g' > $file

exit 0
            print "        rotationAngle   180;";
    if ($1 == "top" || $1 == "bottom")  { type="cyclic"; name=$1 }
