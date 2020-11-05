#!/bin/bash

if [ "$1" = "coarse" ]
then
    echo "creating coarse mesh"
    cp utilities/meshes/coarse_mesh/coanda_mesh_cor.msh .
    gmshToFoam coanda_mesh_cor.msh
    cp utilities/meshes/coarse_mesh/boundary constant/polyMesh/
    rm coanda_mesh_cor.msh
    paraFoam &
else
    if [ "$1" = "refined" ]
    then
        echo "creating refined mesh"
        cp utilities/meshes/fine_mesh/coanda_mesh_ref.msh .
        gmshToFoam coanda_mesh_ref.msh
        cp utilities/meshes/fine_mesh/boundary constant/polyMesh/
        rm coanda_mesh_ref.msh
        paraFoam &
    else
        echo available option: coarse refined
    fi
fi
