#!/bin/bash


cp utilities/U 0/U

case $1 in
"coarse") echo "creating coarse mesh"
          cp utilities/meshes/coarse_mesh/coanda_mesh_cor.msh .
          gmshToFoam coanda_mesh_cor.msh
          cp utilities/meshes/coarse_mesh/boundary constant/polyMesh/
          rm coanda_mesh_cor.msh
          paraFoam &
          ;;
"refined") echo "creating refined mesh"
           cp utilities/meshes/coarse_mesh/coanda_mesh_cor.msh .
           gmshToFoam coanda_mesh_cor.msh
           cp utilities/meshes/coarse_mesh/boundary constant/polyMesh/
           rm coanda_mesh_cor.msh
           paraFoam &
           ;;
"structured") echo "creating structured mesh"
              blockMesh
              ;;
*)  echo "invalid option: possible options= coarse    refined   structured"
    ;;
esac


