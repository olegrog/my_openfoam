PROJECT=spheres-transfinite
BOUNDARY=constant/polyMesh/boundary

all: $(BOUNDARY).gmsh

$(BOUNDARY).gmsh: $(BOUNDARY)
	./updatePolyMesh.sh

constant/polyMesh/boundary: $(PROJECT).msh
	. $(WM_PROJECT_DIR)/bin/tools/RunFunctions && runApplication gmshToFoam $(PROJECT).msh

$(PROJECT).msh: $(PROJECT).geo
	gmsh -3 $(PROJECT).geo > log.gmsh 2>&1

clean:
	. $(WM_PROJECT_DIR)/bin/tools/CleanFunctions && cleanCase
	rm -f 0/wall* 0/grad*
	rm -f $(PROJECT).msh
	rm -rf constant/polyMesh
