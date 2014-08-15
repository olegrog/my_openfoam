PROJECT=noncoaxial
BOUNDARY=constant/polyMesh/boundary

all: $(BOUNDARY).gmsh

$(BOUNDARY).gmsh: $(BOUNDARY)
	./updatePolyMesh.sh

constant/polyMesh/boundary: $(PROJECT).msh
	gmshToFoam $(PROJECT).msh

$(PROJECT).msh: $(PROJECT).geo
	gmsh -3 $(PROJECT).geo

clean:
	. $(WM_PROJECT_DIR)/bin/tools/CleanFunctions && cleanCase
