#-----------------------------------------------------------------------------------
#Compile Options
#-----------------------------------------------------------------------------------
FC := gfortran
EXE := OpenMeshConverter
EXT := .exe

all:			TYPE :=
MPparallel:    TYPE := _MP
debug:    TYPE := _Debug
catch:    TYPE := _Catch
coverage: TYPE := _Coverage

all:      OPT := -O3
MPparallel: OPT := -O3 -fopenmp
debug:    OPT := -O0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow
catch:    OPT := -O0 -g -Wall -W -Wsurprising -Werror
coverage: OPT := -O0 -g --coverage

#-----------------------------------------------------------------------------------
#Add source files as necessary
#If the files are not compiled using the generic rules, add commands under their
#rule declaration. Add these items to FILTER
#-----------------------------------------------------------------------------------
SRC := \
				main.f90 \
				globals.f90 \
				output_thrm.f90 \
				mesh_types.f90 \
				boundary_conditions.f90 \
				read_gmsh.f90

OBJ := $(SRC:.f90=.o)
MOD := $(OBJ:.o=.mod)

FILTER := \

OBJ_FILTER := $(FILTER:.f90=.o)
MOD_FILTER := $(FILTER:.f90=.mod)

#-----------------------------------------------------------------------------------
#Complex external dependencies. Each of these requires a custom build rule.
#Intended for dependencies with their own makefiles or other build system
#-----------------------------------------------------------------------------------
COMP_DEP := \
#Phony targets for cleaning and building

#-----------------------------------------------------------------------------------
#Phony targets for cleaning and building
#-----------------------------------------------------------------------------------
.PHONY: all clean reset

all: $(EXE)

MPparallel: $(EXE)

debug: $(EXE)

catch: $(EXE)

coverage: $(EXE)

#Intended to clean up compilation artifacts but leave executable & coverage
clean:
	rm -f $(OBJ)
	rm -f $(MOD)
	rm -f $(COMP_DEP)

#Intended to reset directory to fresh state with no exe or artifacts
reset: clean
	rm -f *.gcno *.gcda
	rm -f -r $(EXE)*.dSYM
	rm -f ../$(EXE)*

#-----------------------------------------------------------------------------------
#Generics for source files
#-----------------------------------------------------------------------------------
$(filter-out $(OBJ_FILTER), $(OBJ)): %.o:	%.f90
	$(FC) -c $(OPT) $<

$(filter-out $(MOD_FILTER), $(MOD)):	%.mod:	%.f90
	$(FC) -c $(OPT) $<

$(EXE): $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
	$(FC) -o $@$(TYPE)$(EXT) $(OPT) $(OBJ)
	mv ./$(EXE)$(TYPE)$(EXT) ../

#-----------------------------------------------------------------------------------
#Rules for entries in COMP_DEP. Don't forget to add them to make clean / reset
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
#Dependency List
#Use [gfortran -M -cpp *.f90] repeatedly until clean compile to update rules below
#-----------------------------------------------------------------------------------
boundary_conditions.o boundary_conditions.mod: boundary_conditions.f90 \
 globals.mod
globals.o globals.mod: globals.f90 mesh_types.mod
mesh_types.o mesh_types.mod: mesh_types.f90
main.o: main.f90 globals.mod read_gmsh.mod output_thrm.mod \
 boundary_conditions.mod
output_thrm.o output_thrm.mod: output_thrm.f90 globals.mod
read_gmsh.o read_gmsh.mod: read_gmsh.f90 globals.mod