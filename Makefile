export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

PY_INC  = -I/usr/include/python2.7/
GSL_INC = -I/usr/local/include/
#OPT = -fopenmp -O4 -ffast-math 
#OPT = -fopenmp -O2
OPT = -O4 -ffast-math 

LTB: SOLE_DENSITY 
	g++ $(OPT) -o LTB scipy_tools.o splev_bispeu.o \
	tools.o LTB_background.o mass_function.o energy.o geodesics.o \
	model_explorer.o cmb.o sole_density.o -lgsl -lgslcblas -lpython2.7 -L.
SOLE_DENSITY: TOOLS LTB_BACKGROUND MASS_FUNCTION ENERGY GEODESICS MODEL_EXPLORER
	g++ $(OPT) -c $(PY_INC) $(GSL_INC) sole_density.cpp
MODEL_EXPLORER: TOOLS LTB_BACKGROUND MASS_FUNCTION ENERGY GEODESICS
	g++ $(OPT) -c $(PY_INC) $(GSL_INC) model_explorer.cpp

GEODESICS:	
	g++ $(OPT) -c $(PY_INC) $(GSL_INC) geodesics.cpp
ENERGY: TOOLS MASS_FUNCTION
		g++ $(OPT) -c $(PY_INC) $(GSL_INC) energy.cpp
MASS_FUNCTION: TOOLS
		g++ $(OPT) -c $(PY_INC) $(GSL_INC) mass_function.cpp
LTB_BACKGROUND: TOOLS
	g++ $(OPT) -c $(PY_INC) $(GSL_INC) LTB_background.cpp
TOOLS: SPLEV_BISPEU
	g++ $(OPT) -c $(PY_INC) $(GSL_INC) tools.cpp
SPLEV_BISPEU: 
	gfortran $(OPT) -c splev_bispeu.f90  
clean:
	rm *.o

