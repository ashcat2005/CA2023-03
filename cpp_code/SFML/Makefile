# make a rule to compile the program
# 	g++ -Isrc/include -c .\$(MAIN).cpp -fopenmp -Ofast
#	g++ .\$(MAIN).o -o $(MAIN) -Lsrc/lib -lsfml-graphics -lsfml-window -lsfml-system -lblas -fopenmp


# Make file
MAIN = liquid_SFML
MAIN_CPP = $(MAIN).cpp

MAIN_OUT=$(MAIN).out

#OpenMP
NUM_THR = 4

windows : windows_compile windows_run

windows_run :
	$(MAIN).exe	

windows_compile :
	g++ -Isrc/include -c .\$(MAIN).cpp -fopenmp -O3
	g++ .\$(MAIN).o -o $(MAIN) -Lsrc/lib -lsfml-graphics -lsfml-window -lsfml-system