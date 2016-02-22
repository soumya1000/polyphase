
all: phase_sim

phase_sim: main.o emalgo.o particlefilter.o parsedata.o haplo.o globals.o sampling.o
#g++ main.o emalgo.o particlefilter.o sampling.o parsedata.o haplo.o globals.o -o phase_sim  -lgsl -lgslcblas
	g++ main.o emalgo.o particlefilter.o sampling.o parsedata.o haplo.o globals.o  -o phase_sim  -lgmpxx -lgmp  -lmpfr 

main.o: main.cpp
	g++ -std=c++11 -c  main.cpp	
	
sampling.o: sampling.cpp
	    g++ -std=c++11 -c sampling.cpp
	   
emalgo.o: emalgo.cpp
	#g++ -std=c++0x -c  emalgo.cpp -lgmpxx -lgmp
	g++ -std=c++11 -c emalgo.cpp 
	
particlefilter.o: particlefilter.cpp
	#g++ -std=c++0x -c particlefilter.cpp  -lgmpxx -lgmp
	g++ -std=c++11 -c particlefilter.cpp  

haplo.o: haplo.cpp
	g++ -std=c++11 -c haplo.cpp

parsedata.o: parsedata.cpp
	g++ -std=c++11 -c parsedata.cpp 

globals.o: globals.cpp
	#g++ -std=c++0x -c globals.cpp -lgmpxx -lgmp
	g++ -std=c++11 -c globals.cpp 

clean:
	rm -rf *.o phase_sim