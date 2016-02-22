
all: init
#build phase-sim separately
init: init.cpp main_d.o emalgo_d.o particlefilter_d.o parsedata_d.o haplo_d.o globals_d.o sampling_d.o 
	g++ -std=c++11 -g init.cpp -I/run/media/root/System/tbb44_20150728oss/include/ emalgo_d.o particlefilter_d.o sampling_d.o parsedata_d.o haplo_d.o globals_d.o -o init -lmpfr -L /run/media/root/System/tbb44_20150728oss/build/linux_intel64_gcc_cc4.8_libc2.19_kernel3.16.7_debug/ -ltbb

main_d.o: main.cpp
	g++ -std=c++11 -c -g main.cpp -I/run/media/root/System/tbb44_20150728oss/include/ -o main_d.o
	
sampling_d.o: sampling.cpp
	    g++ -std=c++11 -c -g sampling.cpp -I/run/media/root/System/tbb44_20150728oss/include/ -o sampling_d.o
	   
emalgo_d.o: emalgo.cpp
	#g++ -std=c++0x -c  emalgo.cpp -lgmpxx -lgmp
	g++ -std=c++11 -c -g emalgo.cpp -I/run/media/root/System/tbb44_20150728oss/include/ -o emalgo_d.o
	
particlefilter_d.o: particlefilter.cpp
	#g++ -std=c++0x -c particlefilter.cpp  -lgmpxx -lgmp
	g++ -std=c++11 -c -g particlefilter.cpp  -I/run/media/root/System/tbb44_20150728oss/include/ -o particlefilter_d.o

haplo_d.o: haplo.cpp
	g++ -std=c++11 -c -g haplo.cpp -I/run/media/root/System/tbb44_20150728oss/include/ -o haplo_d.o

parsedata_d.o: parsedata.cpp
	g++ -std=c++11 -c -g parsedata.cpp -I/run/media/root/System/tbb44_20150728oss/include/ -o parsedata_d.o

globals_d.o: globals.cpp
	#g++ -std=c++0x -c globals.cpp -lgmpxx -lgmp
	g++ -std=c++11 -c -g globals.cpp -I/run/media/root/System/tbb44_20150728oss/include/ -o globals_d.o

clean:
	rm -rf *.o phase_sim init