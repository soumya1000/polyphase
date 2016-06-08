
CC = g++
CFLAGS = -std=c++11 -c  -I/run/media/root/System/tbb44_20150728oss/include/
LDFLAGS =-lgsl -lgslcblas -lmpfr -ltbb
COMMON_SOURCES = sampling.cpp  emalgo.cpp particlefilter.cpp haplo.cpp parsedata.cpp globals.cpp numoptimisation.cpp
TARGET_SOURCES = main.cpp
COMMON_OBJECTS = $(COMMON_SOURCES:.cpp=.o)
TARGET_OBJECTS = $(TARGET_SOURCES:.cpp=.o)
EXECUTABLE =  phase_sim

all: phase_engine 

phase_engine: $(EXECUTABLE)

$(EXECUTABLE): $(COMMON_OBJECTS) $(TARGET_OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o phase_sim init