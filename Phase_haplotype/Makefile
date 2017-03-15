CC = g++
CFLAGS= -std=c++11  -c -pthread -I  /run/media/root/System/tbb44_20150728oss/include/
LDFLAGS= -lmpfr -ltbb 
SOURCES = globals.cpp parsedata.cpp haplo.cpp emtf.cpp samplingtf.cpp  
TARGET_OBJECTS = $(SOURCES:.cpp=.o)

all: phase_sim 

phase_sim: $(TARGET_OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o phase_sim 