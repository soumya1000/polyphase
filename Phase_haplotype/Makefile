CC = g++
CFLAGS = -std=c++11 -c  -I/TBB_INSTALL_DIR/include/
LDFLAGS = -lmpfr -L /TBB_INSTALL_DIR/build/linux_intel64_gcc_cc4.8_libc2.19_kernel3.16.7_release/ -ltbb
COMMON_SOURCES = sampling.cpp  emalgo.cpp particlefilter.cpp haplo.cpp parsedata.cpp globals.cpp
TARGET_SOURCES = main.cpp
INIT_SOURCES = init.cpp
COMMON_OBJECTS = $(COMMON_SOURCES:.cpp=.o)
TARGET_OBJECTS = $(TARGET_SOURCES:.cpp=.o)
INIT_OBJECTS = $(INIT_SOURCES:.cpp=.o)
EXECUTABLE =  phase_sim
INIT_EXECUTABLE =  init

.PHONY: all dakota_engine initializer

all: dakota_engine initializer

dakota_engine: $(EXECUTABLE)

initializer: $(INIT_EXECUTABLE)

$(EXECUTABLE): $(COMMON_OBJECTS) $(TARGET_OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@

$(INIT_EXECUTABLE): $(COMMON_OBJECTS) $(INIT_OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o phase_sim init