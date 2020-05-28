# Please change this to your g++ location
CC=/usr/bin/g++
# Please change this to your nvcc location
NVCC=/usr/local/cuda/bin/nvcc
# Program name CPU
PROG_CPU:=fractal_CPU
# Program name CPU
PROG_GPU:=fractal_GPU



CPU: ${PROG_CPU} clean

SOURCES_CC:= ./source/general_functions.cc ./source/data_structures.cc ./source/CPU_functions.cc 
SOURCES_CU:= ./source/main.cu 
OBJECTS_CC:= $(SOURCES_CC:./source/%.cc=./bin/%_CPU.o)
OBJECTS_CU:= $(SOURCES_CU:./source/%.cu=./bin/%_CPU.o)


${PROG_CPU}: ${OBJECTS_CC} ${OBJECTS_CU}
	$(CC) -O2 -w -o $@ $^ -lboost_iostreams -lboost_system -lboost_filesystem -lutil

${OBJECTS_CC}: ./bin/%_CPU.o : ./source/%.cc FORCE
	$(CC) -x c++ -o $@ -O2 -w -c $< -lboost_iostreams -lboost_system -lboost_filesystem -lutil

${OBJECTS_CU}: ./bin/%_CPU.o : ./source/%.cu FORCE
	$(CC) -x c++ -o $@ -O2 -w -c $< -lboost_iostreams -lboost_system -lboost_filesystem -lutil


GPU: ${PROG_GPU} clean

SOURCES_CC:=  ./source/general_functions.cc ./source/data_structures.cc 
SOURCES_CU:=  ./source/main.cu ./source/CUDA_functions.cu 
OBJECTS_CC:= $(SOURCES_CC:./source/%.cc=./bin/%_GPU.o)
OBJECTS_CU:= $(SOURCES_CU:./source/%.cu=./bin/%_GPU.o)

fractal_GPU: ${OBJECTS_CU} ${OBJECTS_CC} 
	$(NVCC) --expt-relaxed-constexpr -o $@ $^ -O2 -w  -lboost_iostreams -lboost_system -lboost_filesystem -lutil -D GPU_MODE

$(OBJECTS_CC): ./bin/%_GPU.o : ./source/%.cc FORCE
	$(NVCC) --expt-relaxed-constexpr -o $@ -c $< -O2 -w -c $(SOURCES) -lboost_iostreams -lboost_system -lboost_filesystem -lutil -D GPU_MODE

$(OBJECTS_CU): ./bin/%_GPU.o : ./source/%.cu FORCE
	$(NVCC) --expt-relaxed-constexpr -o $@ -c $< -O2 -w -c $(SOURCES) -lboost_iostreams -lboost_system -lboost_filesystem -lutil -D GPU_MODE


.PHONY: FORCE
FORCE:

clean:
	rm -f ./bin/*.o