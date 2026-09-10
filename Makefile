##################################################################
##################################################################
# User settings: modify these for your use case
##################################################################
##################################################################
# Host compiler choice (needs support for C++20)
HOST_COMPILER ?= g++


##################################################################
##################################################################
# Non-user-settings: You probably won't need to change these.
##################################################################
##################################################################
# File names and file paths for the program
program_NAME := main
src_DIR := src
program_C_SRCS := $(wildcard $(src_DIR)/*.c)
program_CXX_SRCS := $(wildcard $(src_DIR)/*.cpp)
program_H_SRCS := $(wildcard $(src_DIR)/*.h)
program_HPP_SRCS := $(wildcard $(src_DIR)/*.hpp)
program_C_OBJS := ${program_C_SRCS:.c=.o}
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}
program_CXX_ASMS := ${program_CXX_SRCS:.cpp=.s}

program_OBJS := $(program_C_OBJS) $(program_CXX_OBJS)
program_INCLUDE_DIRS := external/eigen external/boost-pfr/include
program_LIBRARY_DIRS :=
program_LIBRARIES := m dl quadmath openblas lapacke lapack


# Compiler flags
CXXFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
CXXFLAGS += -std=c++20 -Wall
#-fext-numeric-literals  	#-DEIGEN_HAS_CONSTEXPR=1 #-DEIGEN_NO_DEBUG
CXXFLAGS += -march=native -pthread -fopenmp
#CXXFLAGS += -march=native -pthread
CXXFLAGS += -O3 -ffast-math
#CXXFLAGS += -g -fno-omit-frame-pointer -fext-numeric-literals
CXXFLAGS += -DNDEBUG

# The standalone relaxation check deliberately excludes -ffast-math so that
# its finiteness and fixed-point checks retain standard IEEE semantics.
check_CXXFLAGS = $(filter-out -ffast-math,$(CXXFLAGS))
check_NAME := check_hydrostatic_relaxation
solver_check_OBJS := test/three_fluid.o test/evolution.o test/reproduction.o
check_OBJS := test/check_hydrostatic_relaxation.o $(solver_check_OBJS)


# Add linker flags
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir)) 
LDLIBS += $(foreach library,$(program_LIBRARIES),-l$(library))


.PHONY: all check clean distclean

all: $(program_NAME)

$(program_NAME): $(program_OBJS)
	$(LINK.cc) $(program_OBJS) -o $(program_NAME) $(LDLIBS)

check: $(check_NAME) check_statler
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./$(check_NAME)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./check_statler
	python3 -m unittest discover -s test -p 'test_*.py'

check_statler: test/check_statler.o $(solver_check_OBJS)
	$(CXX) $(check_CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main-strict: test/main.o $(solver_check_OBJS)
	$(CXX) $(check_CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

test/main.o: src/main.cpp $(program_HPP_SRCS)
	$(CXX) $(CPPFLAGS) $(check_CXXFLAGS) -c $< -o $@

test/check_statler.o: test/check_statler.cpp $(program_HPP_SRCS)
	$(CXX) $(CPPFLAGS) $(check_CXXFLAGS) -c $< -o $@

test/evolution.o: src/evolution.cpp $(program_HPP_SRCS)
	$(CXX) $(CPPFLAGS) $(check_CXXFLAGS) -c $< -o $@

test/reproduction.o: src/reproduction.cpp $(program_HPP_SRCS)
	$(CXX) $(CPPFLAGS) $(check_CXXFLAGS) -c $< -o $@

$(check_NAME): $(check_OBJS)
	$(CXX) $(check_CXXFLAGS) $(LDFLAGS) $(check_OBJS) -o $@ $(LDLIBS)

test/check_hydrostatic_relaxation.o: test/check_hydrostatic_relaxation.cpp $(program_H_SRCS) $(program_HPP_SRCS)
	$(CXX) $(CPPFLAGS) $(check_CXXFLAGS) -c $< -o $@

test/three_fluid.o: src/three_fluid.cpp $(program_H_SRCS) $(program_HPP_SRCS)
	$(CXX) $(CPPFLAGS) $(check_CXXFLAGS) -c $< -o $@

$(program_OBJS): $(program_H_SRCS) $(program_HPP_SRCS) $(program_CUH_SRCS) $(program_GEN_SRCS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.s: %.cpp
	$(CXX) $(CXXFLAGS) -S -fverbose-asm $< -o $@

asm: $(program_CXX_ASMS)

clean:
	$(RM) $(program_NAME)
	$(RM) $(program_OBJS)
	$(RM) $(check_NAME)
	$(RM) $(check_OBJS)
	$(RM) check_statler main-strict test/check_statler.o test/main.o
	$(RM) $(program_CXX_ASMS)
	$(RM) $(wildcard *~)
	$(RM) -r html latex

distclean: clean

show:
	echo $(CXX)
	echo $(GXX)
	echo $(GCC)
	echo $(LINK.cc)
	echo $(CC)
	echo $(CPP)
	echo $(RM)
	echo $(CXXFLAGS)
	echo $(NVCC)
	echo $(program_CXX_SRCS) "\n"
	echo $(program_HPP_SRCS) "\n"
	echo $(program_CXX_OBJS) "\n"
	echo $(program_OBJS) "\n"
	echo $(program_CU_SRCS) "\n"
	echo $(program_CUH_SRCS) "\n"
	echo $(program_CU_OBJS) "\n"
	echo $(device_link_OBJ) "\n"
	echo $(program_CXX_ASMS) "\n"
