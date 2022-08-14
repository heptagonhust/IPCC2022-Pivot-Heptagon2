all: pivot

SHELL := /bin/bash

pivot: pivot.cc defs.hpp
	[ ! -d ./build ] && mkdir build 
	source /public1/soft/modules/module.sh && module load oneAPI/2022.1 && icpx -Wall -Wno-unused-result -O3 -g -mavx2 -lpthread ./pivot.cc -o build/pivot