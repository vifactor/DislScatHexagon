# Simple Makefile

source_dirs:= .  
search_wildcards:= $(addsuffix /*.cpp, $(source_dirs))
source_files:= $(wildcard $(search_wildcards))
source_files:=$(notdir $(source_files))
object_files:=$(patsubst %.cpp, %.o, $(source_files))
include_dirs:=$(addprefix -I, $(source_dirs)) -I/users/victor/libs/include/
release_flags:=  -O3 -Wall 
link_keys:= -lgsl -lcblas -lconfig++ 
link_dirs:= -L/users/victor/libs/lib/

DislScatHexagonBoost: $(object_files)
	mpiicpc $(link_dirs)  $(object_files) -o $@ $(link_keys) 
	mv $@ ../..
	rm *.o *.d
 
clean:
	rm *.o *.d

VPATH:= $(source_dirs)

%.o: %.cpp
	mpiicpc -c -MD $(release_flags)  $(include_dirs) $<

include $(wildcard *.d)
