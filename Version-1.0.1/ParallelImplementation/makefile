# makefile Author: Antonio Carrillo Ledesma 29/12/2010
#                  antonio@mmc.geofisica.unam.mx
#
# usage:
#   First automatic dependency calculation
#     $ make deps
#
#   Compile and link
#     $ make
#
#   Run
#     $ make run
#
#   Clean
#     $ make clean
#
#   Others
# 	target: gz - Store in tar.gz
# 	target: bz2 - Store in tar.bz2
# 	target: tar - Store in tar
# 	target: zip - Store in zip
#	target: ddd - Debuger with ddd
# 	target: val - Debuger with valgrind
# 	target: gpro - Analysis of speedup
# 	target: svnup - svn update
# 	target: svnst - svn status -verbose
# 	target: svnci - svn commit
# 	target: gitup - svn update
# 	target: gitst - svn status
# 	target: gitci - svn commit
# 	target: help - Display callable targets.




# Defining the sources
SOURCES = $(wildcard *.cpp)

# Parameters of executable
PARAMETERS = file LaplaceSinpinxSinpiny -p 3 -m PMF1 -Nx 2 -Mx 256

# Defining the numbers of CPUs
CPUS = -np 3

# Defining the neme of th executable
EXECUTABLE = DVS-DDM

# Defining the debug Flags
#DEBUG = -O0 -g -pg                         # Debuger
#DEBUG += -pedantic-errors -Wall	-Wextra -O # ANSI C++
DEBUG += -O2	                             # Optimization O2

# Defining the compiler Flags
CFLAGS =  $(DEBUG)

# Defining the linker Flags
LDFLAGS =  $(DEBUG)

# Defining the compiler:
#CC = mpiCC.mpich                  # MPICH 1.X
CC = mpic++.mpich2                # MPICH 2.X
#CC = mpic++                       # MPI

# Defining the run the application
#RUN = mpirun.mpich                # MPICH 1.X
RUN =	mpirun.mpich2             # MPICH 2.X
#RUN = mpirun                      # MPI



# Defining the object files:
OBJECTS = $(SOURCES:.cpp=.o)

# Load the automatic dependency
DEPS = $(patsubst %.o,%.d,$(OBJECTS))


# target: all - The default rule - compiling our main program:
all: $(SOURCES)  $(EXECUTABLE) 
	echo all: make complete

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

# Tell make how to build .o files from .cpp files:
$(OBJECTS):
	$(CC) -MD -c $(CFLAGS)  $(@:.o=.cpp) -o $@

#Now make sure that make rebuilds files if included headers change:
-include $(DEPS)




# target: deps - Automatic dependency calculation
deps: $(SOURCES)
	$(CC) -MD -E $(SOURCES) > /dev/null



# The .PHONY rule keeps make from doing something with a file named clean
.PHONY: clean

# target: clean - Removing the executable and the object files
clean:
	rm $(EXECUTABLE) $(OBJECTS) $(DEPS) gmon.out > /dev/null 2>&1
	echo clean: make complete


# target: run - Run the application
run:
	lamboot -v
	$(RUN) $(CPUS) ./$(EXECUTABLE) $(PARAMETERS)
	lamhalt -v


# target: gz - Store in tar.gz
gz:
	tar -zcvpf $(EXECUTABLE)-`date "+%F"`.tar.gz makefile *.?pp
	

# target: bz2 - Store in tar.bz2
bz2:
	tar -jcvpf $(EXECUTABLE)-`date "+%F"`.tar.bz2 makefile *.?pp


# target: tar - Store in tar
tar:
	tar -cvpf $(EXECUTABLE)-`date "+%F"`.tar makefile *.?pp


# target: zip - Store in zip
zip:
	zip  $(EXECUTABLE)-`date "+%F"`.zip makefile *.?pp



# target: val - Debuger with valgrind, required DEBUG = -g 
val:
	lamboot -v
	valgrind -v --tool=memcheck --leak-check=full --show-reachable=yes $(RUN) $(CPUS) ./$(EXECUTABLE) $(PARAMETERS)
	lamhalt -v


# target: gpro - Analysis of speedup, required DEBUG = -g -pg
gpro:
	lamboot -v
	$(RUN) $(CPUS) ./$(EXECUTABLE) $(PARAMETERS)
	gprof -c -z ./$(EXECUTABLE)
	lamhalt -v


# target: svnup - svn update
svnup:
	svn update


# target: svnst - svn status
svnst:
	svn status --verbose


# target: svnci - svn commit
svnci:
	svn ci -m "Update"


# target: gitup - git update
gitup:
	git pull


# target: gitst - git status
gtist:
	git status


# target: gitci - git commit
gitci:
	git commit -a -m 'Update'
	git push
   
# target: help - Display callable targets.
help:
	egrep "^# target:" [Mm]akefile




