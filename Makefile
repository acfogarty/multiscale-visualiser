#-----------------------------------------------------------------------
#  Look for "--->" for things that you need to change for your
#  particular configuration.
#-----------------------------------------------------------------------
#  This simple Makefile which compiles a single C++ program, in the
#  source file $TARGET.cpp (where the variable TARGET is set below)
#  assuming a Unix system for OpenGL, GLUT, and X-windows.  The directories
#  for OpenGL and GLUT will need to be changed depending on your system.
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# ---> INSERT THE NAME OF YOUR SOURCE FILE HERE
#-----------------------------------------------------------------------

TARGET = multiscale-visualiser

#-----------------------------------------------------------------------
#  Variable settings: Set these according to your system configuration
#
#  Compiler/Linker settings:
#  -------------------------
#	C++		Your favorite C++ compiler (or C compiler).
#	CFLAGS		Your compilation flags.
#	    -g		    compile for debugging (e.g., dbx or gdb).
#	    -O2/O3	    optimize code for speed
#	    -Wall	    give all warnings (highly recommended with
#			    OpenGL programs)
#	MISC_LIBS	Any miscellaneous libraries that your program
#			needs.
#			    add -lm here for math library
#			    add -lposix4 here for clock_gettime()
#			    (a real-time clock on Solaris machines).
#
#
#  System/Platform-dependent settings (which you will need to set):
#  ----------------------------------------------------------------
#	OGL_DIR		Directory where OpenGL (or Mesa) resides
#	OGL_LIBS	OpenGL library names
#	X_LIB		X11 lib directory (where libX11, libXext live)
#
#  System/Platform-dependent settings (you should not need to set):
#  ----------------------------------------------------------------
#	OGL_LIB_DIR	OpenGL library directory (where libGL lives)
#	OGL_INC_DIR	OpenGL linclude directory (where GL/gl.h lives)
#	OGL_LFLAGS	Loader flags for OpenGL stuff
#	X_LFLAGS	Loader flags for X11 stuff
#	LD_FLAGS	All loader flags
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# ---> USE THESE DEFINITIONS FOR CSIC Linux Cluster Machines
#-----------------------------------------------------------------------
OGL_DIR = /usr
OGL_LIBS =  -lGLU -lGL
GLUT_DIR = /usr/local/freeglut
GLUT_LIBS = -lglut
X_LIB = /usr/X11R6/lib

#-----------------------------------------------------------------------
# Setting for all UMD platforms
# ---> Set C++ to your C++ compiler
# ---> Set CFLAGS according your desired option (debug/performance)
#-----------------------------------------------------------------------
C++ = g++
# CFLAGS = -g -Wall				(for debugging)
# CFLAGS = -O3 -Wall				(for best performance)
CFLAGS = -O3 -Wall -std=c++11
MISC_LIBS =  -lm

OGL_LIB_DIR  = $(OGL_DIR)/lib
GLUT_LIB_DIR = $(GLUT_DIR)/lib
OGL_INC_DIR  = $(OGL_DIR)/include
GLUT_INC_DIR = $(GLUT_DIR)/include

OGL_LFLAGS   = -L$(OGL_DIR)/lib -L$(GLUT_DIR)/lib -Wl,-rpath=/usr/local/freeglut/lib $(OGL_LIBS) $(GLUT_LIBS) 
INC_FLAGS = -I$(OGL_INC_DIR) -I$(GLUT_INC_DIR) 
LD_FLAGS = $(OGL_LFLAGS) $(X_LFLAGS) $(MISC_LIBS)
X_LFLAGS = -L$(X_LIB) -lX11 -lXext -lXmu -lXi 

#-----------------------------------------------------------------------
# Compilation
#	This is the simplest form assuming a single source file.
#	If the file name is foo.cpp, set TARGET = foo, and the
#	object file will be put in foo.o, and the executable will
#	be put in foo.
#-----------------------------------------------------------------------

default: $(TARGET)

$(TARGET): $(TARGET).o camera.o structs.o
	$(C++) $(CFLAGS) $(TARGET).o camera.o structs.o -o $(TARGET) $(LD_FLAGS)

$(TARGET).o: $(TARGET).cpp camera.h structs.h camera.o structs.o
	$(C++) -c $(CFLAGS) $(INC_FLAGS) $(TARGET).cpp

camera.o: camera.cpp camera.h structs.cpp structs.h
	$(C++) -c $(CFLAGS) $(INC_FLAGS) camera.cpp

structs.o: structs.h
	$(C++) -c $(CFLAGS) $(INC_FLAGS) structs.cpp

#-----------------------------------------------------------------------
# ---> ADD THE FOLLOWING IF YOU HAVE ADDITIONAL SOURCE FILES, xxx.cpp.
#-----------------------------------------------------------------------

# xxx.o: xxx.cpp
# 	$(C++) -c $(CFLAGS) $(INC_FLAGS) xxx.cpp

#-----------------------------------------------------------------------
#  Clean-up
#	make clean		remove .o and core files
#	make realclean		also removes executable
#-----------------------------------------------------------------------
clean:
	rm -f *.o core

realclean: clean
	rm -f $(TARGET)
