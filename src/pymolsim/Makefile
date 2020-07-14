CCcomp=g++
CC=gcc

ifdef VULCAN
	OBJS=  main.o  readinput.o ran3.o init.o forces.o analyse.o md.o print.o
else
	OBJS=  main.o  readinput.o ran3.o init.o forces.o analyse.o md.o draw3D.o print.o
endif

ifdef MAC
	LIBS=-lm -L/usr/X11R6/lib -lX11 -framework Cocoa -framework OpenGL -framework GLUT
	LIBS += -dylib_file \
		/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:\
		/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib
else

#	LIBS=-lm -L/usr/X11R6/lib -lX11 -lGL -lglut 
# for nvidia graphics card on new desktop
	LIBS=-lm -L/usr/lib/nvidia-331/ -lX11 -lGL -lglut 
endif

ifdef VULCAN
	LIBS=-lm -L/usr/X11R6/lib -lX11
endif


#CCFLAGS=-g -O2
#CCFLAGS=-g -O2 -Wall -Wconversion  #-g is for debugging
#CCFLAGS=-pg -O3 -Wall -Wconversion #calls the profiler
CCFLAGS=-O3 -Wall -Wconversion 

EXEflags=-O3
#EXEflags=-O3 -pg		#also needed for profiler

.SUFFIXES: .cpp .o 

.cpp.o:
	$(CCcomp) $(CCFLAGS) -c $<

all: main.exe 

main.exe: $(OBJS) 
	$(CCcomp) -o mdLJ.exe  $(EXEflags) $(OBJS) $(LIBS)
#	mv main.exe /tmp

clean:
	rm -f *.o
	rm -f mdLJ.exe
