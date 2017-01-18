# Compile-time flags
#CFLAGS = -pg -g -Wall -DDMALLOC -DDMALLOC_FUNC_CHECK -ldmalloc
#CFLAGS = -g -Wall -O0
#CFLAGS = -Wall -O3
CFLAGS = -O3

#Epiphany: Intel C compiler 
#CC = /opt/intel_cce_80/bin/icc -O3
CC = /usr/bin/gcc $(CFLAGS)

# Link-time flags
#LFLAGS = -O0

# executable name
TARGET = mendel
OBJECTS = mendel.o init.o sort.o fileio.o selection.o offspring.o \
          diagnostics.o ranlib.o mem.o mpi.o

##########################################
# build rules
##########################################

# target program
$(TARGET) : $(OBJECTS)
	$(CC) $(LFLAGS) -o $(TARGET) $(OBJECTS) -lm -lmpich -lpthread

serial : $(OBJECTS)
	$(CC) $(LFLAGS) -o $(TARGET) $(OBJECTS) -lm 

commit:
	cvs commit ranlib.h selection.c sort.c offspring.c mpi.c fileio.c mendel.c diagnostics.c init.c mem.c ranlib.c mendel.h Makefile

diff:
	cvs diff *.c *.h Makefile

run:
	./mendel

clean:
	\rm -f *.o *.mod $(TARGET) test00.*

distclean:
	\rm -f $(OBJECTS) $(TARGET) $(DEPEND_FILES)


###########################################
# dependencies
###########################################

mendel.o:	mendel.c mendel.h
	$(CC) -c mendel.c

init.o:		init.c mendel.h
	$(CC) -c init.c

mpi.o:		mpi.c mendel.h
	$(CC) -c mpi.c

sort.o:		sort.c
	$(CC) -c sort.c

ranlib.o:	ranlib.c ranlib.h
	$(CC) -c ranlib.c

diagnostics.o:	diagnostics.c mendel.h
	$(CC) -c diagnostics.c

selection.o:	selection.c
	$(CC) -c selection.c

offspring.o:	offspring.c
	$(CC) -c offspring.c

fileio.o:	fileio.c mendel.h
	$(CC) -c fileio.c

mem.o:		mem.c
	$(CC) -c mem.c
