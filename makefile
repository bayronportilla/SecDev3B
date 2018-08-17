############################################################
#
# variables
#
############################################################

CC       = gcc
CFLAGS   = -I/home/bayron/include
CCFLAGS  = -L/home/bayron/lib
TARGET   = main 

#OBJS = FetchInfo.o ModQuad.o Params.o Units.o \
#allvars.o bulk.o interpol.o main.o matrix.o

OBJGS = main generic

#LIBS = -lgsl -lgslcblas -lm -lconfig


main:main.o generic.o matrix.o Units.o bulk.o Params.o FetchInfo.o interpol.o ModOct.o readlines.o allvars.o
	$(CC) main.o generic.o matrix.o Units.o bulk.o Params.o FetchInfo.o interpol.o ModOct.o readlines.o allvars.o $(CCFLAGS) -lgsl -lgslcblas -lm -lconfig -o main.x
	rm *.o






clean:
	rm *~
	rm *.x

run:
	./main.x


