########################################################################
# General Makefile
########################################################################
PROGRAM = cg_test
OBJS = cg_test.o
DEPS = cg.c cg.h

ODIR = .
SDIR = .

CC=gcc
CXX=g++
FC=gfortran

LDFLAGS += -L/opt/local/lib
CFLAGS += -I/opt/local/include
CFLAGS += -I.

CFLAGS += -Wall
CFLAGS += -ggdb
CFLAGS += -O0
#CFLAGS += -O3

#LIBS += -lgsl
#LIBS += -lgslcblas
#LIBS += -llapack -lptcblas -latlas
#LIBS += -lnlopt
LIBS += -lm

.PHONEY: clean

$(PROGRAM): $(patsubst %, $(ODIR)/%, $(OBJS))
	$(CC) $(CFLAGS) $(LDFLAGS) $(LIBS) $^ -o $@

$(ODIR)/%.o: $(SDIR)/%.c $(patsubst %, $(SDIR)/%, $(DEPS))
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm $(ODIR)/*.o $(PROGRAM)
