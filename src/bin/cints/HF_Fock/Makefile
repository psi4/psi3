include ../../../MakeVars

BIN = cints
MAN = cints.1

#Additional options for xlc compiler
ifeq ($(CC),xlc)
CFLAGS += -D__MATH__ -D__STR__
endif

#CFLAGS = -g

INCLUDES += -I../ -I../Tools/

BINSRC = fock.c hash.c
BINOBJ = $(BINSRC:%.c=%.o)

default_ints:   $(BINOBJ)
	echo "Done with HF_Fock"

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	/bin/rm -f $(BIN) *.d *.o *.i *.s core *~
