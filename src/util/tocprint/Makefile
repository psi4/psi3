
include ../../MakeVars

BIN = tocprint

CLIBS = -lpsio -lciomr -lipv1

MAN_FILE = tocprint.1

BINSRC = \
tocprint.c \

BINOBJ = $(BINSRC:%.c=%.o)

$(BIN): $(BINOBJ) $(CLIBS) $(FLIBS)
	$(CC) $(LDFLAGS) $^ -o $(BIN)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

%.d: %.c
	$(SHELL) -ec '$(CC) -M $(CFLAGS) $(INCLUDES) $< \
            | sed '\''s/$*\.o[ :]*/& $@ /g'\'' > $@'

install: $(BIN) $(MAN_FILE)
	$(INSTALL) $(INSTALL_BIN_FLAGS) $(BIN) $(INSTALL_BIN_DIR)
	$(INSTALL) $(INSTALL_FLAGS) $(MAN_FILE) $(INSTALL_MAN_DIR)


clean:
	/bin/rm -f $(BIN) *.o *.d core *~

ifeq ($(CC),gcc)
include $(BINSRC:%.c=%.d)
endif
