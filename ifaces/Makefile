# Makefile for interfaces

# Directory
ODIR = obj
MMAODIR = Mobj
PYODIR = PYobj

# simple C object files
OBJS = test_driver.o
_OBJS = $(patsubst %,$(ODIR)/%,$(OBJS))

# rule for compiling simple C files
$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) 

# Mathlink C object files
MMAOBJS = xs.tm.o xs_mma.o 
_MMAOBJS = $(patsubst %,$(MMAODIR)/%,$(MMAOBJS))

# rules for compiling MathLink C files
$(MMAODIR)/%.tm.c: %.tm
	$(MPREP) -o $@ $< 
$(MMAODIR)/%.tm.o: $(MMAODIR)/%.tm.c $(PYODIR)/%.o
	$(CC) -c -o $@ $< -I$(MLINCDIR) $(CFLAGS) $(MLEXTRA)
$(MMAODIR)/%.o: %.c ifaces.h
	$(CC) -c -o $@ $< -I$(MLINCDIR) $(CFLAGS)


# Python embedding object files
PYOBJS = xs.o 
_PYOBJS = $(patsubst %,$(PYODIR)/%,$(PYOBJS))

# rule for compiling Python embedders
$(PYODIR)/%.o: %.c ifaces.h
	$(CC) -c -o $@ $< `python-config --includes` $(CFLAGS)


# targets

$(TESTTARGETS): %: $(_OBJS) $(_PYOBJS)
	$(CC) -o ../$@ $^ $(CFLAGS)  `python-config --libs`

# mma target

xs_mma.exe: $(_PYOBJS)  $(_MMAOBJS) 
	$(CC) -o ../$@/$(SYS)/$@ $^  -L$(MLLIBDIR) -l$(MLLIB) $(MLEXTRA) `python-config --libs`


.PHONY: clean
clean:
	-rm -f $(ODIR)/*.o $(PYODIR)/*.o $(MMAODIR)/*
	-rm -f ../test_driver
	-rm -f ../xs_mma.exe/$(SYS)/xs_mma.exe
