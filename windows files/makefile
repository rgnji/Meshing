# Compiler and linker definitions
F77=mpiifx
CC=mpiicc
LD=mpiifx

# Flags
FFLAGS=/O2
CFLAGS=
OS_LIB=
FC_LIB=

# Libraries
LIBS=$(OS_LIB) $(FC_LIB)

#
CMDS   = xfdns.exe
FDNS_O = f1.obj f2.obj f3.obj f4.obj f5.obj f6.obj fdns.obj zone.obj io.obj flib.obj

#
all: $(CMDS)
xfdns.exe: $(FDNS_O)
    $(LD) /exe:xfdns.exe $(FDNS_O) $(LIBS)

# Fortran file compilation
.f.obj:
    $(F77) $(FFLAGS) /c $<

# C++ file compilation
.cc.obj:
    $(CC) $(CFLAGS) /c $<

# Clean up
clean:
    del /Q $(FDNS_O) xfdns.exe