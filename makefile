# swig binary to use... important if using a local install
SWIG   = swig
CC     = g++
CFLAGS = -pipe -O3
CFLAGS1 = -pipe -O3
# linker used by the OSX...
LINKER = ld -bundle -flat_namespace -undefined suppress 

# In this folder the numpy include should also be there... if not there, 
# create a symblinc to the numpy include (hint: look for the numpy egg)
PYTHON_INCLUDE = -I /usr/include/python2.5

ISRC           = rbf_solver
SRC       = ../../C/rbf/$(ISRC).c
INTERFACE = $(ISRC).i
WRAP      = $(ISRC)_wrap.c
TARGET    = _$(ISRC).so
OBJS      = $(ISRC).o $(ISRC)_wrap.o

ISRC1           = vorticity_evaluation
SRC1       = ../../C/rbf/$(ISRC1).c
INTERFACE1 = $(ISRC1).i
WRAP1      = $(ISRC1)_wrap.c
TARGET1    = _$(ISRC1).so
OBJS1      = $(ISRC1).o $(ISRC1)_wrap.o

ISRC2           = rbf_solver
SRC2       = ../../C++/rbf/$(ISRC1).cxx
INTERFACE2 = $(ISRC1).i
WRAP2      = $(ISRC1)_wrap.cxx
TARGET2    = _$(ISRC1).so
OBJS2      = $(ISRC1).o $(ISRC1)_wrap.o

ISRC3           = vorticity_evaluation
SRC3       = ../../C++/rbf/$(ISRC1).cxx
INTERFACE3 = $(ISRC1).i
WRAP3      = $(ISRC1)_wrap.cxx
TARGET3    = _$(ISRC1).so
OBJS3      = $(ISRC1).o $(ISRC1)_wrap.o

all:
	make rbfc
	make vorc
	make rbfcxx
	make vorcxx
	make fmm
rbfc: 
	$(RM) $(TARGET) $(ISRC).py
	$(SWIG) -python $(INTERFACE)
	$(CC) -c $(SRC) $(WRAP) $(CFLAGS) $(PYTHON_INCLUDE)
	$(LINKER) $(OBJS) -o $(TARGET)
	$(RM) $(WRAP) $(OBJS)

vorc: 
	$(RM) $(TARGET1) $(ISRC1).py
	$(SWIG) -python $(INTERFACE1)
	$(CC) -c $(SRC1) $(WRAP1) $(CFLAGS) $(PYTHON_INCLUDE)
	$(LINKER) $(OBJS1) -o $(TARGET1)
	$(RM) $(WRAP1) $(OBJS1)

rbfcxx:
	$(RM) $(TARGET2) $(ISRC2).py
	$(SWIG) -c++ -python $(INTERFACE2)
	$(CC) -c $(SRC2) $(WRAP2) $(CFLAGS1) $(PYTHON_INCLUDE)
	$(LINKER) $(OBJS2) -o $(TARGET2)
	$(RM) $(WRAP2) $(OBJS2)

vorcxx:
	$(RM) $(TARGET3) $(ISRC3).py
	$(SWIG) -c++ -python $(INTERFACE3)
	$(CC) -c $(SRC3) $(WRAP3) $(CFLAGS1) $(PYTHON_INCLUDE)
	$(LINKER) $(OBJS3) -o $(TARGET3)
	$(RM) $(WRAP3) $(OBJS3)

fmm:
	f2py -c -m fmm fmm.f
clean:
	$(RM) $(TARGET) $(TARGET1) $(ISRC) $(ISRC1).py *.pyc fmm.so
