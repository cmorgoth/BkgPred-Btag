CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

INC=$(shell pwd)

STDINCDIR :=-I$(INC)/include
STDLIBDIR := 

CPPFLAGS := $(shell root-config --cflags) $(STDINCDIR)
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR)

CPPFLAGS += -g

TARGET1 = FinalBkgPreBTagInc_CP
TARGET2 = FinalClosureBTag_CP

SRC1 = Main/FinalBkgPreBTagInc_CP.cc src/DM_Base.cc src/DM_1DRatio.cc
SRC2 = Main/ClosureBTag_CP.cc src/DM_Base.cc src/DM_1DRatio.cc

OBJ1 = $(SRC1:.cc=.o)
OBJ2 = $(SRC2:.cc=.o)


all : $(TARGET1) $(TARGET2)

$(TARGET1) : $(OBJ1)
	$(LD) $(CPPFLAGS) -o $(TARGET1) $(OBJ1) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET2) : $(OBJ2)
	$(LD) $(CPPFLAGS) -o $(TARGET2) $(OBJ2) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

%.o : %.cc	
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@	
	@echo $<
clean :
	rm -f *.o src/*.o $(TARGET1) *~

