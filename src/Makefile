ARCH  = $(shell root-config --arch)
include $(ROOTSYS)/etc/Makefile.arch

#include $(HDF5ROOT)/config/Makefile.common
ROOTCFLAGS   := $(shell root-config --cflags)
SPMSDIR  := $(HOME)/legend/SPMSRoot
# HDF5 FLAGS
HDF5CFLAGS    = -I$(HDF5CPP)/include
HDF5LIBS      = -L$(HDF5CPP)/lib -lhdf5 -lhdf5_cpp


OBJDIR := $(SPMSDIR)/obj


OBJS  = 

EXES = spmsRoot 

LIBDIR 	:= $(shell pwd)/../lib
SDIRCR := $(shell pwd)


CXXFLAGS     := $(HDF5CFLAGS) $(ROOTCFLAGS) -I$(OBJDIR) -I.
##-I$(LOCAL) 


#LDLAGS := -g  -flto -Wl,-object_path_lto,spms2root.o
ifeq ($(ARCH),macosx64)
				CXXFLAGS += -fPIC
				CXX := clang++ -arch x86_64 -g -O0
				LD :=  clang++ -arch x86_64 -g -flto 
				SOFLAGS := -dynamiclib -single_module

else
				LIBLINK :=
				CXXFLAGS += -fPIC 
endif 


LIBS         := $(ROOTLIBS) $(HDF5LIBS) -L$(OBJDIR) -lSpmsAna
LIBS += $(FREE) 

SHLIBFILE    := lib$(LIB).$(DllSuf)

all: obj bin

obj: $(addsuffix .o,$(OBJS))

lib: $(LIBDIR)/$(SHLIBFILE)

bin: $(addsuffix .e,$(EXES))


%.o: %.cc
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $<

%.e: %.o 
	@echo "Linking $@  $(LDFLAGS) "
	$(LD)  $< $(LIBS) $(H2RLIBS) -o $@ 
	ln -sf $(SDIRCR)/$@ $(SPMSDIR)/$(subst .e,,$@)
	@echo "Linking done  "



$(LIBDIR)/$(SHLIBFILE): $(addsuffix .o,$(OBJS))
	@echo "Building shared library $@ ... $^ $(OBJS)"
	@mkdir -p $(LIBDIR)
	#ln -sf $(LIBDIR)/libhdf2root.dylib $(LIBDIR)/libhdf2root.so 
	@echo "...linking is Done do LD $(LD)  $(LDFLAGS)  $(SOFLAGS) "
	$(LD)   $(SOFLAGS) $^ -o $@ $(LIBS)
	@echo "...LD is Done $(LIBDIR)/$(SHLIBFILE) "


clean:
	@rm -f *.o
	@rm -f *.e
	@rm -f *~
	@rm -f $(LIBDIR)/$(SHLIBFILE)
	@for file in $(EXES); do (rm -f $$file.e); done
	@for file in $(EXES); do (rm -f $(EXESDIR)/$$file); done

# Rule to set up a symbolic links to the created shared library
linklib:
	@if [[ ( ( `root-config --platform` == "macosx" ) && \
			( ! -e lib/lib$(LIB).dylib ) ) || \
			( ! -e lib/lib$(LIB).so ) ]]; then \
		printf "Setting up soft links to the HDF2ROOT library ... "; \
		ln -sf $(SHLIBFILE) $(LIBDIR)/lib$(LIB).so; \
		if [[ `root-config --platform` == "macosx" ]]; then \
			ln -sf $(SHLIBFILE) $(LIBDIR)/lib$(LIB).dylib; \
		fi; \
		echo "Done"; \
	fi

show:	
	@echo "Library: $(LIBDIR)/$(SHLIBFILE)"
	@echo "LIBDIR =$(LIBDIR)"
	@echo "BINDIR =$(BINDIR)"
	@echo "SOFLAGS = $(SOFLAGS)"
	@echo "CURDIR  = $(CURDIR)"
	@echo "SHLIBFILE = $(SHLIBFILE)"
	@echo "FREE  = $(FREE)"
	@echo "CXX  = $(CXX)"
	@echo "LD  = $(LD)"
	@echo "H2RLIBS = $(H2RLIBS)"
	@echo "LIBS = $(LIBS)"
	@echo "INCDIR = $(INCDIR)"
	@echo "ARCH  = $(ARCH)"
	@echo "ROOTLIBS =$(ROOTLIBS)"
	@echo "HDF5LIBS =$(HDF5LIBS)"
	@echo "CXXFLAGS =$(CXXFLAGS)"
	@echo "LDFLAGS = $(LDFLAGS)"
	@echo "HDF5CFLAGS = $(HDF5CFLAGS)" 
	@echo "ROOTCFLAGS = $(ROOTCFLAGS)" 
	@echo "LDFLAGS = $(LDFLAGS)" 
	@echo "SOFLAGS = $(SOFLAGS)" 
	@echo "LIBDIR=   $(SHLIBFILE)"
	@echo "OBJDIR= $(OBJDIR)"
	@echo "SDIRCR = $(SDIRCR)"

.PHONY : clean show
#.PHONY : obj lib bin clean info linklib
