# makefile eventually
TOP	:=	$(shell pwd)

EXEC = rwfd
EXECDIR = $(TOP)/exe 

COMPILER = g++

INCS = -I$(TOP)/src/main
LIBS = -L$(TOP)/lib 
libsys = -lm 


OFLAGS = -O2
DFLAGS = 
WARN = -Wall


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set basic paths
libpaths = src/math
libpaths += src/physics
libpaths += src/utilities

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set paths

# projects = 
# include MyConfig
# libpaths += $(projects)

libpaths += src/main 

libdirs = $(dir $(libpaths))

libnames := $(notdir $(libpaths))

liblist := $(foreach libname,$(libnames),-l$(libname))

liblist := $(subst -lmain,,$(liblist))

LIBS += $(liblist) $(liblist)
LIBS += $(libsys)

libincludes := $(foreach libpath,$(libpaths),\
		$(libpath)/rwfd_$(notdir $(libpath)).h)

autoinclude = src/main/rwfd_automatic_include.h
autotext = \/\* automatically generated from MyConfig \*\/

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set compiler flags
CFLAGS = $(DFLAGS) $(OFLAGS) $(INCS) $(WARN)
export

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set make targes
rwfd: $(autoinclude) 
	for X in $(libnames); do mkdir -p lib/obj/$$X; done
	for X in $(libpaths); do $(MAKE) -C $$X; done

$(autoinclude): 
	echo $(autotext) > $(autoinclude)
	for X in $(libincludes); do \
		echo \#include \"$(TOP)/$$X\" >> $(autoinclude); \
	done 

clean:
	rm -r lib 
	rm $(autoinclude)
	rm $(autoinitial)