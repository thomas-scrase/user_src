# A few file definitions
#-----------------------

# Define the sources
sources =  \
cell_membrane_potential_elements.cc \
monodomain_elements.cc \
refineable_monodomain_elements.cc \
bidomain_elements.cc \
refineable_bidomain_elements.cc

# Define the headers
headers =  \
cell_membrane_potential_elements.h \
monodomain_elements.h \
monodomain_flux_elements.h \
refineable_monodomain_elements.h \
bidomain_elements.h \
bidomain_flux_elements.h \
refineable_bidomain_elements.h \
monodomain_elements_strang_splitting.h \
monodomain_elements_toms_splitting_method.h \
monodomain_elements_integral_splitting_method.h \
monodomain_elements_crank_nicolson_splitting_method.h

# Define name of library
libname = cell_membrane_potential

# Combine headers and sources
headers_and_sources = $(headers) $(sources)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = libcell_membrane_potential.la

# Sources that the library depends on: 
#-------------------------------------
libcell_membrane_potential_la_SOURCES = $(headers_and_sources)

# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
include_HEADERS = $(libname).h 

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/cell_membrane_potential

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers) 


# Required libraries -- [assuming that we want to link against stuff in generic
#-------------------     add other oomph-lib libraries if you need them....]
# Include path for library headers -- need to refer to include files
# in their respective source directories as they will not have been 
# installed yet! 
AM_CPPFLAGS += -I$(top_builddir)/src/generic


# Combined header file
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time)
$(libname).h: dummy_$(libname).h

dummy_$(libname).h:  $(headers)
	echo  $(libname) $(headers) > all_$(libname).aux
	$(AWK) -f $(top_builddir)/bin/headers.awk < \
	       all_$(libname).aux > $(libname).h
	rm all_$(libname).aux


# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
if SYMBOLIC_LINKS_FOR_HEADERS
install-data-hook: 
	(cd $(library_includedir) && rm -f $(headers) )
	(echo "$(headers)" > include_files.list )
	($(top_builddir)/bin/change_headers_to_links.sh `pwd`)
	($(LN_S) `cat include_files.list.aux` $(library_includedir) )
	(rm -r include_files.list.aux include_files.list ) 
else
install-data-hook: 
endif


# Tidy up
#--------
clean-local:
	rm -f $(libname).h


