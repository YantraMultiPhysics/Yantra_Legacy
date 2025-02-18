# -*- makefile -*-
COMPILER = GCC
TEMPLATE = make.inc

include $(TEMPLATE)
.PHONY: clean,clean_all

$(MODULE).so: $(MODULE).pyf
	f2py2 -c $(MODULE).pyf --fcompiler=$(FC) --f90flags=$(F90FLAGS) --opt=$(OPTFLAGS) $(XF2PYFLAGS2) $(MODULE).f90  
$(MODULE).pyf: $(MODULE).f90
	f2py2 -h $(MODULE).pyf -m $(MODULE) $(XF2PYFLAGS1) $(MODULE).f90

clean:
	rm -rf $(MODULE).pyf $(MODULE).so

clean_all :
	rm -rf *.pyf *.so
