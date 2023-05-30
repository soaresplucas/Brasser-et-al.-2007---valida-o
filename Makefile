export OPENGL=0
include ../../REBOUND_main/src/Makefile.defs

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I ../../REBOUND_main/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem2.c ../../forces/drag_forces2.c  -L. -lrebound $(LIB) -o rebound
	@echo ""
	@echo "REBOUND compiled successfully."

librebound:
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C ../../REBOUND_main/src/
	@-rm -f librebound.so
	@ln -s ../../REBOUND_main/src/librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C ../../REBOUND_main/src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
	@-rm -vfr rebound.dSYM/
	@echo "Removing output files..."
	@-rm -vf *out

full:
	make clean && rm -f *png && make && ./rebound && gnuplot plot2.plt
