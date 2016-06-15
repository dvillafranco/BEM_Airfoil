# Make file for BEM2d2.f90 and external modules.  Code and makefile 
#	developed by Trevor Wood, Boston University,  2000.  
#	under the direction of Professor Sheryl Grace
#
#
# Use gmake (GNU make-3.75)
#
# Compilation Options
#
#    MP:       Include OpenMP multi-processing directives
#    DEBUG:    Include debugging info
#    NODEBUG:  No debugging info, Optimized for speed
#    MAIN:     Main program to compile to an executable
#    OPTFLAG:  Use -O1 with OMP, -Ofast=ip27 for serial compilation

#MP         = mp
#OPTFLAG    = -O1
#OPTFLAG    = -Ofast=ip27
#DBGFLAG    = -g #-DEBUG:trap_uninitialized=ON
                #-DEBUG:subscript_check=ON:verbose_runtime=ON

#DEBUG      = deb
NODEBUG    = rel

MAIN       = BEM2d2

.SUFFIXES:	.o .f90 .mod

#BINDIR     = ../bin/
BINDIR     = ./dat/
INCDIR     = 
LIBDIR     = 
LIBS       = #-lscs$(MPSUF) #-lcomplib.sgimath$(MPSUF)

MPFLAG     = $(MP:mp=-mp) $(MP:mp=-mp_schedtype=INTERLEAVE) $(MP:mp=-chunk=1)
MPSUF      = $(MP:mp=_mp)

RELSUF     = $(DEBUG:deb=_deb)$(NODEBUG:rel=_rel)
RELFLAG    = $(DEBUG:deb=$(DBGFLAG)) $(NODEBUG:rel=$(OPTFLAG))
#RELFLAG    = $(DEBUG:deb=-g -check_bounds -DEBUG:verbose_runtime=ON) \
#             $(NODEBUG:rel=-Ofast=ip27 -OPT:roundoff=0)

FFLAGS     = $(RELFLAG) $(MPFLAG) #-64 #-64 -r10000 -mips4
F90FLAGS   = $(FFLAGS)
F90        = pgf90

SRC = TypeKinds.f90 BEM2d_typ.f90 Global2.f90 IO.f90 Geometry.f90 Greens.f90 \
      Vortex.f90 BEM2d_utils.f90 Farfield.f90
      
OBJ = $(SRC:.f90=.o) 

default:
#	@$(MAKE) getmod
	@$(MAKE) $(MAIN)
#	@$(MAKE) storemod

.f90: $(OBJ)
	@echo ... Final compilation of $@$(RELSUF)$(MPSUF)
	$(F90) $(INCDIR) $(F90FLAGS) -o $(BINDIR)$@$(RELSUF)$(MPSUF) $(OBJ) $< $(LIBDIR) $(LIBS)

.f90.o:
	$(F90) $(INCDIR) $(FFLAGS) -c $<

mod: 
	$(MAKE) $(OBJ)

clean:
	rm -f *.o
	rm -f *.mod
          
#storemod:
#	@echo ... Moving modules  to  obj/$(NODEBUG)$(DEBUG)$(MPSUF)/
#	@mv -f *.o *.mod obj/$(NODEBUG)$(DEBUG)$(MPSUF)

#getmod:
#	@echo ... Moving modules from obj/$(NODEBUG)$(DEBUG)$(MPSUF)/
#	@mv -f obj/$(NODEBUG)$(DEBUG)$(MPSUF)/* .

deb2rel:
#	$(MAKE) storemod
	cp -f Makefile.rel Makefile
#	$(MAKE) getmod

rel2deb:
#	$(MAKE) storemod
	cp -f Makefile.deb Makefile
#	$(MAKE) getmod

zip: $(SRC)
	zip kirchhoff $(SRC) Makefile main.f90 ../mat/*.m
	@echo 'kirchhoff.zip generated'
	cp -f kirchhoff.zip ${HOME}/archive/
	@echo 'backup copy stored in $(HOME)/archive'


