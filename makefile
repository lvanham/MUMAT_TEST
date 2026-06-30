# This is only important if MACHINE wasn't set
MACHINE?=viper

#GCC Macports
ifeq ($(MACHINE),macports)
	FC=mpif90
	COMP_OPT = -O2 -fexternal-blas -fallow-argument-mismatch
	LIBS = -L/usr/lib -L/opt/local/lib -lopenblas -lscalapack
	MPI_RUN  = mpiexec
	MPI_RUN_OPTS = -np 2
endif

#Intel Raven/Cobra/Viper
ifeq ($(MACHINE),$(filter $(MACHINE),cobra raven viper))
	FC = mpiifort
	COMP_OPT = -I${MKLROOT}/include/intel64/lp64 -I$(MKL_HOME)/include \
	            -O2 -traceback -assume noold_unit_star \
	            -march=skylake-avx512 -qopt-zmm-usage=high -fp-model strict -ip -fpp -DMPI_OPT
	COMP_OPT_DEBUG = -I${MKLROOT}/include/intel64/lp64 -I$(MKL_HOME)/include \
	           -g -O0 -traceback -check all -fpe0 -assume noold_unit_star \
	           -fp-model strict -fpp -DMPI_OPT
	LIBS = ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a \
	             ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \
	             ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a \
	             -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
	             ${MKLROOT}/lib/intel64/libmkl_sequential.a \
	             ${MKLROOT}/lib/intel64/libmkl_core.a \
	             ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
	             -Wl,--end-group -lpthread -lm -ldl
	MPI_RUN  = srun
	ifeq ($(MACHINE),raven)
		MPI_RUN_OPTS = --nodes=1 --ntasks-per-node=72 --time=0:30:00 -p express
	else ifeq ($(MACHINE),cobra)
		MPI_RUN_OPTS = --nodes=1 --ntasks-per-node=40 --time=0:30:00 -p express
	else ifeq ($(MACHINE),viper)
		MPI_RUN_OPTS = --nodes=1 --ntasks-per-node=128 --time=0:30:00 -p express
	endif
endif


# STELLOPT STUFF
LIBSTELL_INC = -I$(STELLOPT_PATH)/LIBSTELL/Release
LIBSTELL_LIB = $(STELLOPT_PATH)/LIBSTELL/Release/libstell.a
MUMAT_SRC = $(STELLOPT_PATH)/LIBSTELL/Sources/Modules/mumaterial_mod.f90
# FLAGS
FFLAGS= $(LIBSTELL_INC)
LDFLAGS=$(LIBSTELL_LIB)

OBJ=mumaterial_mod.o mumaterial_test.o
mumaterial_test.o : mumaterial_mod.o
EXE=xmumat_test

all: FFLAGS_OPT=$(COMP_OPT)
all: $(EXE)

debug: FFLAGS_OPT=$(COMP_OPT_DEBUG)
debug: $(EXE)
	@echo "Built in DEBUG mode."

mumaterial_mod.o: $(MUMAT_SRC)
	$(FC) -c -o $@ $< $(FFLAGS) $(FFLAGS_OPT)

%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(FFLAGS_OPT)

$(EXE): $(OBJ)
	$(FC) -o $@ $^ $(LDFLAGS) $(LIBS) $(FFLAGS_OPT)

run: $(EXE)
	$(MPI_RUN) $(MPI_RUN_OPTS) $(EXE) -mumat sphere_mu_2625.dat

clean:
	-rm *.o *.mod $(EXE)