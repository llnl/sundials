# ---------------------------------------------------------------
# Programmer: Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Makefile to generate SUNDIALS fortran interfaces with swig
# ---------------------------------------------------------------

SWIG ?= swig

ARKODE=farkode_mod farkode_arkstep_mod farkode_erkstep_mod farkode_sprkstep_mod farkode_mristep_mod farkode_lsrkstep_mod farkode_splittingstep_mod farkode_forcingstep_mod
CVODE=fcvode_mod
CVODES=fcvodes_mod
IDA=fida_mod
IDAS=fidas_mod
KINSOL=fkinsol_mod
CORE=fsundials_core_mod
NVECTOR=openmp pthreads serial parallel manyvector mpiplusx
SUNMATRIX=band dense sparse
SUNLINSOL=band dense lapackdense klu spbcgs spfgmr spgmr sptfqmr pcg
SUNNONLINSOL=newton fixedpoint
SUNADAPTCONTROLLER=imexgus soderlind mrihtol
SUNADJOINTCHECKPOINTSCHEME=fixed
SUNDOMEIGEST=power arnoldi

INCLUDES=-I../include

.PHONY: .SETINT64 .SETINT32 all all32 all64 modules core arkode cvode cvodes ida idas kinsol nvector mpimanyvector sunmatrix sunlinsol sunnonlinsol sunadaptcontroller sunadjointcheckpointscheme sundomeigest clean

.SETINT32:
	$(eval INT_SIZE=32)

.SETINT64:
	$(eval INT_SIZE=64)

all: all32 all64

all32: .SETINT32 modules

all64: .SETINT64 modules

modules: core arkode cvode cvodes ida idas kinsol nvector sunmatrix sunlinsol sunnonlinsol sunadaptcontroller sunadjointcheckpointscheme sundomeigest

core: $(CORE:%:sundials/%.i)
	@for i in ${CORE} ; do \
		set -x; \
		${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sundials/fmod_int${INT_SIZE} -o ../src/sundials/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} sundials/$${i}.i; \
	done

arkode: $(ARKODE:%:arkode/%.i)
	@for i in ${ARKODE} ; do \
		set -x; \
		${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/arkode/fmod_int${INT_SIZE} -o ../src/arkode/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} arkode/$${i}.i; \
	done

cvode: $(CVODE:%:cvode/%.i)
	@for i in ${CVODE} ; do \
		set -x; \
		${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/cvode/fmod_int${INT_SIZE} -o ../src/cvode/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} cvode/$${i}.i; \
	done

cvodes: $(CVODES:%:cvodes/%.i)
	@for i in ${CVODES} ; do \
		set -x; \
		${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/cvodes/fmod_int${INT_SIZE} -o ../src/cvodes/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} cvodes/$${i}.i; \
	done

ida: $(IDA:%:ida/%.i)
	@for i in ${IDA} ; do \
		set -x; \
		${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/ida/fmod_int${INT_SIZE} -o ../src/ida/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} ida/$${i}.i; \
	done

idas: $(IDAS:%:idas/%.i)
	@for i in ${IDAS} ; do \
		set -x; \
		${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/idas/fmod_int${INT_SIZE} -o ../src/idas/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} idas/$${i}.i; \
	done

kinsol: $(KINSOL:%:kinsol/%.i)
	@for i in ${KINSOL} ; do \
		set -x; \
		${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/kinsol/fmod_int${INT_SIZE} -o ../src/kinsol/fmod_int${INT_SIZE}/$${i}.c ${INCLUDES} kinsol/$${i}.i; \
	done

nvector: $(NVECTOR:%:nvector/fnvector_%_mod.i) mpimanyvector
	@for i in ${NVECTOR} ; do \
		set -x; \
		${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/nvector/$${i}/fmod_int${INT_SIZE} -o ../src/nvector/$${i}/fmod_int${INT_SIZE}/fnvector_$${i}_mod.c ${INCLUDES} nvector/fnvector_$${i}_mod.i; \
	done

mpimanyvector: nvector/fnvector_mpimanyvector_mod.i
	set -x; \
	${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/nvector/manyvector/fmod_int${INT_SIZE} -o ../src/nvector/manyvector/fmod_int${INT_SIZE}/fnvector_mpimanyvector_mod.c ${INCLUDES} nvector/fnvector_mpimanyvector_mod.i;

sunmatrix: $(SUNMATRIX:%:sunmatrix/fsunmatrix_%_mod.i)
	@for i in ${SUNMATRIX} ; do \
		set -x; \
		${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sunmatrix/$${i}/fmod_int${INT_SIZE} -o ../src/sunmatrix/$${i}/fmod_int${INT_SIZE}/fsunmatrix_$${i}_mod.c ${INCLUDES} sunmatrix/fsunmatrix_$${i}_mod.i; \
	done

sunlinsol: $(SUNLINSOL:%:sunlinsol/fsunlinsol_%_mod.i)
	@for i in ${SUNLINSOL} ; do \
		set -x; \
		${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sunlinsol/$${i}/fmod_int${INT_SIZE} -o ../src/sunlinsol/$${i}/fmod_int${INT_SIZE}/fsunlinsol_$${i}_mod.c ${INCLUDES} sunlinsol/fsunlinsol_$${i}_mod.i; \
	done

sunnonlinsol: $(SUNNONLINSOL:%:sunnonlinsol/fsunnonlinsol_%_mod.i)
	@for i in ${SUNNONLINSOL} ; do \
		set -x; \
		${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sunnonlinsol/$${i}/fmod_int${INT_SIZE} -o ../src/sunnonlinsol/$${i}/fmod_int${INT_SIZE}/fsunnonlinsol_$${i}_mod.c ${INCLUDES} sunnonlinsol/fsunnonlinsol_$${i}_mod.i; \
	done

sunadaptcontroller: $(SUNADAPTCONTROLLER:%:sunadaptcontroller/fsunadaptcontroller_%_mod.i)
	@for i in ${SUNADAPTCONTROLLER} ; do \
		set -x; \
		${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sunadaptcontroller/$${i}/fmod_int${INT_SIZE} -o ../src/sunadaptcontroller/$${i}/fmod_int${INT_SIZE}/fsunadaptcontroller_$${i}_mod.c ${INCLUDES} sunadaptcontroller/fsunadaptcontroller_$${i}_mod.i; \
	done

sunadjointcheckpointscheme: $(SUNADJOINTCHECKPOINTSCHEME:%:sunadjointcheckpointscheme/fsunadjointcheckpointscheme_%_mod.i)
	@for i in ${SUNADJOINTCHECKPOINTSCHEME} ; do \
		set -x; \
		${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sunadjointcheckpointscheme/$${i}/fmod_int${INT_SIZE} -o ../src/sunadjointcheckpointscheme/$${i}/fmod_int${INT_SIZE}/fsunadjointcheckpointscheme_$${i}_mod.c ${INCLUDES} sunadjointcheckpointscheme/fsunadjointcheckpointscheme_$${i}_mod.i; \
	done

sundomeigest: $(SUNDOMEIGEST:%:sundomeigest/fsundomeigest_%_mod.i)
	@for i in ${SUNDOMEIGEST} ; do \
		set -x; \
		${SWIG} -DGENERATE_INT${INT_SIZE} -fortran -outdir ../src/sundomeigest/$${i}/fmod_int${INT_SIZE} -o ../src/sundomeigest/$${i}/fmod_int${INT_SIZE}/fsundomeigest_$${i}_mod.c ${INCLUDES} sundomeigest/fsundomeigest_$${i}_mod.i; \
	done

clean:
	rm ../src/**/fmod*/*.c; rm ../src/**/fmod*/*.f90
