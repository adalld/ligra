ifdef LONG
INTT = -DLONG
endif

ifdef EDGELONG
INTE = -DEDGELONG
endif

#compilers
ifdef CILK
PCC = g++
PCFLAGS = -fcilkplus -lcilkrts -std=c++11 -O2 -DCILK $(INTT) $(INTE)
PLFLAGS = -fcilkplus -lcilkrts

else ifdef MKLROOT
PCC = icpc
PCFLAGS = -O3 -std=c++11 -DCILKP $(INTT) $(INTE)

else ifdef OPENMP
PCC = g++
PCFLAGS = -fopenmp -O3 -std=c++11 -DOPENMP $(INTT) $(INTE)

else
PCC = g++
PCFLAGS = -std=c++11 -O2 $(INTT) $(INTE)
endif

COMMON= ligra.h graph.h compressedVertex.h vertex.h utils.h IO.h parallel.h gettime.h quickSort.h parseCommandLine.h byte.h byteRLE.h nibble.h byte-pd.h byteRLE-pd.h nibble-pd.h vertexSubset.h encoder.C
LOCAL_COMMON = sweep.h sparseSet.h sampleSort.h ppr_based_sweep.h sweep_prob_first.h ACL-Serial-Opt-Naive.h
INTSORT= blockRadixSort.h transpose.h
SERIAL = ACL-Serial ACL-Serial-Naive ACL-Serial-Opt ACL-Serial-Opt-Naive HeatKernel-Serial HeatKernel-Randomized-Serial Nibble-Serial PrintGraphStats ACL-Serial-Opt-Naive-Weighted
PARALLEL = ACL-Sync-Local ACL-Sync-Local-Opt HeatKernel-Parallel HeatKernel-Randomized-Parallel Nibble-Parallel

all: $(SERIAL) $(PARALLEL)

ser: $(SERIAL)

par: $(PARALLEL)

% : %.C $(COMMON) $(LOCAL_COMMON) $(INTSORT)
	$(PCC) $(PCFLAGS) -o $@ $< 

$(COMMON):
	ln -s ../../ligra/$@ .

$(INTSORT):
	ln -s ../../utils/$@ .

.PHONY : clean

clean :
	rm -f *.o $(SERIAL) $(PARALLEL)

cleansrc :
	rm -f *.o $(SERIAL) $(PARALLEL)
	rm $(COMMON) $(INTSORT)
