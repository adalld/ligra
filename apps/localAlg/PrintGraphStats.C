#include "ligra.h"

template<class vertex>
void Compute(graph<vertex>& GA, commandLine P){
	printf("nodes:%lu\n", GA.n);
	printf("edges:%lu\n", GA.m);
	const int edgeToPrintFrom = P.getOptionIntValue("-node", 0);

    uintE d = GA.V[edgeToPrintFrom].getOutDegree();
    for(long i=0;i<d;i++) {
      uintE ngh = GA.V[edgeToPrintFrom].getInNeighbor(i);
      printf("%iu->%iu\n", edgeToPrintFrom, ngh);
    }
}
