#define WEIGHTED 1
#include "ACL-Serial-Opt-Naive.h"

template<class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
	compute_pprs(GA, P);
}
