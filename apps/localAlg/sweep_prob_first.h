// This code is part of the project "Parallel Local Clustering
// Algorithms".  Copyright (c) 2016 Julian Shun
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights (to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// This file contains both the serial and the parallel implementations
// of sweep cut. Currently only works with uncompressed graphs, and
// not with compressed graphs.

// To print out the set found by the sweep cut, uncomment the
// following line (not to be used when measuring running time).
//#define PRINT_SET

#include "graph.h"
#include <unordered_map>
#include <unordered_set>
#include <math.h>
#include <float.h>
#include "sparseSet.h"
#include "sampleSort.h"
#include "blockRadixSort.h"
using namespace std;

typedef pair<double, double> pairDouble;
typedef pair<uintE, double> pairIF;

template<class vertex>

struct sweep_compare_prob_first {
	vertex* V;
	sweep_compare_prob_first(vertex* _V) :
			V(_V) {
	}
	bool operator ()(pairIF a, pairIF b) {
		return a.second > b.second;
	}
};

struct sweepObjectConductanceToProb {
	double conductance;
	double conductance_to_ppr;
	long sizeS, volS, vol, edgesCrossing;
	sweepObjectConductanceToProb(double _c, double _ctoppr, long _s, long _volS,
			long _vol, long _e) :
			conductance(_c), conductance_to_ppr(_ctoppr), sizeS(_s), volS(
					_volS), vol(_vol), edgesCrossing(_e) {
	}
};

template<class vertex, class fType>
sweepObjectConductanceToProb sweepcut_conductance_to_prob(graph<vertex>& GA,
		pair<uintE, fType> *p, long numNonzeros, long start) {
	double sum_not_zero_ppr = 0;
	for (long i = 1; i < numNonzeros; i++) {
		sum_not_zero_ppr += p[i].second;
	}

	int procs = getWorkers();
	setWorkers(1);
	s1.start();
	sampleSort(p, (uintE) numNonzeros, sweep_compare_prob_first<vertex>(GA.V));
	//find cluster with sweep cut
	unordered_set<uintE> S;
	long volS = 0;
	long edgesCrossing = 0;
	double bestConductance = DBL_MAX;
	double best_conductance_to_ppr = DBL_MAX;
	double best_sum_already_seen_prob = 0;
	long bestCut = 0;
	long bestEdgesCrossing = -1;
	long bestVol = -1;
	double sum_already_seen_ppr = 0;
	for (long i = 0; i < numNonzeros; i++) {
		uintE v = p[i].first;
		S.insert(v);
		volS += GA.V[v].getOutDegree();
		long denom = min(volS, GA.m - volS);
		for (long j = 0; j < GA.V[v].getOutDegree(); j++) {
			uintE ngh = GA.V[v].getOutNeighbor(j);
			if (S.find(ngh) != S.end())
				edgesCrossing--;
			else
				edgesCrossing++;
		}
		//ignore first
		if(i==0)
			continue;
		sum_already_seen_ppr += p[i].second;
		double conductance =
				(edgesCrossing == 0 || denom == 0) ?
						1 : (double) edgesCrossing / denom;
		double conductance_to_ppr =
				(conductance == 0 || sum_already_seen_ppr == 0) ?
						1 : conductance / (sum_already_seen_ppr/sum_not_zero_ppr);

		if (conductance_to_ppr < best_conductance_to_ppr) {
			//printf("\nold cond: %.17g, new cond: %.17g, old cond to ppr  %.17g, new cond to ppr %.17g\n", bestConductance, conductance,best_conductance_to_ppr,conductance_to_ppr );
			best_sum_already_seen_prob = sum_already_seen_ppr;
			bestConductance = conductance;
			best_conductance_to_ppr = conductance_to_ppr;
			bestCut = i;
			bestEdgesCrossing = edgesCrossing;
			bestVol = volS;
		}
	}
#ifdef PRINT_SET
	cout << "S = {";
	{	for(long i=0;i<bestCut;i++) {
			cout << p[i].first << ", ";
		}
	}
	cout << p[bestCut].first << "}" << endl;
#endif
	setWorkers(procs);
	return sweepObjectConductanceToProb(bestConductance, best_conductance_to_ppr,
			bestCut + 1, bestVol, volS, bestEdgesCrossing);
}
