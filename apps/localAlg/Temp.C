#include "ligra.h"
#include <unordered_map>
#include <set>
#include "sweep.h"
#include <string.h>
#include <vector>
using namespace std;

typedef pair<double, double> pairDouble;
typedef pair<uintE, double> pairIF;

struct pq_compare {
	bool operator ()(pairIF a, pairIF b) {
		return a.second > b.second;
	}
};



vector<double> parse_alphas(const char *alphas) {
	vector<double> alphas_d;

	char *str = strdup(alphas);  // We own str's memory now.
	char *token;
	while ((token = strsep(&str, ";"))) {
		alphas_d.push_back(atof(token));
	}
	free(str);
	return alphas_d;
}

template<class vertex>
void Compute(graph<vertex>& GA, commandLine P) {

		const char *alphas = "0.4;0.6";
		vector <double> alphas_d = parse_alphas(alphas);


		for (vector<double>::iterator it = alphas_d.begin(); it != alphas_d.end();
				++it) {
			printf("%.17g\n", *it);
		}

//	const char *alphas = P.getOptionValue("-alphas");
//	vector <double> alphas_d = parse_alphas(alphas);
//
//
//	for (vector<double>::iterator it = alphas_d.begin(); it != alphas_d.end();
//			++it) {
//		printf("%.17g\n", *it);
//	}

//	multiset<pairIF,pq_compare> q;
//	pairIF a = make_pair(1,1.0);
//	q.insert(a);
//	q.insert(make_pair(2,2.0));
//	q.insert(make_pair(3,3.0));
//	for (multiset<pairIF,pq_compare>::iterator it=q.begin(); it!=q.end(); ++it){
//		//pairIF a = *it;
//		printf("(%d, %f)", it->first, it->second);
//		printf("\n\n\n");
//	}

//
//
//	q.erase(uint_pairs_mapping[1]);
//	cout << "nums contains " << uint_pairs_mapping.size() << " elements.\n";
//
//	unordered_map <uintE,pairDouble>::iterator it =uint_pairs_mapping.find(1);
//	printf("%d\n",it!=uint_pairs_mapping.end());

	//uint_pairs_mapping.erase(1);


//	cout << "nums contains " << uint_pairs_mapping.size() << " elements.\n";

}
