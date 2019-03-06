#include "ligra.h"
#include <unordered_map>
#include <queue>
#include "sweep.h"
#include "sweep_prob_first.h"
#include "ppr_based_sweep.h"
#include "sampleSort.h"
#include <set>
#include <limits>
using namespace std;

typedef pair<double, double> pairDouble;
typedef pair<uintE, double> pairIF;
struct pq_compare {
	bool operator ()(pairIF a, pairIF b) {
		if (a.second == b.second) {
			return a.first > b.first;
		}
		return a.second > b.second;
	}
};
typedef struct ppr {
	pairIF* A;
	long num_nonzeros;
	~ppr() {
		if (A != 0)
			free(A);
	}
};

vector<double> parse_alphas(const char *alphas) {
	vector<double> alphas_d;

	char *str = strdup(alphas);  // We own str's memory now.
	char *token;
	while ((token = strsep(&str, ";"))) {
		double candidate = atof(token);
		if (candidate > 0)
			alphas_d.push_back(candidate);
	}
	free(str);
	return alphas_d;
}

typedef struct PushStat {
	double max_pushes_per_node;
	double min_pushes_per_node;
	double sum_all_pushes;
	double number_points_at_least_one_push;

	double max_change_of_prob_mean;
	double min_change_of_prob_mean;
	double sum_change_of_prob_mean;
};

PushStat* initPushStat() {
	PushStat* toRet = (PushStat*) calloc(1, sizeof(PushStat));
	toRet->max_pushes_per_node = 0.;
	toRet->min_pushes_per_node = std::numeric_limits<double>::max();
	toRet->sum_all_pushes = 0.;
	toRet->number_points_at_least_one_push = 0;
	toRet->max_change_of_prob_mean = 0;
	toRet->min_change_of_prob_mean = std::numeric_limits<double>::max();
	toRet->sum_change_of_prob_mean = 0;
	return toRet;
}
void update_push_stat(PushStat* pushStats, double currentTotalPushes,
		double last_prob_changes_mean) {
	if (currentTotalPushes < 1) {
		return;
	}
	pushStats->sum_all_pushes = pushStats->sum_all_pushes + currentTotalPushes;
	pushStats->sum_change_of_prob_mean = pushStats->sum_change_of_prob_mean
			+ last_prob_changes_mean;
	if (pushStats->max_pushes_per_node < currentTotalPushes)
		pushStats->max_pushes_per_node = currentTotalPushes;
	if (pushStats->min_pushes_per_node > currentTotalPushes)
		pushStats->min_pushes_per_node = currentTotalPushes;
	if (pushStats->max_change_of_prob_mean < last_prob_changes_mean)
		pushStats->max_change_of_prob_mean = last_prob_changes_mean;
	if (pushStats->min_change_of_prob_mean > last_prob_changes_mean) {
		pushStats->min_change_of_prob_mean = last_prob_changes_mean;
	}
	pushStats->number_points_at_least_one_push =
			pushStats->number_points_at_least_one_push + 1;
}

void write_ppn_to_file(pairIF* A, uintE currentNode, FILE *result_ppn_file,
		long result_size, bool start_vertice_in_ppr) {
	fprintf(result_ppn_file, "\"%i\":[", currentNode);
	//ignore yourself
//	long start_from_neighbour = 0;
//	if (!start_vertice_in_ppr){
//		start_from_neighbour = 1;
//	}
	bool there_is_already_entry = false;
	for (long i = 0; i < result_size; i++) {
		if (!start_vertice_in_ppr) {
			if (A[i].first == currentNode) {
				continue;
			}
		}
		if (there_is_already_entry) {
			fprintf(result_ppn_file, ", ");
		}
		fprintf(result_ppn_file, "[%i, %.32g]", A[i].first, A[i].second);
		there_is_already_entry = true;
	}
	fprintf(result_ppn_file, "]");
}
template<class vertex>
double compute_conductance(graph<vertex>& GA, pairIF* p, long result_size) {
	unordered_set<uintE> S;
	long volS = 0;
	long edgesCrossing = 0;
	for (long i = 0; i < result_size; i++) {
		uintE v = p[i].first;
		S.insert(v);
		volS += GA.V[v].getOutDegree();
	}
	long denom = min(volS, GA.m - volS);
	for (long i = 0; i < result_size; i++) {
		uintE v = p[i].first;
		for (long j = 0; j < GA.V[v].getOutDegree(); j++) {
			uintE ngh = GA.V[v].getOutNeighbor(j);
			if (S.find(ngh) != S.end())
				edgesCrossing--;
			else
				edgesCrossing++;
		}
	}
	double conductance =
			(edgesCrossing == 0 || denom == 0) ?
					1 : (double) edgesCrossing / denom;
	return conductance;

}

timer t1;
bool to_use_prob_tolerance_to_stop(const double epsilon,
		const long max_number_iterations, const double prob_change_tolerance) {
	bool use_iter = false;
	if ((epsilon > 1 && max_number_iterations < 0)
			|| (epsilon < 1 && max_number_iterations > 0)) {
		fprintf(stderr, "define -e or -iter but not both\n");
		exit(1);
	} else if (max_number_iterations > 0) {
		if (prob_change_tolerance < 0) {
			printf("define -prob_change_tolerance");
			exit(1);
		}
		use_iter = true;
	} else if (epsilon < 1) {
		use_iter = false;
	} else {
		printf("bug selecting between iter and epsilon");
		exit(1);
	}

	return use_iter;
}

template<class vertex>
void compute_pprs(graph<vertex>& GA, commandLine P) {
	bool pr_output = P.getOptionValue("-verbose");
	if (pr_output) {
		cout << "number of vertices, graph = " << GA.n << endl;
		cout << "number of edges, graph = " << GA.m << endl;
	}
	t1.start();
	bool sweepcut = P.getOptionValue("-sweepcut");
	bool ppr_sweepcut = P.getOptionValue("-ppr_sweepcut");
	bool conductance_to_ppr_sweepcut = P.getOptionValue(
			"-conductance_to_ppr_sweepcut");
	const double best_percent = P.getOptionDoubleValue("-best_percent", 0);
	bool start_vertice_in_ppr = P.getOptionValue("-start_vertice_in_ppr");
	int only_one_option_allowed = 0;
	if (sweepcut) {
		only_one_option_allowed = only_one_option_allowed + 1;
	}
	if (ppr_sweepcut) {
		only_one_option_allowed = only_one_option_allowed + 1;
	}
	if (best_percent > 0) {
		only_one_option_allowed = only_one_option_allowed + 1;
	}
	if (conductance_to_ppr_sweepcut) {
		only_one_option_allowed = only_one_option_allowed + 1;
	}
	if (only_one_option_allowed > 1) {
		printf(
				"define one of -sweepcut or -pprsweepcut -best_percent  -conductance_to_ppr_sweepcut or nothing\n");
		exit(1);
	}

	const char* alphas = P.getOptionValue("-as");
	vector<double> alphas_d = parse_alphas(alphas);
	if (alphas_d.size() < 1) {
		printf("define -as! separated as string separated by ;");
		exit(1);
	} else if (alphas_d.size() > 1 && (only_one_option_allowed < 1)) {
		printf(
				"define how to select between different alphas for each node \n");
		exit(1);
	}
	if (pr_output) {
		printf("got following alphas: ");
		for (vector<double>::iterator it = alphas_d.begin();
				it != alphas_d.end(); ++it) {
			if (it != alphas_d.begin())
				printf(", ");
			printf("%.17g ", *it);
		}
		printf("\n");
	}

	const double epsilon = P.getOptionDoubleValue("-e", 10000);
	const long max_number_iterations = P.getOptionLongValue("-max_iter", -1);
	const double prob_change_tolerance = P.getOptionDoubleValue(
			"-prob_change_tolerance", 0.005);

	bool use_iter = to_use_prob_tolerance_to_stop(epsilon,
			max_number_iterations, prob_change_tolerance);
	int last_vertice_for_clustering = P.getOptionIntValue("-last", 5);
	int first_vertice_for_clustering = P.getOptionIntValue("-first", 0);
	PushStat* pstat = initPushStat();
	const intE n = GA.n;
	if (last_vertice_for_clustering < 0) {
		last_vertice_for_clustering = n;
	}

	const char *path_to_store_ppn = P.getOptionValue("-output");
	FILE *result_ppn_file;

	result_ppn_file = fopen(path_to_store_ppn, "w");
	fprintf(result_ppn_file, "{");
	if (pr_output) {
		printf("use iter: %d\n", use_iter);
		if (use_iter) {
			if (pr_output)
				printf("running with -prob_change_tolerance %17g,\n"
						"-max_iter %ld\n"
						"-last %d\n"
						"-sweepcut %d\n"
						"-ppr_sweepcut %d\n"
						"-best_percent %.4g\n"
						"-conductance_to_ppr_sweepcut %d\n",
						prob_change_tolerance, max_number_iterations,
						last_vertice_for_clustering, sweepcut, ppr_sweepcut,
						best_percent, conductance_to_ppr_sweepcut);
		} else {
			printf(
					"running with -e %.17g  -last %d, -sweepcut %d -ppr_sweepcut %d\n",
					epsilon, last_vertice_for_clustering, sweepcut,
					ppr_sweepcut);
		}
	}
	long clusters_size_sum = 0.;
	long cluster_size_min = std::numeric_limits<long>::max();
	long cluster_size_max = 0;
	double best_alpha = -1;

	for (long start = first_vertice_for_clustering;
			start < last_vertice_for_clustering; start++) {

		ppr* best_ppr = NULL;
		double best_smallest_score = std::numeric_limits<double>::max();
		long best_result_size = -1;

		for (vector<double>::iterator it = alphas_d.begin();
				it != alphas_d.end(); ++it) {
			double alpha = *it;
			long result_size;
			double score = 0;

			ppr* ppr_result = run_clustering(start, GA, alpha, use_iter,
					epsilon, max_number_iterations, pstat,
					prob_change_tolerance, start_vertice_in_ppr);

			if (sweepcut) {
				sweepObject sweep = sweepCut(GA, ppr_result->A,
						ppr_result->num_nonzeros, start);
				result_size = sweep.sizeS;
				score = sweep.conductance;
			} else if (ppr_sweepcut) {
				sweepObjectOrdererdPerProb sweep = sweepCutOrderedPerProb(GA,
						ppr_result->A, ppr_result->num_nonzeros, start);
				result_size = sweep.sizeS;
				score = sweep.conductance;
			} else if (conductance_to_ppr_sweepcut) {
				sweepObjectConductanceToProb ppr_sweep =
						sweepcut_conductance_to_prob(GA, ppr_result->A,
								ppr_result->num_nonzeros, start);
				result_size = ppr_sweep.sizeS;
				score = ppr_sweep.conductance_to_ppr;
			} else if (best_percent > 0) {
				double sum_not_zero_ppr = 0;
				for (long i = 1; i < ppr_result->num_nonzeros; i++) {
					sum_not_zero_ppr += ppr_result->A[i].second;
				}
				double sum_ppr_in_result = 0;
				result_size = 1;
				for (long i = 1; i < ppr_result->num_nonzeros; i++) {
					sum_ppr_in_result += ppr_result->A[i].second;
					if (sum_ppr_in_result
							>= (best_percent * sum_not_zero_ppr)) {
						result_size = i;
						result_size = result_size + 1;
						break;
					}
				}
				score = compute_conductance(GA, ppr_result->A, result_size);
			}
			//one alpha, classical version
			else {
				sampleSort(ppr_result->A, (uintE) ppr_result->num_nonzeros,
						pq_compare());
				result_size = ppr_result->num_nonzeros;
			}

			if (best_ppr == NULL || best_smallest_score > score) {
				if (best_ppr != NULL)
					free(best_ppr);
				best_alpha = alpha;
				best_ppr = ppr_result;
				best_smallest_score = score;
				best_result_size = result_size;
//				double sum_not_zero_ppr = 0;
//				double sum_ppr_in_result = 0;
//				for (long i = 1; i < best_ppr->num_nonzeros; i++) {
//					sum_not_zero_ppr += best_ppr->A[i].second;
//					if (i < best_result_size) {
//						sum_ppr_in_result += best_ppr->A[i].second;
//					}
//				}
//				printf(
//						"point %ld, alpha %.5g, better score %.17g, result size %ld, discarded %ld neigh, %.17g prob \n",
//						start, alpha, best_smallest_score, best_result_size,
//						ppr_result->num_nonzeros - result_size, (sum_not_zero_ppr-sum_ppr_in_result)/sum_not_zero_ppr);
			} else {
				free(ppr_result);
			}
		}
		double sum_not_zero_ppr = 0;
		double sum_ppr_in_result = 0;
		for (long i = 1; i < best_ppr->num_nonzeros; i++) {
			sum_not_zero_ppr += best_ppr->A[i].second;
			if (i < best_result_size) {
				sum_ppr_in_result += best_ppr->A[i].second;
			}
		}
//		printf(
//				"point %ld, alpha %.5g, better score %.17g, result size %ld, discarded %ld neigh, %.17g prob \n",
//				start, best_alpha, best_smallest_score, best_result_size,
//				best_ppr->num_nonzeros - best_result_size,
//				(sum_not_zero_ppr - sum_ppr_in_result) / sum_not_zero_ppr);

		clusters_size_sum = clusters_size_sum + best_result_size;
		if (best_result_size > cluster_size_max)
			cluster_size_max = best_result_size;

		if ((best_result_size > 0) && (best_result_size < cluster_size_min)) {
			cluster_size_min = best_result_size;
		}
		write_ppn_to_file(best_ppr->A, start, result_ppn_file, best_result_size,
				start_vertice_in_ppr);
		if (start < last_vertice_for_clustering - 1) {
			fprintf(result_ppn_file, ",\n");
		}
		free(best_ppr);

	}
	fprintf(result_ppn_file, "}");
	fclose(result_ppn_file);
	if (pr_output) {
		t1.reportTotal("computation time");
		cout << "mean cluster size: "
				<< (double) clusters_size_sum / last_vertice_for_clustering
				<< endl;
		cout << "min cluster size: " << cluster_size_min << endl;
		cout << "max cluster size: " << cluster_size_max << "\n\n" << endl;
		printf("min pushes = %.19g\nmax pushes= %.19g \nmean pushes = %.19g\n",
				pstat->min_pushes_per_node, pstat->max_pushes_per_node,
				pstat->sum_all_pushes / pstat->number_points_at_least_one_push);
		printf(
				"last min prob change = %.19g\nmax prob change= %.19g \nmean prob change= %.19g\n\n",
				pstat->min_change_of_prob_mean, pstat->max_change_of_prob_mean,
				pstat->sum_change_of_prob_mean
						/ pstat->number_points_at_least_one_push);
	}
	free(pstat);

}

//double compute_last_iter_prob_change_mean(long times_mass_transfered,
//		int last_iterations_to_consider, double *prob_change) {
//	if (times_mass_transfered < 1) {
//		return 0;
//	}
//	double last_prob_changes_mean = 0;
//	int last_iter_to_cons_to_compute_prob = last_iterations_to_consider;
//	if (times_mass_transfered < (long) last_iter_to_cons_to_compute_prob) {
//		last_iter_to_cons_to_compute_prob = (int) times_mass_transfered;
//	}
//	for (int i = 0; i < last_iter_to_cons_to_compute_prob; i++) {
//		last_prob_changes_mean = last_prob_changes_mean + prob_change[i];
//	}
//	last_prob_changes_mean = last_prob_changes_mean
//			/ last_iter_to_cons_to_compute_prob;
//	return last_prob_changes_mean;
//}

template<class vertex>
ppr* run_clustering(uintE start, graph<vertex>& GA, const double alpha,
		bool use_iter, const double epsilon, const long number_iterations,
		PushStat* pstat, double prob_change_tolerance,
		bool start_vertice_in_ppr) {
	unordered_map<uintE, pairDouble> pr;
	const double twoAlphaOverOnePlusAlpha = 2 * alpha / (1 + alpha);
	const double oneMinusAlphaOverOnePlusAlpha = (1 - alpha) / (1 + alpha);
	pr[start] = make_pair(0.0, 1.0); //starting vertex

	multiset<pairIF, pq_compare> q;
	pairIF start_with = make_pair(start, 1.0);
	unordered_map<uintE, pairIF> queue_mapping;

	q.insert(start_with);
	queue_mapping[start_with.first] = start_with;
	double totalPushes = 0;
	bool continue_to_push = true;
	double prob_change = 0;
	int current_index_prob_change = 0;
	double mass_transfered_to_ppr = 0;
	long times_mass_transfered = 0;
	while (continue_to_push) {
		if (use_iter && queue_mapping.size() != q.size()) {
			printf("bug queue mapping");
			exit(1);
		}
		totalPushes++;

		pairIF top = *q.begin();
		uintE v = top.first;
		q.erase(top);
		queue_mapping.erase(v);

		pairDouble v_pr = pr[v];
		uintE d = GA.V[v].getOutDegree();
		double add_to_prob = twoAlphaOverOnePlusAlpha * v_pr.second;
		//printf("added %d and %.60g \n",add_to_prob>0.0,add_to_prob);
		v_pr.first += add_to_prob;
		if (start_vertice_in_ppr || v != start) {
			times_mass_transfered = times_mass_transfered + 1;
			mass_transfered_to_ppr = mass_transfered_to_ppr + add_to_prob;
			prob_change = add_to_prob / mass_transfered_to_ppr;
		}

		double old_vr = v_pr.second;
		v_pr.second = 0;
		pr[v] = v_pr;
#ifdef WEIGHTED
		int sum_outgoing_weights = 0;
		for (long i = 0; i < d; i++) {
			uintE ngh = GA.V[v].getOutNeighbor(i);
			sum_outgoing_weights += GA.V[v].getOutWeight(i);
		}
#endif
		for (long i = 0; i < d; i++) {
			uintE ngh = GA.V[v].getOutNeighbor(i);
			if (ngh == v) {
				//ignore self-edges
				continue;
			}
			pairDouble ngh_pr = pr[ngh]; //creates default entry if non-existent in table
			double oldRes = ngh_pr.second;
#ifndef WEIGHTED
			ngh_pr.second += old_vr * oneMinusAlphaOverOnePlusAlpha / d;
#else
			double edge_fraction = (double)(GA.V[v].getOutWeight(i))/sum_outgoing_weights;
			ngh_pr.second += old_vr * oneMinusAlphaOverOnePlusAlpha*edge_fraction;
#endif

			pr[ngh] = ngh_pr;
			uintE ngh_d = GA.V[ngh].getOutDegree();
			if (use_iter) {
				if (queue_mapping.find(ngh) != queue_mapping.end()) {
					q.erase(queue_mapping[ngh]);
					queue_mapping.erase(ngh);
				}
				pairIF neigh_residual = make_pair(ngh, ngh_pr.second / ngh_d);
				q.insert(neigh_residual);
				queue_mapping[neigh_residual.first] = neigh_residual;

			} else {
				double epsd = epsilon * ngh_d;
				if (ngh_pr.second > epsd && oldRes < epsd) {
					q.insert(make_pair(ngh, ngh_pr.second / ngh_d));
				}
			}
		}

		if (use_iter) {
			if (q.empty()) {
				continue_to_push = false;
			} else if (totalPushes != 1 && current_index_prob_change == 0) {
				continue_to_push = prob_change_tolerance < prob_change;
			} else {
				continue_to_push = totalPushes < number_iterations;
			}
		} else {
			continue_to_push = !q.empty();
		}
	}

	update_push_stat(pstat, totalPushes, prob_change);
	t1.stop();
	pairIF* A = newA(pairIF, pr.size());

	long numNonzerosQ = 0;
	for (auto it = pr.begin(); it != pr.end(); it++) {
		if (it->second.first > 0) {
			A[numNonzerosQ++] = make_pair(it->first, it->second.first);
		}
	}
	ppr* to_return = (ppr*) malloc(sizeof(struct ppr));
	to_return->A = A;
	to_return->num_nonzeros = numNonzerosQ;
	return to_return;
}

