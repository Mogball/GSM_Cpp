#define EIGEN_NO_DEBUG

#include "GSM_Solver.h"
#include <Eigen>
#include <queue>

typedef Eigen::SparseMatrix<double> sparse_matrix;
typedef Eigen::VectorXd eigen_vector;
typedef Eigen::LeastSquaresConjugateGradient<sparse_matrix> linear_solver;
typedef std::queue<__int32> int_queue;

vector *generate_grid_nodes(__int32 m, __int32 n, __float32 d) {
	vector *nodes = new vector[m * n];
	__int32 node_index = 0;
	for (__int32 y = 0; y < n; y++) {
		for (__int32 x = 0; x < m; x++) {
			nodes[node_index++] = vector(x * d, y * d);
		}
	}
	return nodes;
}

bool **generate_strut_matrix(__int32 m, __int32 n) {
	__int32 num_nodes = m * n;
	bool **strut_matrix = new bool *[num_nodes];
	for (__int32 k = 0; k < num_nodes; k++) {
		strut_matrix[k] = new bool[num_nodes];
	}
	for (__int32 i = 0; i < num_nodes; i++) {
		for (__int32 j = 0; j < num_nodes; j++) {
			strut_matrix[i][j] = i != j;
		}
	}
	return strut_matrix;
}

__int32 **generate_strut_indices(__int32 m, __int32 n) {
	__int32 num_nodes = m * n;
	__int32 strut_index = 0;
	__int32 **strut_indices = new __int32 *[num_nodes];
	for (__int32 k = 0; k < num_nodes; k++) {
		strut_indices[k] = new __int32[num_nodes];
	}
	for (__int32 i = 0; i < num_nodes; i++) {
		for (__int32 j = 0; j < num_nodes; j++) {
			if (i == j) {
				strut_indices[i][j] = -1;
			} else if (i > j) {
				strut_indices[i][j] = strut_indices[j][i];
			} else {
				strut_indices[i][j] = strut_index++;
			}
		}
	}
	return strut_indices;
}

__int32 **generate_strut_inverse(__int32 m, __int32 n) {
	__int32 num_nodes = m * n;
	__int32 num_struts = num_nodes * (num_nodes - 1) / 2;
	__int32 strut_index = 0;
	__int32 **strut_inverse = new __int32 *[num_struts];
	for (__int32 k = 0; k < num_struts; k++) {
		strut_inverse[k] = new __int32[2];
	}
	for (__int32 i = 0; i < num_nodes; i++) {
		for (__int32 j = i + 1; j < num_nodes; j++) {
			strut_inverse[strut_index][0] = i;
			strut_inverse[strut_index][1] = j;
			strut_index++;
		}
	}
	return strut_inverse;
}

__float32 *get_loads(__int32 m, __int32 n, __int32 num_loads, vector *nodes,
                     tuple<__int32, __int32> support_reactions,
                     tuple<__int32, __float32> *applied_loads) {
	__int32 num_nodes = m * n;
	__float32 *loads = new __float32[num_nodes];
	for (__int32 k = 0; k < num_loads; k++) {
		__int32 node_index;
		__float32 load_value;
		tie(node_index, load_value) = applied_loads[k];
		loads[node_index] = load_value;
	}
	__int32 supA_index;
	__int32 supB_index;
	__int32 critical_index = 0;
	tie(supA_index, supB_index) = support_reactions;
	vector supA_pos = nodes[supA_index];
	vector supB_pos = nodes[supB_index];
	__float32 supA_val = 0;
	__float32 supB_val = 0;
	__float32 load_sum = 0;

	printf("Applied loads (%i):\n[", num_loads);
	for (__int32 k = 0; k < num_nodes; k++) {
		if (fabs(loads[k]) < TOL) {
			continue;
		}
		vector *load_pos = &nodes[k];
		__float32 delta_x = supB_pos.x - load_pos->x;
		supA_val += delta_x * loads[k];
		load_sum += loads[k];
		printf("(%.2f, %.2f) : %.2f, ", load_pos->x, load_pos->y, loads[k]);
	}
	printf("]\n\n");

	supA_val /= supA_pos.x - supB_pos.x;
	supB_val = -load_sum - supA_val;
	loads[supA_index] += supA_val;
	loads[supB_index] += supB_val;

	printf("Support reactions (2):\n(%.2f, %.2f) : %.2f\n(%.2f, %.2f) : %.2f\n\n",
	       supA_pos.x, supA_pos.y, supA_val, supB_pos.x, supB_pos.y, supB_val);
	return loads;
}

__int32 *find_critical_nodes(__int32 *truss_start, __int32 *truss_end, __int32 num_loads,
                             tuple<__int32, __int32> support_reactions,
                             tuple<__int32, __float32> *applied_loads) {
	__int32 supA_index;
	__int32 supB_index;
	__int32 *critical_nodes = new __int32[num_loads];
	tie(supA_index, supB_index) = support_reactions;
	*truss_start = supA_index;
	*truss_end = supB_index;
	for (__int32 k = 0; k < num_loads; k++) {
		__int32 load_index;
		__float32 load_value;
		tie(load_index, load_value) = applied_loads[k];
		critical_nodes[k] = load_index;
	}
	return critical_nodes;
}

vector unit_vector(vector *nodes, __int32 i, __int32 j) {
	vector r = nodes[j] - nodes[i];
	return r / r.norm();
}

__float32 node_length(vector *nodes, __int32 i, __int32 j) {
	vector r = nodes[j] - nodes[i];
	return r.norm();
}

void GSM_Solver::compute_solution(__int32 m, __int32 n, __float32 d, __int32 num_loads,
                                  tuple<__int32, __int32> support_reactions,
                                  tuple<__int32, __float32> *applied_loads) {
	printf("Hyperparameters:\nn = %i, m = %i, d = %.2f\n\n", m, n, d);

	__int32 num_nodes = m * n;
	__int32 num_struts = num_nodes * (num_nodes - 1) / 2;
	vector *nodes = generate_grid_nodes(m, n, d);
	bool **struts = generate_strut_matrix(m, n);
	__int32 **strut_indices = generate_strut_indices(m, n);
	__int32 **strut_inverse = generate_strut_inverse(m, n);
	__int32 truss_start = 0;
	__int32 truss_end = 0;
	__float32 *loads = get_loads(m, n, num_loads, nodes, support_reactions, applied_loads);
	__int32 *critical_nodes = find_critical_nodes(&truss_start, &truss_end, num_loads,
	                                              support_reactions, applied_loads);

	for (__int32 it = 0; it < 50; it++) {
		printf("<BEGIN ITERATION %i>\n", it);

		__int32 num_active_nodes = 0;
		__int32 num_active_struts = 0;
		bool visited_nodes[num_nodes]{0};
		for (__int32 i = 0; i < num_nodes; i++) {
			for (__int32 j = i; j < num_nodes; j++) {
				if (!struts[i][j]) {
					continue;
				}
				num_active_struts++;
				if (!visited_nodes[i]) {
					visited_nodes[i] = 1;
					num_active_nodes++;
				}
				if (!visited_nodes[j]) {
					visited_nodes[j] = 1;
					num_active_nodes++;
				}
			}
		}
		if (num_active_struts <= 0) {
			break;
		}

		printf("Active nodes: %i/%i\n", num_active_nodes, num_nodes);
		printf("Active struts: %i/%i\n\n", num_active_struts, num_struts);

		__int32 num_eq = 2 * num_nodes;
		__int32 index_eq = 0;
		sparse_matrix coeff_matrix(num_eq, num_struts);
		eigen_vector const_vector(num_eq);
		eigen_vector forces(num_struts);
		for (__int32 i = 0; i < num_eq; i++) {
			const_vector[i] = 0.0;
		}
		for (__int32 i = 0; i < num_nodes; i++) {
			if (fabs(loads[i]) >= TOL) {
				const_vector[2 * i + 1] = -loads[i];
			}
			for (__int32 j = 0; j < num_nodes; j++) {
				if (!struts[i][j]) {
					continue;
				}
				vector u_vec = unit_vector(nodes, i, j);
				if (fabs(u_vec.x) < TOL) {
					u_vec.x = 0.0;
				}
				if (fabs(u_vec.y) < TOL) {
					u_vec.y = 0.0;
				}
				coeff_matrix.insert(index_eq, strut_indices[i][j]) = u_vec.x;
				coeff_matrix.insert(index_eq + 1, strut_indices[i][j]) = u_vec.y;
			}
			index_eq += 2;
		}
		coeff_matrix.makeCompressed();

		linear_solver solver;
		solver.analyzePattern(coeff_matrix);
		solver.factorize(coeff_matrix);
		forces = solver.solve(const_vector);

		printf("Member forces are:\n");
		bool found_zero = false;
		__int32 zero_count = 0;
		for (__int32 k = 0; k < num_struts; k++) {
			__int32 *strut_index = strut_inverse[k];
			if (fabs(forces[k]) < TOL) {
				if (struts[strut_index[0]][strut_index[1]]) {
					zero_count++;
					found_zero = true;
					struts[strut_index[0]][strut_index[1]] = 0;
					struts[strut_index[1]][strut_index[0]] = 0;
				}
			} else {
				printf("(%i, %i) : %.5f, ", strut_index[0], strut_index[1], forces[k]);
			}
		}
		printf("\n\n");

		__float32 total_cost = num_active_nodes * COST_GUSSET_PLATE;
		__float32 min_efficiency = -1;
		__int32 min_index = -1;
		for (__int32 k = 0; k < num_struts; k++) {
			__float32 abs_force = fabs(forces[k]);
			if (abs_force < TOL) {
				continue;
			}
			__int32 *strut_index = strut_inverse[k];
			__float32 strut_length = node_length(nodes, strut_index[0], strut_index[1]);
			__float32 strut_efficiency = abs_force / strut_length;
			total_cost += strut_length * COST_MEMBER_LENGTH;
			if (min_efficiency == -1 || strut_efficiency < min_efficiency) {
				min_efficiency = strut_efficiency;
				min_index = k;
			}
			if (forces[k] < MIN_FORCE || forces[k] > MAX_FORCE) {
				printf("(%i, %i) : %.5f, ", strut_index[0], strut_index[1], forces[k]);
			}
		}
		printf("\n\nTruss cost: $%.2f\n\n", total_cost);
		if (found_zero) {
			printf("Found and removed %i zero-members\n\n", zero_count);
		} else {
			__int32 *remove_index = strut_inverse[min_index];
			struts[remove_index[0]][remove_index[1]] = 0;
			struts[remove_index[1]][remove_index[0]] = 0;
			printf("Removing strut %i : (%i, %i)\nWith efficiency %f\n\n",
			       min_index, remove_index[0], remove_index[1], min_efficiency);
		}

		int_queue next_nodes;
		bool found_nodes[num_nodes]{0};
		__int32 found_nodes_count = 1;
		next_nodes.push(truss_start);
		while(next_nodes.size()) {
			__int32 node_ptr = next_nodes.front();
			next_nodes.pop();
			found_nodes[node_ptr] = 1;
			for (__int32 k = 0; k < num_nodes; k++) {
				if (struts[node_ptr][k] && !found_nodes[k]) {
					found_nodes_count++;
					found_nodes[k] = 1;
					next_nodes.push(k);
				}
			}
		}
		if (found_nodes_count < num_active_nodes) {
			printf("Truss contains disjoint nodes (%i / %i)\n", found_nodes_count, num_active_nodes);
		}
		if (!found_nodes[truss_end]) {
			printf("Truss does not connect supported ends (%i, %i)\n", truss_start, truss_end);
		}
		for (__int32 k = 0; k < num_loads; k++) {
			if (!found_nodes[critical_nodes[k]]) {
				printf("Critical node %i unconnected\n", critical_nodes[k]);
			}
		}
		printf("\n\n");
	}

	delete[] nodes;
	delete[] loads;
	for (__int32 k = 0; k < num_nodes; k++) {
		delete[] struts[k];
		delete[] strut_indices[k];
	}
	for (__int32 k = 0; k < num_struts; k++) {
		delete[] strut_inverse[k];
	}
	delete[] struts;
	delete[] strut_indices;
	delete[] strut_inverse;
}