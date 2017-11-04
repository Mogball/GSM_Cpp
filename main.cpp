#include <iostream>
#include <ctime>
#include <queue>
#include <Eigen>

#include "GSM_Solver.h"

typedef Eigen::SparseMatrix<double> sparse_matrix;
typedef Eigen::VectorXd eigen_vector;
typedef Eigen::LeastSquaresConjugateGradient<sparse_matrix> linear_solver;
typedef std::queue<__int32> int_queue;
typedef GSM_Solver truss_solver;

#define strut(i, j) \
    struts[i][j] = true; \
    struts[j][i] = true; \
    strut_inverse[strut_count][0] = i; \
    strut_inverse[strut_count][1] = j; \
    strut_indices[i][j] = strut_count; \
    strut_indices[j][i] = strut_count; \
    strut_count++;
#define pos(k, tx, ty) nodes[k].x = tx; nodes[k].y = ty; node_count++;

int main() {
	/*std::freopen(("logs_" + std::to_string(std::time(NULL)) + ".txt").c_str(), "w", stdout);

	__int32 m = 15;
	__int32 n = 5;
	__int32 num_loads = 5;
	__float32 d = 1.0;
	tuple<__int32, __int32> support_reactions = make_tuple(60, 74);
	tuple<__int32, __float32> applied_loads[]{
			make_tuple(31, 12.0),
			make_tuple(34, 12.0),
			make_tuple(37, 12.0),
			make_tuple(40, 12.0),
			make_tuple(43, 12.0)
	};

	truss_solver solver;
	solver.compute_solution(m, n, d, num_loads, support_reactions, applied_loads);*/

	__int32 num_nodes = 16;
	__int32 num_struts = 29;
	__float32 loads[num_nodes]{0};
	__float32 load = 56 / 9.0;
	__float32 sup = 56 / 2.0;
	loads[0] = load;
	loads[1] = load;
	loads[2] = load - sup;
	loads[3] = load;
	loads[4] = load;
	loads[5] = load;
	loads[6] = load - sup;
	loads[7] = load;
	loads[8] = load;
	bool struts[num_nodes][num_nodes]{{0}};
	__int32 strut_inverse[num_struts][2];
	__int32 strut_indices[num_nodes][num_nodes];
	__int32 strut_count = 0;
	strut(0, 1)
	strut(1, 2)
	strut(2, 3)
	strut(3, 4)
	strut(4, 5)
	strut(5, 6)
	strut(6, 7)
	strut(7, 8)
	strut(9, 0)
	strut(9, 1)
	strut(9, 2)
	strut(11, 2)
	strut(11, 3)
	strut(11, 4)
	strut(10, 9)
	strut(10, 2)
	strut(10, 11)
	strut(12, 4)
	strut(12, 5)
	strut(12, 6)
	strut(14, 6)
	strut(14, 7)
	strut(14, 8)
	strut(13, 12)
	strut(13, 6)
	strut(13, 14)
	strut(15, 11)
	strut(15, 4)
	strut(15, 12)
	assert(strut_count == num_struts);
	vector nodes[num_nodes];
	__int32 node_count = 0;
	for (__int32 b = 0; b <= 8; b++) {
		pos(b, b * 1.75, 0)
	}
	pos(9, 1.849609, 3.756348)
	pos(10, 3.616211, 5.171875)
	pos(11, 5.857422, 4.156250)
	pos(12, 10.383789, 6.016113)
	pos(13, 12.387695, 6.501709)
	pos(14, 13.750000, 3.746826)
	pos(15, 7.026367, 3.937256)
	assert(node_count == num_nodes);

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

	double total_cost = num_nodes * COST_GUSSET_PLATE;
	printf("Member forces are:\n");
	bool found_zero = false;
	__int32 zero_count = 0;
	vector node_sum[num_nodes];
	for (__int32 k = 0; k < num_nodes; k++) {
		node_sum[k].x = 0.0;
		node_sum[k].y = loads[k];
	}
	for (__int32 k = 0; k < num_struts; k++) {
		__int32 *strut_index = strut_inverse[k];
		__int32 i = strut_index[0];
		__int32 j = strut_index[1];
		vector r_vec = nodes[j] - nodes[i];
		double length = r_vec.norm();
		vector u_vec = r_vec / length;
		vector f_ij = u_vec * forces[k];
		vector f_ji = f_ij * -1;
		node_sum[i].x += f_ij.x;
		node_sum[i].y += f_ij.y;
		node_sum[j].x += f_ji.x;
		node_sum[j].y += f_ji.y;
		total_cost += length * COST_MEMBER_LENGTH;
		printf("(%i, %i) : %.5f\n", i, j, forces[k]);
	}
	printf("\n\n");
	for (__int32 k = 0; k < num_nodes; k++) {
		printf("%i : %.2f\n", k, node_sum[k]);
	}
	printf("\n\n");
	printf("Total cost: $%.2f\n", total_cost);
}