#define EIGEN_NO_DEBUG

#include <Eigen>

#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>

typedef __int32 _us;
typedef __int32 _s;

#define COST_GUSSET_PLATE 6
#define COST_PER_LENGTH 12

const double tol = 1e-6;

struct vector {

	double x;
	double y;

	vector() : x(0.0), y(0.0) {}

	vector(double x, double y) : x(x), y(y) {}

	double norm() const {
		double r = x * x + y * y;
		return sqrt(r);
	}

	bool nonzero() const {
		return fabs(x) < tol && fabs(y) < tol;
	}

	vector &operator=(const vector &v) {
		x = v.x;
		y = v.y;
		return *this;
	}

	vector &operator+=(const vector &v) {
		this->x += v.x;
		this->y += v.y;
		return *this;
	}

	std::string to_string() const {
		std::stringstream ss;
		ss << '(' << x << ", " << y << ')';
		return ss.str();
	}
};

static inline vector operator-(const vector &a, const vector &b) {
	return vector(a.x - b.x, a.y - b.y);
}

static inline vector operator/(const vector &v, const double &s) {
	return vector(v.x / s, v.y / s);
}

static inline vector operator*(const vector &v, const double &s) {
	return vector(v.x * s, v.y * s);
}

static inline std::ostream &operator<<(std::ostream &os, const vector &v) {
	return os << v.to_string();
}

vector *generate_positions(_us m, _us n, double d) {
	vector *pos = new vector[m * n];
	_us i = 0;
	for (_us y = 0; y < n; y++)
		for (_us x = 0; x < m; x++, i++)
			pos[i] = vector(x * d, y * d);
	return pos;
}

_us node_index(_us i, _us j, _us m, _us n) {
	return i + j * m;
}

bool *generate_struts(_us m, _us n) {
	_us mn = m * n;
	_us n_struts = mn * (mn - (_us) 1) / (_us) 2;
	bool *struts = new bool[n_struts];
	for (_us i = 0; i < n_struts; i++) struts[i] = 1;
	return struts;
}

bool **generate_strut_matrix(_us m, _us n) {
	_us mn = m * n;
	bool **struts = new bool *[mn];
	for (_us a = 0; a < mn; a++) struts[a] = new bool[mn];
	// Fully connect the nodes except for the same nodes
	for (_us i = 0; i < mn; i++)
		for (_us j = 0; j < mn; j++)
			struts[i][j] = i != j;
	return struts;
}

vector unit_vector(vector *pos, _us i, _us j) {
	// Unit vector pointing from position i to position j
	vector r = pos[j] - pos[i];
	return r / r.norm();
}

double node_length(vector *pos, _us i, _us j) {
	vector r = pos[j] - pos[i];
	return r.norm();
}

int main() {
	std::freopen(("logs_" + std::to_string(std::time(NULL)) + ".txt").c_str(), "w", stdout);

	_us m = 16;
	_us n = 5;
	_us mn = m * n;
	_us n_struts = mn * (mn - (_us) 1) / (_us) 2;
	vector *pos = generate_positions(m, n, 1);
	bool **struts = generate_strut_matrix(m, n);
	_s **strut_indices = new _s *[mn];
	std::pair<_us, _us> *strut_indices_inverse = new std::pair<_us, _us>[n_struts];
	_us strut_index = 0;
	for (_us k = 0; k < mn; k++) strut_indices[k] = new _s[mn];
	for (_us i = 0; i < mn; i++)
		for (_us j = 0; j < mn; j++)
			if (i == j) strut_indices[i][j] = -1;
			else if (i > j) strut_indices[i][j] = strut_indices[j][i];
			else {
				strut_indices_inverse[strut_index].first = i;
				strut_indices_inverse[strut_index].second = j;
				strut_indices[i][j] = strut_index++;
			}

	_us n_vert_sup = 2;
	bool *vert_sup = new bool[mn];
	double *vert_sup_val = new double[mn];
	double *vert_load_val = new double[mn];
	for (_us i = 0; i < mn; i++) vert_sup[i] = 0;
	for (_us i = 0; i < mn; i++) vert_sup_val[i] = 0.0;
	for (_us i = 0; i < mn; i++) vert_load_val[i] = 0.0;
	vert_sup[32] = 1;
	vert_sup[47] = 1;
	vert_load_val[32] = 10;
	vert_load_val[35] = 10;
	vert_load_val[38] = 10;
	vert_load_val[41] = 10;
	vert_load_val[44] = 10;
	vert_load_val[47] = 10;

	// Find position of support A
	vector sup_A_pos;
	vector sup_B_pos;
	_us sup_A_index;
	_us sup_B_index;
	for (sup_A_index = 0; sup_A_index < mn; ++sup_A_index)
		if (vert_sup[sup_A_index]) {
			sup_A_pos = pos[sup_A_index];
			break;
		}
	// Find position of support B
	for (sup_B_index = sup_A_index + (_us) 1; sup_B_index < mn; ++sup_B_index)
		if (vert_sup[sup_B_index]) {
			sup_B_pos = pos[sup_B_index];
			break;
		}
	// Moment about point B and sum of y forces
	// (xA - xB) * fy = (xB - x0) * p0 + (xB - x1) * p1 + ... + (xB - xn) * pn
	double sup_A_val = 0;
	double sup_B_val = 0;
	double vert_load_sum = 0;
	std::cout << "Applied loads are:" << std::endl << "[";
	for (_us i = 0; i < mn; i++) {
		if (!vert_load_val[i]) continue;
		vector *load_pos = &pos[i];
		double delta_x = sup_B_pos.x - load_pos->x;
		sup_A_val += delta_x * vert_load_val[i];
		vert_load_sum += vert_load_val[i];
		std::cout << "(" << load_pos->x << ", " << load_pos->y << ") : " << vert_load_val[i] << ", ";
	}
	std::cout << "]" << std::endl << std::endl;
	sup_A_val /= sup_A_pos.x - sup_B_pos.x;
	sup_B_val = -vert_load_sum - sup_A_val;
	// Set the support values as loads
	vert_load_val[sup_A_index] += sup_A_val;
	vert_load_val[sup_B_index] += sup_B_val;

	std::cout << "Support reactions are" << std::endl;
	std::cout << sup_A_pos << " : " << sup_A_val << std::endl;
	std::cout << sup_B_pos << " : " << sup_B_val << std::endl << std::endl;

	// Build the coefficient matrix
	// Matrix contains 2 * mn equations and a total of mn * (mn - 1) / 2 unknowns
	// The linear system is sparse: each equation has max (mn - 1) unknowns
	for (int iteration = 0; iteration < n_struts; iteration++) {
		std::cout << "***BEGIN ITERATION " << iteration << "***" << std::endl;
		// Count the number of active nodes and active struts
		_us n_active_struts = 0;
		_us n_active_nodes = 0;
		bool marked_nodes[mn]{0};
		for (_us i = 0; i < mn; i++) {
			for (_us j = i; j < mn; j++) {
				if (!struts[i][j]) continue;
				n_active_struts++;
				if (!marked_nodes[i]) {
					marked_nodes[i] = 1;
					n_active_nodes++;
				}
				if (!marked_nodes[j]) {
					marked_nodes[j] = 1;
					n_active_nodes++;
				}
			}
		}
		if (n_active_struts <= 0) break;

		std::cout << "Active nodes: " << n_active_nodes << "/" << mn << std::endl;
		std::cout << "Active struts: " << n_active_struts << "/" << n_struts << std::endl;

		Eigen::MatrixXd coeff_matrix(2 * mn, n_struts);
		Eigen::VectorXd const_vector(2 * mn), forces(n_struts);
		for (_us i = 0; i < 2 * mn; i++) {
			const_vector[i] = 0.0;
			for (_us j = 0; j < n_struts; j++)
				coeff_matrix(i, j) = 0.0;
		}
		//std::cout << std::endl << "Unit vectors" << std::endl;
		_us eq_index = 0;
		for (_us i = 0; i < mn; i++) {
			// Set the load value in the constant vector
			if (fabs(vert_load_val[i]) >= tol)
				const_vector[2 * i + 1] = -vert_load_val[i];
			for (_us j = 0; j < mn; j++) {
				if (!struts[i][j]) continue;
				// Find unit vector between node i and node j
				vector u = unit_vector(pos, i, j);
				// Apply tolerance
				if (fabs(u.x) < tol) u.x = 0.0;
				if (fabs(u.y) < tol) u.y = 0.0;
				// Add the coefficents
				//std::cout << u << std::endl;
				coeff_matrix(eq_index, strut_indices[i][j]) = u.x;
				coeff_matrix(eq_index + 1, strut_indices[i][j]) = u.y;
			}
			eq_index += 2;
		}

		// Move values into the sparse matrix
		Eigen::SparseMatrix<double> sparse_coeff_matrix(2 * mn, n_struts);
		std::vector<double> coeff_vector;
		for (_us i = 0; i < 2 * mn; i++)
			for (_us j = 0; j < n_struts; j++)
				if (fabs(coeff_matrix(i, j)) > tol)
					sparse_coeff_matrix.insert(i, j) = coeff_matrix(i, j);
		sparse_coeff_matrix.makeCompressed();

		// Solve for the forces in the members
		//forces = coeff_matrix.fullPivHouseholderQr().solve(const_vector);
		Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
		solver.analyzePattern(sparse_coeff_matrix);
		solver.factorize(sparse_coeff_matrix);
		forces = solver.solve(const_vector);

		//std::cout << coeff_matrix << std::endl;

		// Eliminate zero_member beams
		bool hasZero = false; // indicate whether a zero member was found
		_us nonzeroCount = 0;
		std::vector<_us> removed;
		std::cout << "Member forces are:" << std::endl;
		for (_us k = 0; k < n_struts; k++) {
			std::pair<_us, _us> strut_i = strut_indices_inverse[k];
			if (fabs(forces[k]) < tol) {
				if (struts[strut_i.first][strut_i.second]) {
					hasZero = true;
					struts[strut_i.first][strut_i.second] = 0;
					struts[strut_i.second][strut_i.first] = 0;
					removed.push_back(k);
				}
			} else {
				std::cout << '(' << strut_i.first << ", " << strut_i.second << ')';
				std::cout << " : " << forces[k] << ", ";
			}
		}
		std::cout << std::endl << std::endl;
		//for (_us k = 0; k < removed.size(); k++)
		//std::cout << removed[k] << ", ";
		//std::cout << std::endl;
		//for (_us k = 0; k < mn; k++) {
		//    if (fabs(force_sum[k].x) < tol) force_sum[k].x = 0.0;
		//    if (fabs(force_sum[k].y) < tol) force_sum[k].y = 0.0;
		//    std::cout << k << " : " << force_sum[k] << std::endl;
		//}
		std::cout << std::endl;
		// If we did not find any new nonzero struts, remove the least efficent
		// Strut efficiency = abs(force in strut) / strut length
		double cost = n_active_nodes * COST_GUSSET_PLATE;
		if (!hasZero) {
			double min_efficiency = -1;
			_us min_efficiency_index = -1;
			for (_us k = 0; k < n_struts; k++) {
				if (fabs(forces[k]) < tol) continue;
				std::pair<_us, _us> strut_i = strut_indices_inverse[k];
				double member_length = node_length(pos, strut_i.first, strut_i.second);
				cost += member_length * COST_PER_LENGTH;
				double efficiency = fabs(forces[k]) / member_length;
				if (min_efficiency == -1 || efficiency < min_efficiency) {
					min_efficiency = efficiency;
					min_efficiency_index = k;
				}
			}
			std::pair<_us, _us> remove_strut_i = strut_indices_inverse[min_efficiency_index];
			struts[remove_strut_i.first][remove_strut_i.second] = 0;
			struts[remove_strut_i.second][remove_strut_i.first] = 0;

			std::cout << "Truss cost: " << cost << std::endl << std::endl;

			std::cout << "Removing strut " << min_efficiency_index << " : (" << remove_strut_i.first;
			std::cout << ", " << remove_strut_i.second << ")" << std::endl;
			std::cout << "With efficiency " << min_efficiency << " F/l" << std::endl << std::endl;
		}
	}

	delete[] pos;
	delete[] vert_sup;
	delete[] vert_sup_val;
	delete[] vert_load_val;
}