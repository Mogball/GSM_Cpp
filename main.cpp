#include <iostream>
#include <ctime>

#include "GSM_Solver.h"

typedef GSM_Solver truss_solver;

// m  = 16
// n = 5
// 32, 47
// 32, 35, 38, 41, 44, 47 (10 kN each)
int main() {
	std::freopen(("logs_" + std::to_string(std::time(NULL)) + ".txt").c_str(), "w", stdout);

	__int32 m = 16;
	__int32 n = 5;
	__int32 num_loads = 6;
	__float32 d = 1;
	tuple<__int32, __int32> support_reactions = make_tuple(32, 47);
	tuple<__int32, __float32> applied_loads[]{
			make_tuple(32, 10.0),
			make_tuple(35, 10.0),
			make_tuple(38, 10.0),
			make_tuple(41, 10.0),
			make_tuple(44, 10.0),
			make_tuple(47, 10.0)
	};

	truss_solver solver;
	solver.compute_solution(m, n, d, num_loads, support_reactions, applied_loads);
}