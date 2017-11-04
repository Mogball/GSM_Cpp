#ifndef GSM_GSM_SOLVER_H
#define GSM_GSM_SOLVER_H

#define TOL 1e-6
#define COST_GUSSET_PLATE 6.0;
#define COST_MEMBER_LENGTH 12.0;
#define MIN_FORCE -8.0
#define MAX_FORCE 13.0

#include <tuple>
#include "vector.h"

using std::tuple;
using std::tie;
using std::make_tuple;

inline vector unit_vector(vector *nodes, __int32 i, __int32 j) {
	vector r = nodes[j] - nodes[i];
	return r / r.norm();
}

class GSM_Solver {
public:
	virtual void compute_solution(__int32 m, __int32 n, __float32 d, __int32 num_loads,
	                              tuple<__int32, __int32> support_reactions,
	                              tuple<__int32, __float32> *applied_loads);
};

#endif //GSM_GSM_SOLVER_H
