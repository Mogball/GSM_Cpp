#ifndef GSM_VECTOR_H
#define GSM_VECTOR_H

#define __float32 double

#include <cmath>
#include <iostream>

class vector {
public:
	__float32 x;
	__float32 y;

	inline vector() : x(0.0), y(0.0) {}

	inline explicit vector(__float32 x, __float32 y) : x(x), y(y) {}

	inline __float32 norm() const {
		return sqrt(x * x + y * y);
	}
};

inline vector operator-(const vector &a, const vector &b) {
	return vector(a.x - b.x, a.y - b.y);
}

inline vector operator/(const vector &a, const __float32 &s) {
	return vector(a.x / s, a.y / s);
}

inline vector operator*(const vector &a, const __float32 &s) {
	return vector(a.x * s, a.y * s);
}

inline std::ostream &operator<<(std::ostream &os, const vector &v) {
	os << '(' << v.x << ',' << ' ' << v.y << ')';
}

#endif //GSM_VECTOR_H
