#pragma once
#include "AbstractFiniteElement.hpp"

class Triangle_P1_Lagrange 
	: public TriangularScalarFiniteElement {
	// singleton
	Triangle_P1_Lagrange() : AbstractFiniteElement(1.) {}
	Triangle_P1_Lagrange(Triangle_P1_Lagrange const &);
	Triangle_P1_Lagrange& operator=(Triangle_P1_Lagrange const &);
public:
	static auto& instance() {
		static Triangle_P1_Lagrange single;
		return single;
	}
	std::vector<SmartScalarField2D> getShapesOf(Triangle const & t) const final {
		return {
			[&](Node2D const & p) { 
				return ((p[1] - t[1][1])*t[2][0] + p[0] * (t[1][1] - t[2][1]) + t[1][0] * (-p[1] + t[2][1])) / 2. / area(t); 
			},
			[&](Node2D const & p) {
				return ((-p[1] + t[0][1])*t[2][0] + t[0][0] * (p[1] - t[2][1]) + p[0] * (-t[0][1] + t[2][1])) / 2. / area(t);
			},
			[&](Node2D const & p) {
				return ((p[1] - t[0][1])*t[1][0] + p[0] * (t[0][1] - t[1][1]) + t[0][0] * (-p[1] + t[1][1])) / 2. / area(t);
			}
		};
	}
	std::vector<SmartVectorField2D> getSGradsOf(Triangle const & t) const final {
		return {
			[&](Node2D const &) { 
				return Node2D { { t[1][1] - t[2][1], -t[1][0] + t[2][0] } } / 2. / area(t);
			},
			[&](Node2D const & p) {
				return Node2D { { -t[0][1] + t[2][1], t[0][0] - t[2][0] } } / 2. / area(t);
			},
			[&](Node2D const & p) {
				return Node2D { { t[0][1] - t[1][1], -t[0][0] + t[1][0] } } / 2. / area(t);
			}
		};
	}
};