#pragma once
#include "SingletonLogger.hpp"
#include "Mapping.hpp"
#include "Element.hpp"

// D := dimension of mesh (nodes) 
// N := numb of nodes in an element

template <LocalIndex D, LocalIndex N>
class AbstractQuadratureRule {
protected:
	Element<D, N> _domainOfIntegration;
	std::vector<std::vector<Node<D>>> _nodes;
	std::vector<std::vector<double>>  _weights;
public:
	AbstractQuadratureRule(
		Element<D, N> const & d, 
		std::vector<std::vector<Node<D>>> const & n, 
		std::vector<std::vector<double>>  const & w
	) : _domainOfIntegration(d), _nodes(n), _weights(w)
	{}
	virtual ~AbstractQuadratureRule() {}
	// usually we invent quadrature rule for so-called master elements
	auto domainOfIntegration() const {
		return _domainOfIntegration;
	}
	// max degree of polynomials we can integrate exactly
	auto maxDeg() const {
		return _nodes.size() - 1;
	}
	auto getQuadratureNodes(LocalIndex deg) const {
		return _nodes[checkDeg(deg)];
	}
	double computeQuadrature(ScalarField<D> const & f, LocalIndex deg = 1) const {
		deg = checkDeg(deg);
		std::vector<double> images(_nodes[deg].size());
		std::transform(_nodes[deg].begin(), _nodes[deg].end(), images.begin(), [&](auto const & p) {
			return f(p);
		});
		return _weights[deg] * images;
	}
	LocalIndex checkDeg(LocalIndex deg) const {
		static auto wrn = true;
		if (wrn && deg > maxDeg()) {
			auto& logger = SingletonLogger::instance();
			logger.buf << "max polynomial degree for this quadrature rule is " << maxDeg() << " < " << deg << '\n'
			           << "max degree was changed to " << maxDeg();
			logger.wrn();
			wrn = false;
		}
		return min(deg, maxDeg());
	}
};
