#pragma once
#include "SingletonLogger.hpp"
#include "Mapping.hpp"

template <LocalIndex NODE_DIM>
class AbstractQuadratureRule {
protected:
	std::vector<std::vector<Node<NODE_DIM>>> const _nodes;
	std::vector<std::vector<double>> const _weights;
public:
	AbstractQuadratureRule(
		std::vector<std::vector<Node<NODE_DIM>>> const & nodes,
		std::vector<std::vector<double>>  const & weights
	) : _nodes(nodes), _weights(weights)
	{}
	virtual ~AbstractQuadratureRule() {}
	// max degree of polynomials we can integrate exactly
	auto maxDeg() const {
		return _nodes.size() - 1;
	}
	auto getQuadratureNodes(LocalIndex deg) const {
		return _nodes[checkDeg(deg)];
	}
	double computeQuadrature(ScalarField<NODE_DIM> const & f, LocalIndex deg = 1) const {
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
