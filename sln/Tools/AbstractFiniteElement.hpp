#pragma once
// to get DOFs numn for given mesh
#include "AbstractMesh.hpp"
// in order to store computed images of quadrature nodes
#include "SmartMapping.hpp"

// prolongation
template <typename T> class CSCMatrix;

// D := dimension of mesh (nodes) 
// N := numb of nodes in an element
// M := dimension of shapes

template <LocalIndex D, LocalIndex N, LocalIndex M>
class AbstractFiniteElement {
protected:
	double _deg;
	bool   _isConform;
public:
	AbstractFiniteElement(double deg, bool isConform) : _deg(deg), _isConform(isConform) {}
	virtual ~AbstractFiniteElement() {}
	// U_H := FE-space on the coarse mesh, U_h := ″ fine; return true if U_H in U_h
	bool isConform() const { return _isConform; }
	// get degree of polynomial space of shapes
	double deg() const { return _deg; };
	// get shape funcs of an element
	virtual std::vector<SmartMapping<D, M>>     getShapesOf(Element<D, N> const &) const = 0;
	// get gradients of ″
	virtual std::vector<SmartMapping<D, D * M>> getSGradsOf(Element<D, N> const &) const = 0;
	// DOFs
	virtual Index numbOfDOFs(AbstractMesh<D, N> const & mesh) const = 0;
	virtual std::vector<Index>   getDOFsNumeration(AbstractMesh<D, N> const & mesh, Index e) const = 0;
	virtual std::vector<Node<D>> getDOFsNodes     (AbstractMesh<D, N> const & mesh, Index e) const = 0;
	virtual std::vector<LocalIndex> getBndryDOFsLocalIndicies(AbstractMesh<D, N> const & mesh, LocalIndex b) const = 0;
	// default implementation of constructing prolongation matrix
	// from coarse mesh to fine mesh
	// suitable e.g. for Lagrange elements
	virtual void prolongate(CSCMatrix<double>& P, AbstractMesh<D, N> const & cMesh, AbstractMesh<D, N> const & fMesh) const;
	CSCMatrix<double> prolongate(AbstractMesh<D, N> const & cMesh, AbstractMesh<D, N> const & fMesh) const;
	virtual std::vector<double> getNodalValues(Mapping<D, M> const & f, AbstractMesh<D, N> const & mesh, Index e) const {
		// if (M != 1) throw std::logic_error("getNodalValues() should be overriden");
		auto nodes = getDOFsNodes(mesh, e);
		std::vector<double> values(nodes.size());
		std::transform(nodes.begin(), nodes.end(), values.begin(), [&](Node<D> const & p) { return f(p); });
		return values;
	}
};

#include "CSCMatrix.hpp"

template <LocalIndex D, LocalIndex N, LocalIndex M>
CSCMatrix<double> AbstractFiniteElement<D, N, M>::prolongate(AbstractMesh<D, N> const & cMesh, AbstractMesh<D, N> const & fMesh) const {
	CSCMatrix<double> P(numbOfDOFs(cMesh), numbOfDOFs(fMesh));
	P.generatePatternFrom(createDOFsConnectivityList(cMesh, fMesh, *this));
	prolongate(P, cMesh, fMesh);
	return P;
}

template <LocalIndex D, LocalIndex N, LocalIndex M>
void AbstractFiniteElement<D, N, M>::prolongate(CSCMatrix<double>& P, AbstractMesh<D, N> const & cMesh, AbstractMesh<D, N> const & fMesh) const {
	auto chop = [](double x) { return fabs(x) < 1e-10 ? 0. : x; };
	for (Index ci = 0; ci < cMesh.numbOfElements(); ++ci) {
		auto shapes = getShapesOf(cMesh.getElement(ci));
		auto nodes = getDOFsNodes(fMesh, ci);
		auto rows = getDOFsNumeration(cMesh, ci),
		     cols = getDOFsNumeration(fMesh, ci);
		for (LocalIndex j = 0; j < cols.size(); ++j)
			for (LocalIndex i = 0; i < rows.size(); ++i)
				P(rows[i], cols[j]) = chop(shapes[i](nodes[j]));
		for (Index fi : fMesh.getFineNeighborsIndicies(ci)) {
			cols = getDOFsNumeration(fMesh, fi);
			nodes = getDOFsNodes(fMesh, fi);
			for (LocalIndex j = 0; j < cols.size(); ++j)
				for (LocalIndex i = 0; i < rows.size(); ++i)
					P(rows[i], cols[j]) = chop(shapes[i](nodes[j]));
		}
	}
}

	using SegmentFiniteElement = AbstractFiniteElement<1, 2, 1>;

	using TriangularScalarFiniteElement = AbstractFiniteElement<2, 3, 1>;
	using TriangularVectorFiniteElement = AbstractFiniteElement<2, 3, 2>;


