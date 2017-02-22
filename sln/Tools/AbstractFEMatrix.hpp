#pragma once
#include <set>
#include "AbstractSparseMatrix.hpp"
#include "AbstractFiniteElement.hpp"

/*
	Alexander Žilyakov, Sep 2016
*/

// in order to build patterns of sparse matrices,
// we should be able to describe how DOFs are connected
using DOFsConnectivityList = std::vector<std::set<Index>>;

template <LocalIndex D, LocalIndex N, LocalIndex M>
inline auto createDOFsConnectivityList(
	AbstractMesh<D, N> const & mesh,
	// FE to get rows numeration from
	AbstractFiniteElement<D, N, M> const & colFE,
	// FE to get rows numeration from
	AbstractFiniteElement<D, N, M> const & rowFE,
	// predicate defining when to connect col and row, e.g.
	// for matrices w/ symmetric pattern it is sufficient to add only those row indicies 
	// which are more than col index
	std::function<bool(Index, Index)> const & pred = [](Index, Index) { return true; }
) {
	DOFsConnectivityList list(colFE.numbOfDOFs(mesh));
	for (Index e = 0; e < mesh.numbOfElements(); ++e) {
		auto cols = colFE.getDOFsNumeration(mesh, e);
		auto rows = rowFE.getDOFsNumeration(mesh, e);
		for (Index col : cols)
			for (Index row : rows)
				if (pred(row, col)) list[col].insert(row);
	}
	return list;
}

// simpler version
template <LocalIndex D, LocalIndex N, LocalIndex M>
inline auto createDOFsConnectivityList(
	AbstractMesh<D, N> const & mesh,
	AbstractFiniteElement<D, N, M> const & FE,
	std::function<bool(Index, Index)> const & pred = [](Index row, Index col) { return row > col; }
) {
	return createDOFsConnectivityList(mesh, FE, FE, pred);
}

#include "Triangulation.hpp"
inline auto createDOFsConnectivityList(
	Triangulation const & cMesh,
	Triangulation const & fMesh,
	TriangularScalarFiniteElement const & FE
) {
	DOFsConnectivityList list(FE.numbOfDOFs(fMesh));
	for (Index ci = 0; ci < cMesh.numbOfElements(); ++ci) { // ci for “coarse index”
		auto rows = FE.getDOFsNumeration(cMesh, ci);
		for (Index j : FE.getDOFsNumeration(fMesh, ci))
			for (Index i : rows)
				list[j].insert(i);
		for (Index fi : fMesh.getNeighborsIndicies(ci)) // fi for “fine index”
			for (Index j : FE.getDOFsNumeration(fMesh, fi))
				for (Index i : rows)
					list[j].insert(i);
	}
	return list;
}

// for strong (Dirichlet) BCs
// key = DOF index where strong BC should be enforced, value = value of DOF 
template <typename T>
using Index2Value = std::unordered_map<Index, T>;

template <typename T>
class AbstractFEMatrix
	: virtual public AbstractSparseMatrix<T> {
public:
	// construct matrix pattern from adjacency list 
	// of finite element mesh verticies
	virtual AbstractFEMatrix& generatePatternFrom(DOFsConnectivityList const &) = 0;
	// enforce strong BCs
	virtual AbstractFEMatrix& enforceDirichletBCs(Index2Value<T> const &, T*) = 0;
};