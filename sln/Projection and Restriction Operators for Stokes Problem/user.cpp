#include <map>
#include "SingletonLogger.hpp"
#include "Triangulation.hpp"
#include "CSCMatrix.hpp"

/*
	Alexander Žilyakov, Oct 2016
*/

//
//// Crouzeix–Raviart shapes
//
//double s1(Node const & p, array<Node, 3> const & t) { 
//	return (
//		t[2].x() * ( -2. * p.y() + t[0].y() + t[1].y() ) +
//		t[1].x() * (  2. * p.y() - t[0].y() - t[2].y() ) - 
//		( 2. * p.x() - t[0].x() ) * (t[1].y() - t[2].y() )
//	) / (
//		t[2].x() * (  t[0].y() - t[1].y() ) +
//		t[0].x() * (  t[1].y() - t[2].y() ) +
//		t[1].x() * ( -t[0].y() + t[2].y() )
//	);
//}
//
//double s2(Node const & p, array<Node, 3> const & t) {
//	return (
//		t[2].x() * ( 2. * p.y() - t[0].y() - t[1].y() ) +
//		( 2. * p.x() - t[1].x() ) * ( t[0].y() - t[2].y() ) +
//		t[0].x() * ( -2. * p.y() + t[1].y() + t[2].y() )
//	) / (
//		t[2].x() * (t[0].y() - t[1].y()) +
//		t[0].x() * (t[1].y() - t[2].y()) +
//		t[1].x() * (-t[0].y() + t[2].y())
//	);
//}
//
//double s3(Node const & p, array<Node, 3> const & t) {
//	return 1. + 2. * (
//		t[1].x() * ( -p.y() + t[0].y() ) +
//		t[0].x() * (  p.y() - t[1].y() ) +
//		p.x() * ( -t[0].y() + t[1].y() )
//	) / (
//		t[2].x() * (t[0].y() - t[1].y()) +
//		t[0].x() * (t[1].y() - t[2].y()) +
//		t[1].x() * (-t[0].y() + t[2].y())
//	);
//}

// Lagrange shapes

double l1(Node2D const & p, Triangle const & t) {
	return (
		t[2][0] * ( p[1]    - t[1][1] ) + 
		p[0]    * ( t[1][1] - t[2][1] ) +
		t[1][0] * ( t[2][1] - p[1]    )
	) / (
		t[0][0] * ( t[1][1] - t[2][1] ) +
		t[1][0] * ( t[2][1] - t[0][1] ) +
		t[2][0] * ( t[0][1] - t[1][1] )
	);
}

double l2(Node2D const & p, Triangle const & t) {
	return (
		t[2][0] * ( t[0][1] - p[1]    ) +
		p[0]    * ( t[2][1] - t[0][1] ) +
		t[0][0] * ( p[1]    - t[2][1] )
	) / (
		t[0][0] * ( t[1][1] - t[2][1] ) +
		t[1][0] * ( t[2][1] - t[0][1] ) +
		t[2][0] * ( t[0][1] - t[1][1] )
	);
}

double l3(Node2D const & p, Triangle const & t) {
	return (
		t[1][0] * ( p[1]    - t[0][1] ) +
		p[0]    * ( t[0][1] - t[1][1] ) +
		t[0][0] * ( t[1][1] - p[1]    )
	) / (
		t[0][0] * ( t[1][1] - t[2][1] ) +
		t[1][0] * ( t[2][1] - t[0][1] ) +
		t[2][0] * ( t[0][1] - t[1][1] )
	);
}

using std::array;
using std::string;
using std::to_string;
using std::vector;
using std::map;
using std::ifstream;
using std::ofstream;

int main() {
	array<double(*)(Node2D const &, Triangle const &), 3> shapes = { l1, l2, l3 };
	// path
	string iPath("Mathematica/Generate Mesh/"),
	       oPath("Mathematica/Projection and Restriction Visualization/");
	// logger
	SingletonLogger& logger = SingletonLogger::instance();
	try {
		logger.beg("load initial mesh from\n" + iPath);
			Triangulation Omega;
			Omega.import(iPath + "mesh.ntn");
			Omega.enumerateRibs();
			auto numbOfNodesCoarse = Omega.numbOfNodes(),
			     numbOfTrianglesCoarse = Omega.numbOfElements(),
			     numbOfRibsCoarse = Omega.numbOfRibs();
			auto ribsNumnCoarse = Omega.getRibsNumeration();
			// … compute system matrix …
			Omega.export(oPath + "meshes/0_mesh.ntr", { {"format", "NTR"} });
			//vector<double> xi;
			//for (Node const & p : Omega.computeMiddleNodes(ribsNumnCoarse)) xi.push_back(u(p));
			//ofstream oXi(oPath + "0_xi.dat");
			//oXi << xi;
		logger.end();
		Index numbOfMeshLevels;
		logger.inp("enter numb of mesh levels", numbOfMeshLevels);
		ofstream oNumbOfMeshLevels(oPath + "numbOfMeshLevels.dat");
		oNumbOfMeshLevels << numbOfMeshLevels;
		vector<CSCMatrix<double>> Restrictions; // restriction matrices
		for (Index currentMeshLevel = 1; currentMeshLevel <= numbOfMeshLevels; ++currentMeshLevel) {
			logger.beg("mesh level #" + to_string(currentMeshLevel));
				logger.beg("refine mesh");
					Omega.refine();
				logger.end();
				logger.beg("save fine mesh to\n" + oPath + "meshes");
					Omega.export(oPath + "meshes/" + to_string(currentMeshLevel) + "_mesh.ntr", {{ "format", "NTR" }});
				logger.end();
				// … compute system matrix …
				logger.beg("build pattern and fill in values of Restriction matrix");
					vector<Index> rowptr(3 * numbOfTrianglesCoarse + 1),
								  colind;
					vector<double> values(6 * numbOfTrianglesCoarse);
					// (1) for DOFs located inside coarse triangles
					fill(values.begin(), values.end(), .5);
					for (Index i = 0; i < rowptr.size(); ++i) rowptr[i] = 2 * i;
					for (Index t = 0; t < numbOfTrianglesCoarse; ++t) {
						colind.insert(colind.end(), { 
							ribsNumnCoarse[t][1], ribsNumnCoarse[t][2], 
							ribsNumnCoarse[t][0], ribsNumnCoarse[t][2],
							ribsNumnCoarse[t][0], ribsNumnCoarse[t][1]
						});
					}
					// (2) for cross-element DOFs
					for (Index t1 = 0; t1 < numbOfTrianglesCoarse; ++t1)
						for (SignedIndex n : Omega.getNeighborsIndicies(t1)) 
							// for (LocalIndex k : excludeIndex(Omega.getCommonRibLocalIndicies(n, t1).value()[0]))
							// for (LocalIndex k : { nextIndex(Omega.getCommonRibLocalIndicies(n, t1).value()[0]), nextIndex(nextIndex(Omega.getCommonRibLocalIndicies(n, t1).value()[0])) })
							for (LocalIndex k : { 1, 2 })
								if (Omega.getNeighborsIndicies(n)[k] < 0) // boundary
									rowptr.push_back(rowptr.back()); // zero row
								else {
									Index t2;
									for (SignedIndex i : Omega.getNeighborsIndicies(Omega.getNeighborsIndicies(n)[k]))
										if (0 <= i && i < numbOfTrianglesCoarse) {
											t2 = i;
											break;
										}
									if (t1 < t2) {
										rowptr.push_back(rowptr.back() + 5);
										map<Index, double> colind2value;
										for (LocalIndex i : {0, 1, 2}) {
											colind2value[ribsNumnCoarse[t2][i]] += .5 * shapes[i](midNodes(Omega.getElement(n))[k], Omega.getElement(t2));
											colind2value[ribsNumnCoarse[t1][i]] += .5 * shapes[i](midNodes(Omega.getElement(n))[k], Omega.getElement(t1));
										}
										for (auto const & kvp : colind2value) {
											colind.push_back(kvp.first);
											values.push_back(kvp.second);
										}
									}
								}
					Restrictions.emplace_back(rowptr, colind, values, numbOfRibsCoarse);
				logger.end();
				logger.beg("save Restriction matrix to\n" + oPath + "operators");
					Restrictions.back().saveHarwellBoeing(oPath + "operators/" + to_string(currentMeshLevel) + "_R.rra");
				logger.end();
				// save before the next refinement
				numbOfNodesCoarse = Omega.numbOfNodes();
				numbOfTrianglesCoarse = Omega.numbOfElements();
				numbOfRibsCoarse = Omega.numbOfRibs();
				ribsNumnCoarse = Omega.getRibsNumeration();
			logger.end();
		}
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}