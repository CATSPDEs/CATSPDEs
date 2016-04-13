#include"DenseMatrix.hpp"
#include"SymmetricContainer.hpp"
#include"Triangulation.hpp"
#include"BC.hpp"

Node dummy;

inline double a(Node const& p){
	return 1;
}

inline double f(Node const& p){
	return 1;
}

inline double g_D(Node const& p){
	return 0;
}

inline double g_N(Node const& p){
	return 0;
}

inline double kappa(Node const& p) {
	return 1;
}

//double kappa = 1.;

int main() {
	Triangulation Omega(Node(-1,-1),Node(1,1),.99);
	DenseMatrix K(Omega.numbOfNodes());
	SymmetricContainer<double> K_loc(3);
	SymmetricContainer<double> robinLocalMatrix(2);
	array<double,3> b,c;

	Node leftNode, rightNode;

	array<Node, 3> nodes;
	
	array<size_t, 3> loc2glob;
	vector<double> L(Omega.numbOfNodes(), 0);
	array<double, 3> l_loc;
	array<double, 2> robinLocalMatrix;
	LocalIndicies boundaryIndicies;
	for (size_t i = 0; i < Omega.numbOfTriangles(); i++) {

		nodes = Omega.getNodes(i);
		
		b = { (nodes[1] - nodes[2]).y() / 2 / Omega.area(i),(nodes[2] - nodes[0]).y() / 2 / Omega.area(i),(nodes[0] - nodes[1]).y() / 2 / Omega.area(i) };
		c = { (nodes[2] - nodes[1]).x() / 2 / Omega.area(i),(nodes[0] - nodes[2]).x() / 2 / Omega.area(i),(nodes[2] - nodes[1]).x() / 2 / Omega.area(i) };
		l_loc = { f(nodes[0]) / 3 * Omega.area(i),f(nodes[1]) / 3 * Omega.area(i),f(nodes[2]) / 3 * Omega.area(i) };
		for (localIndex j = 0;j < 3;j++)
			for (localIndex k = j;k < 3;k++)
				K_loc(j, k) = a(1 / 3 * (nodes[0] + nodes[1] + nodes[2]))*(b[j] * b[k] + c[j] * c[k])*Omega.area(i);
		boundaryIndicies = Omega.getBoundaryIndicies(i);
		for (auto edgeIndex : boundaryIndicies)
		{
			// left and right nodes of the edge
			leftNode = nodes[(edgeIndex + 1) % 3];
			rightNode = nodes[(edgeIndex + 2) % 3];
			// assemble Robin local matrix
			robinLocalMatrix(0, 0) = Omega.length(i, edgeIndex) / 12 * (kappa(leftNode) + 3 * kappa(rightNode));
			robinLocalMatrix(0, 1) = Omega.length(i, edgeIndex) / 12 * (kappa(leftNode) + kappa(rightNode));
			robinLocalMatrix(1, 1) = Omega.length(i, edgeIndex) / 12 * (3 * kappa(leftNode) + kappa(rightNode));
			// assemble Robin local vector
			robinLocalVector = { Omega.length(i, edgeIndex) / 12 };
		}
		for (localIndex j = 0;j < 3;j++){
			for (localIndex k = 0;k < 3;k++)
				K(Omega.loc2glob(i)[j], Omega.loc2glob(i)[k]) += K_loc(j, k)+robinLocalMatrix(j,k);
			L[Omega.loc2glob(i)[j]] += l_loc[j]+robinLocalMatrix[j];
		}
	}
	K.save();
	cout << L << endl;
	cout << L / K << endl;
}