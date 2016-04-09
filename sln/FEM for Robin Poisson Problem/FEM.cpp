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

double kappa = 1.;

int main() {
	Triangulation Omega(Node(-1,-1),Node(1,1),.99);
	DenseMatrix K(Omega.numbOfNodes());
	SymmetricContainer<double> K_loc(3);
	SymmetricContainer<double> R_loc(3);
	array<double,3> b,c;
	array<Node, 3> nodes;
	array<size_t, 3> loc2glob;
	vector<double> L(Omega.numbOfNodes(), 0);
	array<double, 3> l_loc;
	array<double, 3> r_loc;
	LocalIndicies boundary;
	for (size_t i = 0;i < Omega.numbOfTriangles();i++) {
		nodes = Omega.getNodes(i);
		b = { (nodes[1] - nodes[2]).y() / 2 / Omega.area(i),(nodes[2] - nodes[0]).y() / 2 / Omega.area(i),(nodes[0] - nodes[1]).y() / 2 / Omega.area(i) };
		c = { (nodes[2] - nodes[1]).x() / 2 / Omega.area(i),(nodes[0] - nodes[2]).x() / 2 / Omega.area(i),(nodes[2] - nodes[1]).x() / 2 / Omega.area(i) };
		l_loc = { f(nodes[0]) / 3 * Omega.area(i),f(nodes[1]) / 3 * Omega.area(i),f(nodes[2]) / 3 * Omega.area(i) };
		for (localIndex j = 0;j < 3;j++)
			for (localIndex k = j;k < 3;k++)
				K_loc(j, k) = a(1 / 3 * (nodes[0] + nodes[1] + nodes[2]))*(b[j] * b[k] + c[j] * c[k])*Omega.area(i);
		boundary = Omega.getBoundaryIndicies(i);
		for (auto edge : boundary)
		{
			for (localIndex j = 0;j < 3;j++)
				for (localIndex k = j;k < 3;k++)
					R_loc(j, k) = kappa / 6 * Omega.length(i, edge);
			for (localIndex j = 0;j < 3;j++)
				R_loc(j, j) *= 2;
			r_loc = { 1 / 2 * (kappa*g_D(dummy) + g_N(dummy))*Omega.length(i,edge),1 / 2 * (kappa*g_D(dummy) + g_N(dummy))*Omega.length(i,edge) ,1 / 2 * (kappa*g_D(dummy) + g_N(dummy))*Omega.length(i,edge) };
		}
		for (localIndex j = 0;j < 3;j++){
			for (localIndex k = 0;k < 3;k++)
				K(Omega.loc2glob(i)[j], Omega.loc2glob(i)[k]) += K_loc(j, k)+R_loc(j,k);
			L[Omega.loc2glob(i)[j]] += l_loc[j]+r_loc[j];
		}
	}
	K.save();
	cout << L << endl;
	cout << L / K << endl;
}