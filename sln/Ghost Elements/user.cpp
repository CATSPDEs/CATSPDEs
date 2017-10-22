#include "SingletonLogger.hpp"
#include "Triangulation.hpp"

using std::vector;
using std::string;

int main() {
	string oPath("Mathematica/postprocessing/");
	auto& logger = SingletonLogger::instance();
	try {
		logger.beg("create a simple triangulation of right triangles");
			Triangulation Omega {
				Node2D { 0., 0. }, // left bottom node
				Node2D { 8., 4.}, // right top node
				4, 4 // numb of horizontal and vertical intervals
			};
			Omega.export(oPath + "rect.nt");
		logger.end();
		logger.beg("truncate mesh");
			// our region is a union of region1, region2, region3, and region4 
			Triangle2D region1 {
				{ { 0., 3. }, { 2., 3. }, { 0., 4.} }
			};
			Quadrilateral2D region2 {
				{ { 0., 2. }, { 2., 2. }, { 2., 3. }, { 0., 3. } }
			};
			Quadrilateral2D region3 {
				{ { 0., 0. }, { 8., 0. }, { 8., 2. }, { 0., 2. } }
			};
			Quadrilateral2D region4 {
				{ { 6., 2. }, { 8., 2. }, { 8., 3. }, { 6., 3. } }
			};
			auto regionPredicate = [&](Node2D const & p) {
				return
					nodeInElement(region1, p) ||
					nodeInElement(region2, p) ||
					nodeInElement(region3, p) ||
					nodeInElement(region4, p);
			};
			// create ghost elements
			Omega.truncate(regionPredicate);
		logger.end();
		logger.beg("uniform refinement");
			Omega.refine(2);
			Omega.export(oPath + "uniform.nt");
		logger.end();
		logger.beg("non-uniform refinement");
			// create list of elements that are in the union of region1 and region2
			Indicies elems2refine;
			for (Index t = 0; t < Omega.numbOfElements(); ++t) {
				auto p = centroid(Omega.getElement(t));
				if (nodeInElement(region1, p) || nodeInElement(region2, p)) elems2refine.push_back(t);
			}
			// refine these elements
			Omega.refine(elems2refine);
			Omega.export(oPath + "nonuniform.nt");
		logger.end();
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}