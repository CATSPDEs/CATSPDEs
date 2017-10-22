#include "SingletonLogger.hpp"
#include "SegmentElement.hpp"
#include "TriangleElement.hpp"
#include "QuadrilateralElement.hpp"

int main() {
	auto& logger = SingletonLogger::instance();
	try {
		logger.log("hello");

		SegmentElement1D I1;
		SegmentElement2D I2 { {0., 0.}, {1., 1.} };
		logger.buf << I1.measure() << " " << I2.measure();
		logger.log();

		TriangleElement2D T2;
		TriangleElement3D T3{ {0,0,2},{1,0,0},{0,1,0} };
		logger.buf << T3.measure() << " " << T3.diameter();
		logger.log();

		QuadrilateralElement2D Q2;
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}

	return 0;
}