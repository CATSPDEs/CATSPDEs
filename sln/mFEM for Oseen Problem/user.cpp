#include "SingletonLogger.hpp"
#include "Node.hpp"

// logger
SingletonLogger& logger = SingletonLogger::instance();

int main() {
	try {
		Node2D a = {1., 2.};
		Node<1> b = 3.;
		logger.buf << midNode(a,a);
		logger.log();
	}
	catch (std::exception const & e) {
		logger.err(e.what());
	}
	return 0;
}