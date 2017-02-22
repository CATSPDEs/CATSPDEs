#include "Triangulation.hpp"

bool Triangulation::_makeNeighbors(Index t1, Index t2) {
	// make triangles t1 and t2 neighbors
	// if they are not adjacent, return false
	if (t1 == t2) return false;
	auto u = getCommonNodesLocalIndicies(t1, t2);
	if (u.size() > 2) // triangles share… more than 2 nodes?
		throw std::logic_error("invalid mesh: check out tringles #" + std::to_string(t1) + " and #" + std::to_string(t2));
	else if (u.size() == 2) {
		(*_neighbors)[t1][excludeIndicies(u[0][0], u[1][0])] = t2;
		(*_neighbors)[t2][excludeIndicies(u[0][1], u[1][1])] = t1;
		return true; // we are gold
	}
	return false;
}

Triangulation& Triangulation::import(std::istream& from) {
	std::string meshType;
	from >> meshType;
	Index n, t; // numb of nodes and ″ triangles
	from >> n >> t;
	_nodes.resize(n);
	_elements.resize(t);
	if (meshType == "NT") // nodes and triangles
		from >> _nodes >> _elements;
	else if (meshType == "NTN") { // nodes, triangles, and neighbors
		from >> _nodes >> _elements;
		std::vector<std::array<SignedIndex, 3>> neighbors(t);
		from >> neighbors;
		_neighbors.emplace(neighbors);
	}
	else if (meshType == "NTR") { // nodes, triangles, and ribs
		Index r;
		decltype(_ribs->second) ribs(t);
		from >> r >> _nodes >> _elements >> ribs;
		_ribs.emplace(std::make_pair(r, ribs));
	}
	else if (meshType == "NTNR") { // nodes, triangles, neighbors, and ribs
		Index r;
		std::vector<std::array<SignedIndex, 3>> neighbors(t);
		decltype(_ribs->second) ribs(t);
		from >> r >> _nodes >> _elements >> neighbors >> ribs;
		_neighbors.emplace(neighbors);
		_ribs.emplace(std::make_pair(r, ribs));
	}
	else throw std::invalid_argument("unknown mesh type");
	return *this;
}

void Triangulation::export(std::ostream& to, Parameters const & params) const {
	std::string format = "NT";
	if (params.find("format") != params.end()) format = params.at("format");
	if (format == "NT") // nodes and triangles
		to << "NT\n" << numbOfNodes() << ' ' << numbOfElements() << '\n' << _nodes << _elements;
	else if (format == "NTN") // nodes, triangles, and neighbors
		to << "NTN\n" << numbOfNodes() << ' ' << numbOfElements() << '\n' << _nodes << _elements << *_neighbors;
	else if (format == "NTNR") // nodes, triangles, neighbors, and ribs
		to << "NTNR\n" << numbOfNodes() << ' ' << numbOfElements() << ' ' << numbOfRibs() << '\n' << _nodes << _elements << *_neighbors << (*_ribs).second;
	else if (format == "NTR") // nodes, triangles, and ribs
		to << "NTR\n" << numbOfNodes() << ' ' << numbOfElements() << ' ' << numbOfRibs() << '\n' << _nodes << _elements << (*_ribs).second;
	else
		throw std::invalid_argument("unknown mesh format: try NT, NTN, NTNR, or NTR");
}

Triangulation& Triangulation::refine(Index numbOfRefinements) {
	/*
	author:
	Alexander Žilyakov, Jun 2016
	edited:
	Oct 2016
	comments:
	here we implement uniform refinement
	works just like red green refinement, but
	ALL triangles will be refined in red fashion
	so we do not need stuff like redList, greenMap etc. here
	*/
	Indicies redNeighborsIndicies;
	Index i;
	auto getRedNodeIndexFrom = [&](Index t) {
		LocalIndex m;
		for (m = 0; m < 3; ++m)
			if (getNeighborsIndicies(getNeighborsIndicies(t)[m])[1] != i &&
				getNeighborsIndicies(getNeighborsIndicies(t)[m])[2] != i) break;
		redNeighborsIndicies.push_back(getNeighborsIndicies(t)[nextIndex(m)]);
		redNeighborsIndicies.push_back(getNeighborsIndicies(t)[nextIndex(nextIndex(m))]);
		return getNodesIndicies(t)[m];
	};
	for (Index k = 0; k < numbOfRefinements; ++k) {
		Index n = numbOfElements();
		for (i = 0; i < n; ++i) {
			auto nodesIndicies = getNodesIndicies(i);
			auto neighborsIndicies = getNeighborsIndicies(i);
			std::array<Index, 3> redNodesIndicies;
			for (LocalIndex j : {0, 1, 2})
				if (0 <= neighborsIndicies[j] && neighborsIndicies[j] < i) redNodesIndicies[j] = getRedNodeIndexFrom(neighborsIndicies[j]);
				else {
					redNodesIndicies[j] = numbOfNodes();
					if (neighborsIndicies[j] < -1) {
						//m = _neighbor2edge(t[j]);
						//_nodes.push_back(_curves[_curvilinearEdges[m].curveIndex()](_curvilinearEdges[m].thetaMiddle()));
						//_curvilinearEdges.push_back(CurvilinearEdge(_curvilinearEdges[m].thetaMiddle(), _curvilinearEdges[m].thetaEnd(), _curvilinearEdges[m].curveIndex()));
						//_curvilinearEdges[m].thetaEnd() = _curvilinearEdges[m].thetaMiddle();
					}
					else _nodes.push_back(midNodes(getElement(i))[j]);
				}
				_elements[i] = redNodesIndicies;
				SignedIndex dummy = numbOfElements();
				(*_neighbors)[i] = { dummy, dummy + 1, dummy + 2 };
				_elements.insert(_elements.end(), {
					{ nodesIndicies[0], redNodesIndicies[2], redNodesIndicies[1] },
					{ nodesIndicies[1], redNodesIndicies[0], redNodesIndicies[2] },
					{ nodesIndicies[2], redNodesIndicies[1], redNodesIndicies[0] }
				});
				dummy = i;
				(*_neighbors).insert((*_neighbors).end(), {
					{ dummy, neighborsIndicies[1], neighborsIndicies[2] },
					{ dummy, neighborsIndicies[2], neighborsIndicies[0] },
					{ dummy, neighborsIndicies[0], neighborsIndicies[1] }
				});
				//j = _curvilinearEdges.size();
				//for (Index newTriangle : { _triangles.size() - 2, _triangles.size() - 3, _triangles.size() - 1 })
				//	if (_triangles[newTriangle].neighbors(1) < -1)
				//		_triangles[newTriangle].neighbors(1) = _neighbor2edge(--j);
				for (Index j : { numbOfElements() - 3, numbOfElements() - 2, numbOfElements() - 1 })
					for (Index m : redNeighborsIndicies)
						_makeNeighbors(j, m);
				redNeighborsIndicies.clear();
		}
	}
	if (_ribs) enumerateRibs();
	return *this;
}

boost::optional<std::array<LocalIndex, 2>> Triangulation::getCommonRibLocalIndicies(Index t1, Index t2) const {
	/*
	author:
	Alexander Žilyakov, Oct 2016
	edited:
	_
	comments:
	if triangles #t1 and #t2 are neighbors
	we want to get local indicies of common rib:
	____________
	/2\0         2/
	/   \   t2    /
	/     \       /
	/  t1   \     /
	/         \   /
	/0_________1\1/

	for the stencil above, we get {0, 2}
	*/
	for (LocalIndex i = 0; i < 3; ++i)
		if (getNeighborsIndicies(t1)[i] == t2)
			for (LocalIndex j = 0; j < 3; ++j)
				if (getNeighborsIndicies(t2)[j] == t1)
					return std::array<LocalIndex, 2>{ i, j };
	return boost::none;
}

std::vector<std::array<LocalIndex, 2>> Triangulation::getCommonNodesLocalIndicies(Index t1, Index t2) const {
	/*
	author:
	Alexander Žilyakov, Oct 2016
	edited:
	_
	comments:
	we want to get local indicies of common nodes triangle #t1 and #t2 share:

	/2\
	/   \
	/     \
	/  t1   \
	/         \
	/0_________1\_____________
	\0         2/
	\         /
	\  t2   /
	\     /
	\   /
	\1/

	for the stencil above, we get { {1, 0} }
	*/
	std::vector<std::array<LocalIndex, 2>> res;
	for (LocalIndex i = 0; i < 3; ++i)
		for (LocalIndex j = 0; j < 3; ++j)
			if (getNodesIndicies(t1)[i] == getNodesIndicies(t2)[j]) res.push_back({ i, j });
	return res;
}

Triangulation& Triangulation::computeNeighbors() {
	_neighbors.emplace(std::vector<std::array<SignedIndex, 3>>(numbOfElements(), { -1, -1, -1 }));
	for (Index i = 0; i < numbOfElements(); ++i)
		for (Index j = i + 1; j < numbOfElements(); ++j)
			_makeNeighbors(i, j);
	return *this;
}

Triangulation& Triangulation::enumerateRibs() {
	std::vector<std::array<Index, 3>> numeration(numbOfElements());
	Index currentRibIndex = 0;
	for (Index i = 0; i < numbOfElements(); ++i)
		for (LocalIndex j : { 0, 1, 2 }) {
			auto n = getNeighborsIndicies(i)[j];
			// if neighbor of ith triangle has not been reached yet or does not exist, then we should numerate
			if (i < n || n < 0)
				numeration[i][j] = currentRibIndex++;
			// otherwise, rib already has an index so we need to use it
			else {
				auto k = getCommonRibLocalIndicies(n, i);
				if (!k) throw std::logic_error("invalid mesh: check out tringles #" + std::to_string(i) + " and #" + std::to_string(n));
				numeration[i][j] = numeration[n][(*k)[0]];
			}
		}
	std::pair<Index, std::vector<std::array<Index, 3>>> ribs { currentRibIndex, numeration };
	_ribs.emplace(ribs);
	return *this;
}