#include "Triangulation.hpp"

bool Triangulation::_makeNeighbors(Index t1, Index t2) {
	// make triangles t1 and t2 neighbors
	// if they are not adjacent, return false
	if (t1 == t2) return false;
	auto u = getCommonNodesLocalIndicies(t1, t2);
	if (u.size() > 2) // triangles share… more than 2 nodes?
		throw std::logic_error("invalid mesh: check out tringles #" + std::to_string(t1) + " and #" + std::to_string(t2));
	else if (u.size() == 2) {
		_neighbors[t1][excludeIndicies(u[0][0], u[1][0])] = t2;
		_neighbors[t2][excludeIndicies(u[0][1], u[1][1])] = t1;
		return true; // we are gold
	}
	return false;
}

SignedIndex Triangulation::_neighbor2edge(SignedIndex i) const {
	// convert element index to
	// index in _edges vector
	return - i - 2;
}

Triangulation::Triangulation(Node2D const & lb, Node2D const & rt, Index ndx, Index ndy) {
	/*
		author: 
			Alexander Zhiliakov, May 2016
		edited:
			Oct 2017
		comments:
			here we construct dummy rect triangulation
			w/ left bottom point @lb and right top point @rt
			elements will be right-angled triangles
			@ndx = numb of INTERVALS on x-axis and
			@ndy = ″ on y-axis,
			it will look something like this:
					 ________(3,2)
					| /| /| /|
					|/_|/_|/_|
					| /| /| /|
					|/_|/_|/_|
				(0,0)  
	*/
	Node2D diag = rt - lb;
	// so x-projection of size is width of our rect and y-projection is height
	if (diag[0] <= 0 || diag[1] <= 0) throw std::invalid_argument("invalid rect");
	if (ndx * ndy == 0) throw std::invalid_argument("invalid numb of elements");
	double dx = diag[0] / ndx, // legs sizes 
	       dy = diag[1] / ndy;
	Index nx = ndx + 1, // numb of NODES on x-axis and
	      ny = ndy + 1; // on y-axis, 
	_nodes.resize(nx * ny); // so numb of points is nx * ny
	for (Index i = 0; i < ny; ++i)
		for (Index j = 0; j < nx; ++j) 
			_nodes[i * nx + j] = lb + Node2D { j * dx, i * dy }; // compute nodes
	_elements.resize(2 * ndx * ndy);

	for (Index i = 0; i < ndy; ++i)
		for (Index j = 0; j < ndx; ++j) {
			Index k = i * ndx + j; // index of current rectangle
			Index n = i * nx + j; // index of left bottom node of this rectangle
			_elements[2 * k] = { n, n + 1, n + nx };
			_elements[2 * k + 1] = { n + 1, n + 1 + nx, n + nx };
		}
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
		_neighbors.resize(t);
		from >> _neighbors;
	}
	else if (meshType == "NTR") { // nodes, triangles, and ribs
		_ribs.second.resize(t);
		from >> _ribs.first >> _nodes >> _elements >> _ribs.second;
	}
	else if (meshType == "NTNR") { // nodes, triangles, neighbors, and ribs
		_neighbors.resize(t);
		_ribs.second.resize(t);
		from >> _ribs.first >> _nodes >> _elements >> _neighbors >> _ribs.second;
	}
	else throw std::invalid_argument("unknown mesh type");
	_fineNeighborsIndicies.resize(numbOfElements());
	return *this;
}

void Triangulation::export(std::ostream& to, Parameters const & params) const {
	std::string format = "NT";
	if (params.find("format") != params.end()) format = params.at("format");
	if (format == "NT") // nodes and triangles
		to << "NT\n" << numbOfNodes() << ' ' << numbOfElements() << '\n' << _nodes << _elements;
	else if (format == "NTN") // nodes, triangles, and neighbors
		to << "NTN\n" << numbOfNodes() << ' ' << numbOfElements() << '\n' << _nodes << _elements << _neighbors;
	else if (format == "NTNR") // nodes, triangles, neighbors, and ribs
		to << "NTNR\n" << numbOfNodes() << ' ' << numbOfElements() << ' ' << numbOfRibs() << '\n' << _nodes << _elements << _neighbors << _ribs.second;
	else if (format == "NTR") // nodes, triangles, and ribs
		to << "NTR\n" << numbOfNodes() << ' ' << numbOfElements() << ' ' << numbOfRibs() << '\n' << _nodes << _elements << _ribs.second;
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
		_fineNeighborsIndicies.resize(n);
		for (i = 0; i < n; ++i) {
			if (_ghostElements.find(i) != _ghostElements.end()) continue;
			auto nodesIndicies = getNodesIndicies(i);
			auto neighborsIndicies = getNeighborsIndicies(i);
			std::array<Index, 3> redNodesIndicies;
			for (LocalIndex j : {0, 1, 2})
				if (0 <= neighborsIndicies[j] && neighborsIndicies[j] < i) redNodesIndicies[j] = getRedNodeIndexFrom(neighborsIndicies[j]);
				else {
					redNodesIndicies[j] = numbOfNodes();
					_nodes.push_back(getRibNode(i, j, .5));
					SignedIndex edgeIndex = _neighbor2edge(neighborsIndicies[j]);
					if (edgeIndex > -1) {
						_edges.push_back({ _edges[edgeIndex].thetaMiddle(), _edges[edgeIndex].thetaEnd(), _edges[edgeIndex].curveIndex() });
						_edges[edgeIndex].thetaEnd() = _edges[edgeIndex].thetaMiddle();
					}
				}
			_elements[i] = redNodesIndicies;
			_fineNeighborsIndicies[i] = { numbOfElements(), numbOfElements() + 1, numbOfElements() + 2 };
			SignedIndex dummy = numbOfElements();
			_neighbors[i] = { dummy, dummy + 1, dummy + 2 };
			_elements.insert(_elements.end(), {
				{ nodesIndicies[0], redNodesIndicies[2], redNodesIndicies[1] },
				{ nodesIndicies[1], redNodesIndicies[0], redNodesIndicies[2] },
				{ nodesIndicies[2], redNodesIndicies[1], redNodesIndicies[0] }
			});
			dummy = i;
			_neighbors.insert(_neighbors.end(), {
				{ dummy, neighborsIndicies[1], neighborsIndicies[2] },
				{ dummy, neighborsIndicies[2], neighborsIndicies[0] },
				{ dummy, neighborsIndicies[0], neighborsIndicies[1] }
			});
			Index j = _edges.size();
			for (Index newTriangle : { _elements.size() - 2, _elements.size() - 3, _elements.size() - 1 })
				if (_neighbors[newTriangle][1] < -1)
					_neighbors[newTriangle][1] = _neighbor2edge(--j);

			for (Index j : { numbOfElements() - 3, numbOfElements() - 2, numbOfElements() - 1 })
				for (Index m : redNeighborsIndicies)
					_makeNeighbors(j, m);
			redNeighborsIndicies.clear();
		}
	}
	if (_ribs.first) enumerateRibs();
	return *this;
}

Triangulation& Triangulation::refine(Indicies& redList) {
	/*
		author: 
			Alexander Žilyakov, Jun 2016
		edited:
			Oct 2016
		comments:
			here we implement red–green refinement of our mesh
			@redList is a vector of indicies in _triangles to be red–refined
	*/
	std::unordered_map<Index, LocalIndex> greenMap;
	// hash table of indicies of triangles w/ hanging nodes
	// we must refine them green later
	// key is index in _triangles vector,
	// value is numb of hanging nodes in (1, 2, or 3)
	Indicies redNeighborsIndicies;
	// n <= 6 unknown neighbors of three new red triangles
	// that were added on previous iterations
	Index i;
	auto getRedNodeIndexFrom = [&](Index t) { // …from triangle _triangles[t]
		// we will need this function in order to 
		// organize red refinement if @redList contains neighbor triangles
		// because we do not want to add same nodes to _nodes vector several times
		LocalIndex m;
		for (m = 0; m < 3; ++m)
			if (getNeighborsIndicies(getNeighborsIndicies(t)[m])[1] != i &&
				getNeighborsIndicies(getNeighborsIndicies(t)[m])[2] != i) break;
		redNeighborsIndicies.push_back(getNeighborsIndicies(t)[nextIndex(m)]); // add red neighbors from
		redNeighborsIndicies.push_back(getNeighborsIndicies(t)[nextIndex(nextIndex(m))]); // previous refinements
		return getNodesIndicies(t)[m];
	};
	// RED PART
	Indicies::iterator redListIter = redList.begin(),
					   searchStartIter = prev(redList.end()), 
					   searchIter;

	// start with empty vector of fine neighbors
	_fineNeighborsIndicies = std::vector<std::vector<Index>>(numbOfElements());
	while (redListIter != redList.end()) {
		if (_ghostElements.find(i = *redListIter) != _ghostElements.end()) continue;
		// well, since this iteration
		// _elements[i] sholuld not ever be added to redList again
		greenMap[i] = 3; 
		// we will need ith nodes and neighbors later
		auto nodesIndicies = getNodesIndicies(i);
		auto neighborsIndicies = getNeighborsIndicies(i);
		// let’s deal w/ red points
		std::array<Index, 3> redNodesIndicies;
		Index k = 0; // numb of red nodes
		for (LocalIndex j : {0, 1, 2})
			if (0 <= neighborsIndicies[j] && !_makeNeighbors(i, neighborsIndicies[j])) {
				// if our triangle has a neighbor which has already been red-refined,
				// then no need in creating a new node
				// we just have to find out its index
				redNodesIndicies[j] = getRedNodeIndexFrom(neighborsIndicies[j]);
				++k; // how many nodes we do not need to create
			}
			else {
				redNodesIndicies[j] = numbOfNodes();
				if (neighborsIndicies[j] < -1) {
					//Index m = _neighbor2edge(t[j]);
					//_nodes.push_back(_curves[_curvilinearEdges[m].curveIndex()](_curvilinearEdges[m].thetaMiddle()));
					//_curvilinearEdges.push_back(CurvilinearEdge(_curvilinearEdges[m].thetaMiddle(), _curvilinearEdges[m].thetaEnd(), _curvilinearEdges[m].curveIndex()));
					//_curvilinearEdges[m].thetaEnd() = _curvilinearEdges[m].thetaMiddle();
				}
				else _nodes.push_back(midNodes(getElement(i))[j]);
			}	
		if (k > 1 && (searchIter = find(next(searchStartIter), redList.end(), i)) != redList.end()) 
			// if k equals 2 or 3
			// and we have also 3 red triangles yet to be added
			redList.erase(searchIter);
		// our ith triangle splits into 4 new ones
		// central one will take place of the old one
		_elements[i] = redNodesIndicies;
		_fineNeighborsIndicies[i] = { numbOfElements(), numbOfElements() + 1, numbOfElements() + 2 };
		SignedIndex dummy = numbOfElements();
		// and its neighbors are known (3 other triangles) and will be added soon
		_neighbors[i] = { dummy, dummy + 1, dummy + 2 };
		// now we have to add 3 other triangles
		// you have to draw them not to get confused w/ numeration
		// or just TRUST ME I AM A DOCTOR
		_elements.insert(_elements.end(), {
			{ nodesIndicies[0], redNodesIndicies[2], redNodesIndicies[1] },
			{ nodesIndicies[1], redNodesIndicies[0], redNodesIndicies[2] },
			{ nodesIndicies[2], redNodesIndicies[1], redNodesIndicies[0] }
		});
		dummy = i;
		_neighbors.insert(_neighbors.end(), {
			{ dummy, neighborsIndicies[1], neighborsIndicies[2] },
			{ dummy, neighborsIndicies[2], neighborsIndicies[0] },
			{ dummy, neighborsIndicies[0], neighborsIndicies[1] }
		});

		// fix curvilinear edges
		//j = _curvilinearEdges.size();
		//for (Index newTriangle : { _triangles.size() - 2, _triangles.size() - 3, _triangles.size() - 1 })
		//	if (_triangles[newTriangle].neighbors(1) < -1)
		//		_triangles[newTriangle].neighbors(1) = _neighbor2edge(--j);
		
		// we set neighbors of our new triangles to be neighbors of 
		// our old (refined) triangle
		// so there n <= 6 neighbors to be found
		// if ith neighbors were red-refined previously 
		for (Index j : { numbOfElements() - 3, numbOfElements() - 2, numbOfElements() - 1 })
			for (Index m : redNeighborsIndicies)
				_makeNeighbors(j, m);
		redNeighborsIndicies.clear(); // we are done w/ neighbors
		// if ith has no neighbors or its neighbors were refined earlier, we are gold 
		// otherwise we have to deal w/ hanging nodes
		for (LocalIndex j : {0, 1, 2})
			if (neighborsIndicies[j] > -1 && ++greenMap[neighborsIndicies[j]] == 2) 
					// if we have 2 (or 3 actually) hanging nodes,
					// we will not refine _triangle[t[j]] green
					redList.push_back(neighborsIndicies[j]); // we will refine it red instead!
		// we do not want to search from the very begining
		// because O(n^2) is too slow
		// so we use this smart hack
		if (redListIter++ == searchStartIter) searchStartIter = prev(redList.end());
	}
	// GREEN PART
	for (auto const & keyValue : greenMap) {
		if (keyValue.second > 1) continue; // we should do green refinement iff there’s only one hanging node
		i = keyValue.first; // index of triangle to be green-refined
		SignedIndex redTriangleIndex;
		LocalIndex j;
		for (j = 0; j < 3; ++j) {
			redTriangleIndex = getNeighborsIndicies(i)[j];
			if (0 <= redTriangleIndex && !_makeNeighbors(i, redTriangleIndex)) break; // so this is red triangle w/ a handing node
		}
		auto redNodeIndex = getRedNodeIndexFrom(redTriangleIndex); // so this is our hagning node
		// we want to split our ith triangle into two triangles
		// one of them we will add… 
		_elements.push_back({ getNodesIndicies(i)[j], getNodesIndicies(i)[nextIndex(j)], redNodeIndex });
		auto neighborIndex = getNeighborsIndicies(i)[nextIndex(nextIndex(j))];
		_neighbors.push_back({ -1, (SignedIndex)i, neighborIndex });
		// …(fix green neighbor)…
		if (neighborIndex > -1)
			_makeNeighbors(numbOfElements() - 1, neighborIndex);
		// …and another one will take place of the old one
		_fineNeighborsIndicies[i] = { numbOfElements() - 1 };
		_elements[i] = { getNodesIndicies(i)[j], redNodeIndex, getNodesIndicies(i)[nextIndex(nextIndex(j))] };
		_neighbors[i] = { -1, getNeighborsIndicies(i)[nextIndex(j)], (SignedIndex)numbOfElements() - 1 };
		// finally, lets fix neighbors
		_makeNeighbors(i, redNeighborsIndicies.front());
		_makeNeighbors(numbOfElements() - 1, redNeighborsIndicies.back());
		redNeighborsIndicies.clear();
	}
	if (_ribs.first) enumerateRibs();
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
	_neighbors = std::vector<std::array<SignedIndex, 3>>(numbOfElements(), { -1, -1, -1 });
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
	_ribs = { currentRibIndex, numeration };
	return *this;
}

Node2D Triangulation::getRibNode(Index e, LocalIndex r, double t) const {
	SignedIndex edgeIndex = _neighbor2edge(getNeighborsIndicies(e)[r]);
	if (edgeIndex < 0) {
		auto
			A  = getElement(e)[nextIndex(r)],
			AB = getElement(e)[nextIndex(nextIndex(r))] - A;
		return A + t * AB;
	}
	auto s = t * (_edges[edgeIndex].thetaEnd() - _edges[edgeIndex].thetaStart());
	return _curves[_edges[edgeIndex].curveIndex()](_edges[edgeIndex].thetaStart() + s);
}