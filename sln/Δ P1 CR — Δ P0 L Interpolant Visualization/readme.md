This is how current project is organized:

(1) First, we generate symbolic mesh and convert it to CATSPDEs format w/ 
    ./Mathematica/Generate Mesh/generateMesh.nb

(2) Then we import generated mesh w/ this unit, ./user.cpp, (uniformly) refine it if necessary, add ribs’ (mid nodes’) 
    numeration (i.e. convert mesh from “Nodes and Triangles” to “Nodes, Simple Ribs, and Triangles” format).
	One needs to enumerate ribs since their numeration is, in fact, numeration of DOFs for Crouzeix–Raviart 
	finite elements (Δ P1 CR). We will use them for velocity components.
	DOFs enumeration for Δ P0 L finite elements (which are used for pressure distribution) are simple and is 
	induced by triangles numeration.
	./user.cpp also provides vectors u1Vec, u2Vec, and pVec for further processing.

(3) ./Mathematica/Draw Interpolant/drawInterpolant.nb then loads resulting from (2) mesh and vectors.
    It provides all post–processing stuff:
	* visualazing resulting mesh,
	* drawing FE interpolants P [ u1Vec ], P [u2Vec], P [pVec] 
	  (also as a vector velocity field w/ background density plot for pressure distribution),
	* computing L2–norms of errors,
	* etc.
	 
— Žilyzkov Alexander, Oct 2016