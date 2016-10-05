# Δ **P**<sup>1</sup> CR — Δ P<sup>0</sup> L Interpolant Visualization for Velocity Field and Pressure Distribution

This is how current project is organized:

1. First, we generate symbolic mesh and convert it to CATSPDEs format w/ [`./Mathematica/Generate Mesh/generateMesh.nb`](https://github.com/CATSPDEs/CATSPDEs/tree/master/sln/%CE%94%20P1%20CR%20%E2%80%94%20%CE%94%20P0%20L%20Interpolant%20Visualization/Mathematica/Generate Mesh/generateMesh.nb).

2. Then we import generated mesh w/ cpp unit, [`./user.cpp`](https://github.com/CATSPDEs/CATSPDEs/tree/master/sln/%CE%94%20P1%20CR%20%E2%80%94%20%CE%94%20P0%20L%20Interpolant%20Visualization/user.cpp), (uniformly) refine it if necessary, and add ribs’ (middle nodes’) numeration (i.e. convert mesh from “Nodes and Triangles” to “Nodes, Simple Ribs, and Triangles” format).
   
   One needs to enumerate ribs since their numeration, in fact, is a numeration of DOFs for Crouzeix–Raviart finite elements (Δ P<sup>1</sup> CR). We will use them for velocity components.
   
   DOFs enumeration for Δ P<sup>0</sup> L finite elements (which are used for pressure distribution) is simple and is induced by triangles’ numeration.
   
   [`./user.cpp`](https://github.com/CATSPDEs/CATSPDEs/tree/master/sln/%CE%94%20P1%20CR%20%E2%80%94%20%CE%94%20P0%20L%20Interpolant%20Visualization/user.cpp) also provides vectors `u1Vec`, `u2Vec`, and `pVec` for further processing.

3. [`./Mathematica/Draw Interpolant/drawInterpolant.nb`](https://github.com/CATSPDEs/CATSPDEs/tree/master/sln/%CE%94%20P1%20CR%20%E2%80%94%20%CE%94%20P0%20L%20Interpolant%20Visualization/Mathematica/Draw Interpolant/drawInterpolant.nb) then loads resulting from (1 – 2) mesh and vectors.

   It provides all post–processing stuff:
   * visualazing resulting mesh,
   * drawing FE interpolants *P* [ `u1Vec` ], *P* [`u2Vec`], *P* [`pVec`] (also as a vector velocity field w/ background density plot for pressure distribution),
   * computing L<sub>2</sub>–norms of errors,
   * etc.
   
## Some Post–processing Results

![ΔP1CR shape funcs](https://github.com/CATSPDEs/CATSPDEs/tree/master/sln/%CE%94%20P1%20CR%20%E2%80%94%20%CE%94%20P0%20L%20Interpolant%20Visualization/img/ΔP1CRshapes.png)

![ΔP0L shape func](https://github.com/CATSPDEs/CATSPDEs/tree/master/sln/%CE%94%20P1%20CR%20%E2%80%94%20%CE%94%20P0%20L%20Interpolant%20Visualization/img/ΔP0Lshape.png)
	 
— Žilyzkov Alexander, Oct 2016
