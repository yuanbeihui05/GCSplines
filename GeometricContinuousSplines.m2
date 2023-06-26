newPackage(
    "GeometricContinuousSplines",
    Version => "0.1",
    Date => "4 May 2023",
    Headline => "This is a package for geometrically continuous splines.",
    Authors => {{ Name => "Beihui Yuan", Email => "by238@cornell.edu", HomePage => "https://sites.google.com/view/beihuiyuan/home"}},
    AuxiliaryFiles => false,
    DebuggingMode => false,
    PackageExports => {"Cyclotomic"}
    )

export {
    "computeBaseField",
    "createPyFile",
    "createStarVertexPatch",
    "exportXYZcoeff",
    "generateAmbientRing",
    "gSplineBasis",
    "monomialBasisBiDegree",
    "monomialBasisT"
    }

-* Todo *-
-- Creating an external file for visualization (05/10/2023)--
-- Give G splines a Type (?)--
-- Field Extenstion --
-- Symmetric gluing data --


-* Code section *-
---------------------------
createPyFile = method()
---------------------------
---------------------------
--This method creates .py file for visualization
--At present, it is only for square patches 
---------------------------
--Inputs:
---------------------------
-- mP = a (3xn) matrix, 
--      giving the x,y,z coordinates, 
--      each entry is a polynomial
-- deg = a non-negative integer, the degree of the spline space 
-- uvrange = a list of two numbers.
-- fileName = a string, the name of output file
---------------------------
--Outputs:
---------------------------
--A .py file
--usage: createPyFile(mP, fileName)
---------------------------
createPyFile(Matrix, ZZ, List, String) := (mP,deg,uvrange,fileName) -> (
    S := ring mP;
    nump := numColumns mP;
    f := concatenate{fileName, ".py"} << ""; --initial a file
    f << "import numpy as np" << endl;
    f << "import pyvista as pv" << endl;
    f << "" << endl;
    f << "xyz_poly_coefs = (" << endl;
    for k from 0 to (nump-1) do(
        f << "                  np.array([" << endl;
        for polynf in entries mP_k do(
            coeffm := for i from 0 to deg list
            for j from 0 to deg list coefficient(S_(2*k)^i*S_(2*k+1)^j,polynf);
            stringPolycoeff := replace("{","[",toString coeffm);
            stringPolycoeff = replace("}","]",stringPolycoeff);
            f << "                    " << stringPolycoeff << "," << endl;
        );
        f << "                    ]).transpose(1,2,0)," << endl;
        f << "" << endl;
    );
    f << "                  )" << endl;
    stringuvrange := replace("{","[", toString uvrange);
    stringuvrange = replace("}","]",stringuvrange);
    f << "uv_ranges = [[" << stringuvrange << "," << stringuvrange << "]]*" << nump << endl;
    f << "" << endl;
    f << "if __name__ == '__main__':" << endl;
    f << "    pv.set_plot_theme('paraview')" << endl;
    f << "    plotter = pv.Plotter()" << endl;
    f << "    plotter.set_color_cycler([" << endl;
    f << "        \"#e60049\", \"#0bb4ff\", \"#50e991\", \"#e6d800\", \"#9b19f5\", " << endl;
    f << "        \"#ffa300\", \"#dc0ab4\", \"#b3d4ff\", \"#00bfa0\"" << endl;
    f << "        ])" << endl;
    --f << "    plotter.set_color_cycler(['magenta', 'seagreen', 'aqua', 'orange'])" << endl;
    f << "    plotter.show_axes_all()" << endl;
    f << "    plotter.show_grid()" << endl;
    f << "    plotter.enable_anti_aliasing()" << endl;
    f << "" << endl;
    f << "    for (u_range,v_range),coefs in zip(uv_ranges, xyz_poly_coefs): " << endl;
    f << "        u = np.linspace(u_range[0], u_range[1], 20)" << endl;
    f << "        v = np.linspace(v_range[0], v_range[1], 20)" << endl;
    f << "        value = np.polynomial.polynomial.polygrid2d(u, v, coefs)" << endl;
    f << "" << endl;
    f << "        mesh = pv.StructuredGrid(value[0], value[1], value[2])   " << endl;
    f << "        plotter.add_mesh(mesh, smooth_shading=True, split_sharp_edges=True)" << endl;
    f << "" << endl;
    f << "    plotter.show()" << endl << close;
    f
)

---------------------------
computeBaseField = method()
---------------------------
---------------------------
--This method compute a field extension of QQ such that 
--the transition map can be defined over. 
---------------------------
--Inputs: 
---------------------------
--valences = list of positive integers.
---------------------------
--Outputs:
---------------------------
--A field of extension over QQ.
--usage: computeBaseField(valences)
---------------------------
computeBaseField(List):= Ring =>(valences) ->(
    n := lcm(valences);
    KK0 := cyclotomicField(n);
    ksi := KK0_0;
    e := symbol e;
    R := QQ[for i in valences list e_i];
    phi := map(KK0, R, matrix {for i in valences list ksi^(sub(n/i,ZZ))/2+ksi^(sub(n-n/i,ZZ))/2});
    KK := R/ker phi;
    KK
    )
---------------------------
createStarVertexPatch = method()--not done yet
---------------------------
---------------------------
--This method compute a basis for G spline spaces over a star of a vertex, 
--assuming that transition maps are from symmetric gluing data. 
---------------------------
--Inputs: 
---------------------------
--valences = list of positive integers. May have repetition.
-------------input the valences at each BOUNDARY vertex clockwisely 
--deg = a non-negative integer, the degree bound for the spline space
---------------------------
--Outputs:
---------------------------
--A nested list. A basis for the spline space over base filed RR
--usage: createStarVertexPatch(valences, degree)
---------------------------
---------------------------
--Function dependence:
---------------------------
-- gSplineBasis
-- computeBaseField
-- generateAmbientRing
---------------------------
createStarVertexPatch(List, ZZ):= (valences, deg) ->(
    n := #valences;
    E := for i from 1 to n-1 list {i,i+1};
    E = {{n,1}}|E;
    uniValences := unique ({n}|valences);
    --uniValences = {n}|uniValences;
    posValences := for m in ({n}|valences) list position(uniValences,i->i==m);
    KK := computeBaseField(uniValences);
    S := generateAmbientRing(E,KK);
    SRR := generateAmbientRing(E,RR);
    gluea := (f, i)->(2*KK_(posValences_0)*(1-f)^2-2*KK_(posValences_i)*f^2); 
    glueb := -1;
    um := S_(2*n-2);
    vm := S_(2*n-1);
    up := S_0;
    vp := S_1;
    ideals := {ideal (um^2,vp^2,um+vp,vm-up-vp*gluea(up,n))};
    for i from 0 to n-2 do (
        um = S_(2*i);
        vm = S_(2*i+1);
        up = S_(2*i+2);
        vp = S_(2*i+3);
        ideals = ideals|{ideal(um^2,vp^2,um+vp,vm-up-vp*gluea(up,i+1))};);
    starPatchVertexBasis := (gSplineBasis(E,ideals,deg))#"basis";
    sub(starPatchVertexBasis, 
    (for i from 0 to #uniValences-1 list (KK_i => cos(2*pi/uniValences_i)))|
    for j from 0 to 2*n-1 list (S_j => SRR_j))
    )---complete on June 26, 2023.

---------------------------
exportXYZcoeff = method()
---------------------------
---------------------------
--This method creates .txt file containing coefficient matrices
---------------------------
--Inputs:
---------------------------
-- mP = a (3xn) matrix, 
--      giving the x,y,z coordinates, 
--      each entry is a polynomial.
-- uvrange = a list of two numbers.
-- fileName = a string, the name of output file.
---------------------------
--Outputs:
---------------------------
--A .txt file
--usage: createPyFile(mP, fileName)
---------------------------
exportXYZcoeff(Matrix, List, String) := (mP, uvrange, fileName) -> (
    S := ring mP;
    nump := numColumns mP;--number of patches
    f := concatenate{fileName, ".txt"} << ""; --initial a file
    polycoef := for k from 0 to (nump-1) list
      for polynf in entries mP_k list
      for i from 0 to degree(S_(2*k),polynf) list
      for j from 0 to degree(S_(2*k+1),polynf) list coefficient(S_(2*k)^i*S_(2*k+1)^j,polynf); 
    f << nump << endl;
    stringPolycoeff := replace("{","[",toString polycoef);
    stringPolycoeff = replace("}","]",stringPolycoeff);
    f << stringPolycoeff << endl;
    f << toString uvrange << endl << close;
    f
)

---------------------------
generateAmbientRing = method()
---------------------------
---------------------------
--This method creates an ambient ring for G spline space. 
---------------------------
--Inputs:
---------------------------
--E = list of edges. Each edge is a list with two vertices
-- kk = base field
---------------------------
--Outputs:
---------------------------
--A polynomial ring kk[u_sigma,v_sigma:sigma in vert]
--usage: generateAmbientRing(E)
---------------------------
generateAmbientRing(List,Ring):= Ring => (E,kk) ->(
    vert := sort unique flatten E;
    u := symbol u;
    v := symbol v;
    S := kk[flatten for sigma in vert list {u_sigma,v_sigma}];
    S
    )

generateAmbientRing(List,InexactFieldFamily):= Ring => (E,kk) ->(
    vert := sort unique flatten E;
    u := symbol u;
    v := symbol v;
    S := kk[flatten for sigma in vert list {u_sigma,v_sigma}];
    S
    )

---------------------------
monomialBasisBiDegree = method()
---------------------------
--This method creates a list of monomial basis (for bi-degree)
--for each coordinate ring with bi-degree no more than (d,d). 
---------------------------
--Inputs:
---------------------------
--E = list of edges. Each edge is a list with two vertices
--S = ambient ring
--deg = bi-degree
---------------------------
--Outputs:
---------------------------
--A hash table sigma => a list of monomials over vertex sigma
--usage: monomialBasisBiDegree(E,S,deg)
---------------------------
monomialBasisBiDegree(List,Ring,ZZ) := (E,S,deg) ->(
    vert := sort unique flatten E;
    monB:= hashTable for sigma in vert list 
    sigma => flatten entries monomials((S_(2*sigma-2)+1)^deg*(S_(2*sigma-1)+1)^deg);--flatten for i from 0 to deg list for j from 0 to deg-i list varsList_(2*sigma-2)^i*varsList_(2*sigma-1)^j;
    monB
    )

---------------------------
monomialBasisT = method()
---------------------------
--This method creates a list of monomial basis (for total degree)
--for each coordinate ring with total degree no more than d. 
---------------------------
--Inputs:
---------------------------
--E = list of edges. Each edge is a list with two vertices
--S = ambient ring
--deg = the total degree
---------------------------
--Outputs:
---------------------------
--A hash table sigma => a list of monomials over vertex sigma
--usage: monomialBasisT(E,S,deg)
---------------------------
monomialBasisT(List,Ring,ZZ) := (E,S,deg) ->(
    vert := sort unique flatten E;
    monB:= hashTable for sigma in vert list 
    sigma => flatten entries monomials((S_(2*sigma-2)+S_(2*sigma-1)+1)^deg);--flatten for i from 0 to deg list for j from 0 to deg-i list varsList_(2*sigma-2)^i*varsList_(2*sigma-1)^j;
    monB
    )

---------------------------
gSplineBasis = method()
---------------------------
---------------------------
--This method computes the geometrically continuous spline spaces 
--associated with a graph whose edges are labeled by ideals. 
--Also see generalizedSplines in AlgebraicSplines.m2 package.
---------------------------
--Inputs:
---------------------------
--E = list of edges. Each edge is a list with two vertices
--ideals = list of ideals that labeled by the edges
--deg = an integer, which is the degree of the spline space
---------------------------
--Outputs:
---------------------------
--A Hash Table of dimension and basis in the given geometrically continuous spline spaces, 
--up to the given degree
--usage: gSplineBasis(E, L, d)#"dimension", gSplineBasis(E,L,d)#"basis"
---------------------------
gSplineBasis(List,List,ZZ) := HashTable => (E,ideals,deg) ->(
    vert := sort unique flatten E;
    S:= ring first ideals;
    --one has to input the underlying ring before using this function
    --make sure ideals all lie in the same ring
    ideals = apply(ideals, I->sub(I,S));
    --Possible future development: check if the generators of ideals are
    ------------------in variables of the two faces, using function support(f)
    --hashTable E=>ideals, label ideals with the correponding edge
    labelIdeals := hashTable for i from 0 to #E-1 list E#i=>ideals#i;
    --Boundary map from Edges to vertices, 
    --its rows are labeled by edges and columns are labeled by vertices
    boundaryEV:= matrix apply(E,
	e->apply(vert,
	    sigma->if(sigma===first e) then 1
	    else if(sigma===last e) then -1
	    else 0));
    boundaryEV = sub(boundaryEV,S);
    --Generate a monomial basis for each vertex sigma and degree no more than deg
    --monB is a hash table sigma => a list of monomial basis over sigma
    varsList := flatten entries vars S;
    monB:= hashTable for sigma in vert list 
    sigma => flatten for i from 0 to deg list for j from 0 to deg-i list varsList_(2*sigma-2)^i*varsList_(2*sigma-1)^j;
    --hashTable remainders, {e,sigma, m} => remainder
    hTrem:= hashTable flatten flatten for e in E list
    for sigma in vert list
    for m in monB#sigma list {e,sigma,m}=>
    (if (sigma===first e) then m%(labelIdeals#e)
    else if (sigma===last e) then -m%(labelIdeals#e)
    else 0);
    --generate a matrix corresponding to delta_2 from the hashTable hTEvm
    --of which rows labeled by edges, columns labeled by faces
    bdryM := matrix for e in E list {
	(coefficients matrix{flatten for sigma in vert list 
		for m in monB#sigma list hTrem#{e,sigma,m}})_1
	};
    --compute the kernel of bdryM, denoted by kerBdryM.
    --columns of which form a basis to the spline space.
    kerBdryM := ker bdryM;
    dimBasis := numColumns gens kerBdryM;
    sourceB := directSum(for sigma in vert list
	matrix {for m in monB#sigma list m});
    GSplineSpaceB := sourceB * gens kerBdryM;
    --output: hashTable of dimension and basis
    hTresults := new HashTable from { 
	"dimension" => dimBasis, 
	"basis" => GSplineSpaceB};
    hTresults
    --coding completes, yet to be tested. --Mar. 12, 2023
)


-* Documentation section *-
beginDocumentation()

doc ///
Key
  GeometricContinuousSplines
Headline
 A package for G splines
Description
  Text
   Still trying to figure this out.
  --Example
  --CannedExample
Acknowledgement
Contributors
References
Caveat
SeeAlso
Subnodes
///

doc ///
Key
 gSplineBasis
Headline
 Headline for gSplineBasis
Usage
 gSplineBasis(edges,ideals,deg)
Inputs
 edges:List
     list of edges
 ideals:List
     list of ideals
 deg:ZZ
     desired degree
Outputs
 LB:List
      a list of basis
--Consequences
--  Item
Description
  Text
      This method works for...
  --Example
  --CannedExample
  --Code
  --Pre
--ExampleFiles
--Contributors
--References
--Caveat
SeeAlso
///

-* Test section *-
-- test code and assertions here
-- may have as many TEST sections as needed
TEST /// -* The 2-patches case *-
E = {{1,2}};
R = generateAmbientRing(E,QQ);
r=1;
mathfraka=(f) -> 2*f-1;--(n_0,n_1)=(3,3)
mathfraka=(f)-> f^2; -- (n_0,n_1)=(4,3)
mathfraka=(f)->-f^2+2*f-1; --(n_0,n_1)=(3,4)
--mathfraka=(f)-> 0--(n_0,n_1)=(4,4)
gtm=(uminus,vminus,uplus,vplus)->{uminus+vplus, vminus-(uplus+vplus*mathfraka(uplus))};
J = ideal gtm(R_0,R_1,R_2,R_3);
I = J+ideal(R_0^(r+1),R_3^(r+1));
d = 4;
hTtwopatches = gSplineBasis(E,{I}, d);
monomialBasisT(E,R,3)
-- test exportXYZcoeff and createPyFile --
mBasis = hTtwopatches#"basis";
dimS = hTtwopatches#"dimension";
mP = transpose (mBasis*random(QQ^dimS,QQ^3));
--exportXYZcoeff(mP,{0,1},"2patches")
createPyFile(mP,d,{0,1},"2patches")
///

TEST ///-* a star of a vertex (3 faces) *-
E = {{1,2},{2,3},{3,1}};
R = generateAmbientRing(E,QQ);
r=1;
mathfraka = (f)->-f^2+2*f-1; 
mathfrakb = (f) -> -1;
Lgtm = (uplus,vplus,uminus,vminus)->({uminus-mathfrakb(uplus)*vplus, vminus-uplus-vplus*mathfraka(uplus)});
generateItau = (indexface1,indexface2)->(
    ideal Lgtm(R_(indexface1*2-2),R_(indexface1*2-1),R_(indexface2*2-2),R_(indexface2*2-1))+ideal(R_(indexface1*2-1)^(r+1),R_(indexface2*2-2)^(r+1))
    );
hTIdeals = hashTable for e in E list e=>generateItau(first e, last e);
d = 4;
hTthreeSplit = gSplineBasis(E, for e in E list hTIdeals#e, d);
--interpolation
mBasis = hTthreeSplit#"basis";
dimS = hTthreeSplit#"dimension";
mP = transpose (mBasis*random(QQ^dimS,QQ^3));
createPyFile(mP,d,{0,1},"starVertex")
///

TEST ///-* The cube case *-
---Step 1: generating the corresponding cell complex, in this case a cube--
needsPackage "SimplicialComplexes"
numfaces = 6;
positionfaces={{0,0,1},{1,0,0},{0,1,0},{-1,0,0},{0,-1,0},{0,0,-1}};
T=ZZ[sigma_1..sigma_numfaces,Degrees=>positionfaces];
Octahedral=simplicialComplex monomialIdeal(sigma_1*sigma_6,sigma_2*sigma_4,sigma_3*sigma_5);
facesOfCube = (faces(Octahedral))#0;
edgesOfOct = (faces(Octahedral))#1;
verticesOfCube = (faces(Octahedral))#2;
E = for m in edgesOfOct list (indices m)+{1,1};
vert = unique flatten E;

--Step 2: generating the list of ideals--
r=1;
R=generateAmbientRing(E,QQ);
basisOnChart={{R_0,R_1,1},{1,R_2,R_3},{R_4,1,R_5},{-1,R_6,R_7},{R_8,-1,R_9},{R_10,R_11,-1}};

mathfraka = (f)->(2*f-1);--(n_0,n_1)=(3,3)
mathfrakb = (f)-> -1;
Lgtm = (uplus,vplus,uminus,vminus)->({uminus-mathfrakb(uplus)*vplus, vminus-uplus-vplus*mathfraka(uplus)});

endpoints = (indexface1,indexface2) -> (faces star(Octahedral,sigma_indexface1*sigma_indexface2))#2
theOtherEndPoint = (indexface1, indexface2, indexface3) -> (
    if member(sigma_indexface1*sigma_indexface2*sigma_indexface3,endpoints(indexface1,indexface2)) then 
    for x in support(sum endpoints(indexface1,indexface2)-sigma_indexface1*sigma_indexface2*sigma_indexface3) list index(x)+1)  

orientationOfFaces = (indexface1, indexface2, indexface3)->(det matrix{degree sigma_indexface1, degree sigma_indexface2, degree sigma_indexface3})
putIntoClockwise = (indexface1, indexface2, indexface3)->(
    if orientationOfFaces(indexface1,indexface2,indexface3) == 1 then {indexface2, indexface1, indexface3}
    else if orientationOfFaces(indexface1,indexface2,indexface3) == -1 then {indexface1,indexface2,indexface3}
    else "N/A"
    )
verticesOfFace = (indexface) -> (faces star (Octahedral,sigma_indexface))#2

coordPlus = (indexface1,indexface2,indexface3)->(
    Lindexfaces=putIntoClockwise(indexface1,indexface2,indexface3);
    coordp00=product for x in Lindexfaces list sigma_x;
    coordp10=product for x in theOtherEndPoint(Lindexfaces_0,Lindexfaces_1,Lindexfaces_2) list sigma_x;
    coordp01=product for x in theOtherEndPoint(Lindexfaces_0,Lindexfaces_2,Lindexfaces_1) list sigma_x;
    coordp11= (flatten delete(coordp00,delete(coordp01,delete(coordp10,verticesOfFace(Lindexfaces_0)))))_0;
    matrixcoordp =promote(matrix {degree coordp00, degree coordp01, degree coordp10, degree coordp11},QQ);
    coordup = solve(matrixcoordp,promote(matrix{{0},{0},{1},{1}},QQ));
    coordvp = solve(matrixcoordp,promote(matrix{{0},{1},{0},{1}},QQ));
    transpose matrix{flatten entries coordup, flatten entries coordvp}
    )
coordMinus = (indexface1,indexface2,indexface3)->(
    Lindexfaces=putIntoClockwise(indexface1,indexface2,indexface3);
    coordm00=product for x in Lindexfaces list sigma_x;
    coordm01=product for x in theOtherEndPoint(Lindexfaces_1,Lindexfaces_0,Lindexfaces_2) list sigma_x;
    coordm10=product for x in theOtherEndPoint(Lindexfaces_1,Lindexfaces_2,Lindexfaces_0) list sigma_x;
    coordm11= (flatten delete(coordm00,delete(coordm01,delete(coordm10,verticesOfFace(Lindexfaces_1)))))_0;
    matrixcoordm =promote(matrix {degree coordm00, degree coordm01, degree coordm10, degree coordm11},QQ);
    coordum = solve(matrixcoordm,promote(matrix{{0},{0},{1},{1}},QQ));
    coordvm = solve(matrixcoordm,promote(matrix{{0},{1},{0},{1}},QQ));
    transpose matrix{flatten entries coordum, flatten entries coordvm}
    )
affineTransformation = (indexface1,indexface2,indexface3)->(
    Lindexfaces=putIntoClockwise(indexface1,indexface2,indexface3);
    indexp = Lindexfaces_0;
    indexm = Lindexfaces_1;
    pmatrix = matrix {basisOnChart_(indexp-1)}*coordPlus(indexface1,indexface2,indexface3);
    mmatrix = matrix {basisOnChart_(indexm-1)}*coordMinus(indexface1,indexface2,indexface3);
    {flatten entries pmatrix, flatten entries mmatrix}--output {{up,vp},{um,vm}}
    )

generateItau = (indexface1,indexface2)->(
    endpointsL = endpoints(indexface1,indexface2)/(sigma_indexface1*sigma_indexface2);
    indexf3 = (index endpointsL_0)+1;
    indexf4 = (index endpointsL_1)+1;
    matrixOfuvAround123=matrix affineTransformation(indexface1,indexface2, indexf3);
    up123 = matrixOfuvAround123_(0,0);
    vp123 = matrixOfuvAround123_(0,1);
    um123 = matrixOfuvAround123_(1,0);
    vm123 = matrixOfuvAround123_(1,1);
    matrixOfuvAround124=matrix affineTransformation(indexface1,indexface2, indexf4);
    up124 = matrixOfuvAround124_(0,0);
    vp124 = matrixOfuvAround124_(0,1);
    um124 = matrixOfuvAround124_(1,0);
    vm124 = matrixOfuvAround124_(1,1);
    ideal Lgtm(up123,vp123,um123,vm123)+ideal Lgtm(up124,vp124,um124,vm124)+ideal(um123^(r+1),vp123^(r+1),um124^(r+1),vp124^(r+1))
    )
--hashTable edge => ideals
hTIdeals = hashTable for e in E list e=>generateItau(first e, last e);
d = 4;
hTcube = gSplineBasis(E, for e in E list hTIdeals#e, d);
--interpolation
mBasis = hTcube#"basis";
mP = mBasis*inverse sub(mBasis, for i from 0 to 11 list R_i=>0)*matrix{{0,0,1/1},{1,0,0},{0,1,0},{-1,0,0},{0,-1,0},{0,0,-1}};
mP = transpose mP;
exportXYZcoeff(mP,{-1,1},"cube")
///
end--

-* Development section *-
restart
path = append(path,"./")
debug needsPackage "GeometricContinuousSplines"
check "GeometricContinuousSplines"

uninstallPackage "GeometricContinuousSplines"
restart
installPackage "GeometricContinuousSplines"
viewHelp "GeometricContinuousSplines"

createStarVertexPatch({4,4,4},3)
createStarVertexPatch({3,5,4,5},4)
createStarVertexPatch({4,4,4,4,4},4)
