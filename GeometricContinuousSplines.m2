newPackage(
    "GeometricContinuousSplines",
    Version => "0.1",
    Date => "28 March 2023",
    Headline => "This is a package for geometrically continuous splines.",
    Authors => {{ Name => "Beihui Yuan", Email => "by238@cornell.edu", HomePage => "https://sites.google.com/view/beihuiyuan/home"}},
    AuxiliaryFiles => false,
    DebuggingMode => false,
    PackageImports => {"RationalPoints2"}
    )

export {
    "computeBaseField",
    "generateAmbientRing",
    "gSplineBasis",
    "monomialBasisBiDegree",
    "monomialBasisT"
    }

-* Code section *-
---------------------------
computeBaseField = method()
---------------------------
---------------------------
--This method compute a field extension of QQ such that 
--the transition map can be defined over. 
---------------------------
--Inputs:
---------------------------
--valences = list of valences. Each valence is a positive integer.
---------------------------
--Outputs:
---------------------------
--A field of extension over QQ.
--usage: computeBaseField(valences)
---------------------------
computeBaseField(List):=(valences) ->(
    n := lcm(valences);
    x := symbol x; 
    QQ[x];
    f := sum for i from 0 to n-1 list x^i;
    F := extField(QQ,f)
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
generateAmbientRing(List,Ring):=(E,kk) ->(
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
gSplineBasis(List,List,ZZ) := (E,ideals,deg) ->(
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
TEST /// -* The 2-patches case *-
-- test code and assertions here
-- may have as many TEST sections as needed
E = {{1,2}};
R = generateAmbientRing(E,QQ)
r=1;
mathfraka=(f) -> 2*f-1;--(n_0,n_1)=(3,3)
mathfraka=(f)-> f^2; -- (n_0,n_1)=(4,3)
mathfraka=(f)->-f^2+2*f-1; --(n_0,n_1)=(3,4)
--mathfraka=(f)-> 0--(n_0,n_1)=(4,4)
gtm=(uminus,vminus,uplus,vplus)->{uminus+vplus, vminus-(uplus+vplus*mathfraka(uplus))};
J = ideal gtm(R_0,R_1,R_2,R_3);
I = J+ideal(R_0^(r+1),R_3^(r+1));
gSplineBasis(E,{I},2)
monomialBasisT(E,R,3)
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
gSplineBasis(E, for e in E list hTIdeals#e, d)
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
