newPackage(
    "GeometricContinuousSplines",
    Version => "0.1",
    Date => "28 March 2023",
    Headline => "This is a package for geometrically continuous splines.",
    Authors => {{ Name => "Beihui Yuan", Email => "by238@cornell.edu", HomePage => "https://sites.google.com/view/beihuiyuan/home"}},
    AuxiliaryFiles => false,
    DebuggingMode => false
    )

export {
    "generateAmbientRing",
    "gSplineBasis"
    }

-* Code section *-
---------------------------
generateAmbientRing = method()
---------------------------
---------------------------
--This method creates an ambient ring for G spline space. 
---------------------------
--Inputs:
---------------------------
--E = list of edges. Each edge is a list with two vertices
---------------------------
--Outputs:
---------------------------
--A polynomial ring QQ[u_sigma,v_sigma:sigma in vert]
--usage: generateAmbientRing(E)
---------------------------
generateAmbientRing(List):=(E) ->(
    vert := sort unique flatten E;
    u := symbol u;
    v := symbol v;
    S := QQ[flatten for sigma in vert list {u_sigma,v_sigma}];
    S
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
--ideals = list of ideals that label the edges
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
TEST /// -* [insert short title for this test] *-
-- test code and assertions here
-- may have as many TEST sections as needed
E = {{1,2}};
R = generateAmbientRing(E)
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
