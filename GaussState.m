(* ::Package:: *)

(* Mathematica Raw Program 

add 
- CV-to-Fock-state function
- Gaussian quantum discord 1.2358   -0.0868   -1.2209   -0.1559
   -0.0868    1.3874   -0.0844    1.3874
   -1.2209   -


*)

(* 
Author: Jason Hoelscher-Obermaier
Last Modified: 2012-10-12

This package contains tools for the calculation of covariance matrices of Gaussian states - 
as well as the evaluation of certain entanglement measures on the level of these covariance matrices.
*)

BeginPackage["GaussState`",{"GraphUtilities`"}]

Print["Type \"?aboutGaussState\ for general info about the package and \"?GAUSSSTATE`GaussState`*\" for a list of the built-in functions. \n
By default, the mode operator convention is the canonical with [a, \!\(\*SuperscriptBox[\"a\", \"\[Dagger]\"]\)]=1 and [q,p]=\[ImaginaryI].\n
 But it can be changed to the quantum-optical convention by changing the value of the variable \[Kappa]1, \[Kappa]2 and \[Kappa]3 
at the beginning of the private environment of the package file GaussState.m

Ferraro et al. [2005] provide an excellent reference for this purpose: They list many formulas explicitly as functions of the convention-dependent variables \[Kappa] 1 , \[Kappa] 2 , \[Kappa] 3
which are defined in [Ferraro et al., section 1.5.2].

Reference:
Ferraro, A., Olivares, S., & Paris, M. G. A. [2005, Mar.] Gaussian states in continuous variable
quantum information. quant-ph/0503237."];

(* Begins USAGE description of Commands *)
aboutGaussState::usage = 
" This package contains tools for the calculation of covariance matrices of Gaussian states - 
\n as well as the evaluation of certain entanglement measures on the level of these covariance matrices. 
\n All commands are named in a Mathematica-like fashion, except for a \"gs\" in front of every command. "

(*  ================================================================= *)
(*  START USAGE DOCUMENTATION  *)
gsDirectSum::usage=
" gsDirectSum[m1,m2] returns the direct sum of two matrices m1 and m2: m1\[CirclePlus]m2.
\n The direct sum of two matrices is the matrix representation of the direct sum of the
\n linear operators corresponding to m1 and m2. "

gsVacuumState::usage=
" gsVacuumState[n] returns the covariance matrix of the vacuum state of n modes."
gsVac::usage=
" gsVacuumState[n] returns the covariance matrix of the vacuum state of n modes."

gsSqueezedState::usage=
" gsSqueezedState[N,{r_1, ..., r_N}] returns the covariance matrix of N modes, 
\n with the i-th mode subject to momentum squeezing by the factor r_i."
gsTMSS::usage=
"gsTMSS[r] returns the two-mode squeezed state with squeezing parameter r."

gsStandardFormTMGS::usage=
"gsStandardFormTMGS[a,b,cp,cm] returns the standard form of a two-mode Gaussian state with parameters a,b,cp,cm.
If no other values are assigned to a,b,cp,cm these values are chosen by default."

gsStandardFormTMGSParameters::usage=
"gsStandardFormTMGSParameters[M] returns the parameters {a, b, cp, cm} of the standard form (gsStandardFormTMGS)of a two-mode Gaussian state\n
described by the covariance matrix M."
gsSFParams::usage=
"gsStandardFormTMGSParameters[M] returns the parameters {a, b, cp, cm} of the standard form (gsStandardFormTMGS)of a two-mode Gaussian state\n
described by the covariance matrix M."

gsSymplecticForm::usage=
" gsSymplecticForm[N] returns the symplectic form for a N-mode system."

gsSymplecticVector::usage=
" gsSymplecticVector[N] returns a vector (q1,p1,...,qn,pn) containing the 2N conjugate variable of a N-mode system."

gsWignerFunction::usage=
" gsWignerFunction[M] returns the Wigner function of the state corresponding to the covariance matrix M.
\n The variables are called q1,p1,....,qN,pN. "
gsSymplecticTrafo::usage=
" gsSymplecticTrafo[N,{{R_i,R_j,H_{i,j}},...}] returns the symplectic trafo corresponding to the unitary trafo
\n U=exp[\[ImaginaryI]/2 \[Sum]_{k,l} H_{k,l} R_k R_l], hence to a homogeneous quadratic Hamiltonian in the mode operators."

gsModeInterchange::usage=
"gsModeInterchange[M,m1,m2] returns the covariance matrix M after interchange of the modes m1 and m2." 

gsQPOrdering::usage=
"gsQPOrdering[M] returns the covariance matrix in the ordering (\!\(\*SubscriptBox[\"Q\", \"1\"]\),...,\!\(\*SubscriptBox[\"Q\", \"N\"]\),\!\(\*SubscriptBox[\"P\", \"1\"]\),...,\!\(\*SubscriptBox[\"P\", \"N\"]\))."

gsBeamSplitter::usage=
" gsBeamSplitter[N,m1,m2,t] returns the symplectic transformation on N modes
\n corresponding to a beam splitter operation of transmittivity t acting on modes m1 and m2."

gsPhaseShifter::usage=
" gsPhaseShifter[N,m,\[Phi]] returns the symplectic transformation on N modes
\n corresponding to a phase shifter operation by the angle \[Phi] acting on mode m."

gsSqueezer::usage=
"gsSqueezer[N,{\!\(\*SubscriptBox[\"r\", \"1\"]\),...,\!\(\*SubscriptBox[\"r\", \"N\"]\)}] returns the squeezing transformation on N moddes with squeezing parameter {\!\(\*SubscriptBox[\"r\", \"1\"]\),...,\!\(\*SubscriptBox[\"r\", \"N\"]\)}.
The vector of the squeezing parameters can be omitted; its default value is {\!\(\*SubscriptBox[\"r\", \"1\"]\),...,\!\(\*SubscriptBox[\"r\", \"N\"]\)}."

gsPTTrafo::usage=
" gsPTTrafo[N,{m_1,..,m_k}] returns the symplectic trafo on N modes corresponding to a partial transposition
\n with respect to a bipartition, where one party is comprised of the modes m_1 to m_k. "

gsSymplecticSpectrum::usage=
" gsSymplecticSpectrum[M] returns a list of the symplectic eigenvalues of a covariance matrix M.
\n According to Williamson's theorem every positive and symmetric matrix (hence every covariance matrix and its partial transpose)
\n admits a Williamson normal-mode decomposition of the form M = S^T \[Nu] S, where S is symplectic, \[Nu] is a diagonal matrix of the form  diag(\[Nu]_1,\[Nu]_1,..,\[Nu]_N,\[Nu]_N)
\n and the N (N=number of modes) numbers \[Nu]_i are called the symplectic eigenvalues of M.
\n\n The symplectic spectrum of M can be computed as the eigenvalues of the matrix |\[ImaginaryI] \[CapitalOmega] M|.
\n\n For any positive matrix M the Heisenberg condition for \"physical states\" M + \[ImaginaryI] \[CapitalOmega] \[GreaterEqual] 0 (\[CapitalOmega] is the Symplectic Form) 
\n is completely equivalent to \[Nu]_i\[GreaterEqual] (2 \[Kappa]2\!\(\*SuperscriptBox[\")\", 
RowBox[{\"-\", \"2\"}]]\) (canonical convention: \[GreaterEqual] 1/2; quantum-optical convention: \[GreaterEqual] 1/4) for all i\[Element]{1,..,N}.
\n The purity condition reads \[Nu]_i = (2 \[Kappa]2\!\(\*SuperscriptBox[\")\", 
RowBox[{\"-\", \"2\"}]]\) (canonical convention: = 1/2; quantum-optical convention: = 1/4)  for all i\[Element]{1,..,N}."

gsPTSymplecticSpectrum::usage=
" gsPTSymplecticSpectrum[M,{m_1,..,m_k}] returns the symplectic spectrum of the partial transpose of a covariance matrix M
\n with respect to a bipartition, where one party is comprised of the modes m_1 to m_k. 
\n If the partially transposed CM still represents a physical state (= if M represents a separable state wrt to this bipartition)
\n then all Symplectic eigenvalues are greater than or equal to 1.
\n Note that - different from the symplectic spectrum of a physical covariance matrix - the elements PTSymplecticSpectrum do
\n in general not appear twice! "

gsHeisenbergCheck::usage=
" gsHeisenbergCheck[M] returns true if and only if the covariance matrix M fulfills the 
\n conditions a physical state has to fulfill due to the Heisenberg uncertainty relations."

gsPPTSeparabilityCheck::usage=
" gsPPTSeparabilityCheck[M,{m_1,...,m_k}] returns true if the covariance matrix M fulfills the 
\n necessary conditions for being separable with respect to a bipartition, where one party is comprised of the modes m_1 to m_k.
\n If it returns false (or a relation which cannot be fulfilled!), the state is therefore entangled with respect to this bipartition.
\n However \"false\" is only sufficient AND NECESSARY if we are looking at 1xN-mode Gaussian States. "

gsSymplecticityCheck::usage=
"gsSymplecticityCheck[M] returns true, if the matrix M is symplectic."

gsPurity::usage=
"gsPurity[M] returns the purity \[Mu] \[Congruent] tr[ \!\(\*SuperscriptBox[\"\[Rho]\", \"2\"]\) ] of the Gaussian state \[Rho], which is characterized by the covariance matrix M."
gsPurityCheck::usage=
"gsPurityCheck[M] returns True if the CM M corresponds to a pure state."

gsNegativity::usage=
" gsNegativity[M,{m_1,..,m_k}] returns the negativity of the covariance matrix M with respect to the bipartition, 
\n where one party is comprised of the modes m_1 to m_k. It only works with covariance matrices of numbers."

gsLogNegativity::usage=
" gsLogNegativity[M,{m_1,..,m_k}] returns the logarithmic negativity of the covariance matrix M with respect to the bipartition, 
\n where one party is comprised of the modes m_1 to m_k. \n
The logarithmic negativity can be calculated as follows: \"Renormalize\" the symplectic spectrum of the partially transposed CM \n
by multiplying it by (2Private`\[Kappa]2)^2. Take the -Log of all those \"renormalized\" sympl. eigenvalues, which are less than 1. Sum everything up.\n
If all \"renormalized\" sympl. eigenvalues are greater than or equal to 1, the state is not entangled and the log. Negativity is zero. \n\n
 gsLogNegativity only works with covariance matrices of numbers (numerical input)."

gsEoF::usage=
"gsEoF[{a,b,cp,cm}] returns the Entanglement of Formation of an arbitrary two-mode Gaussian state with standard form parameters a, b, cp and cm.\n
gsEoF determines the value numerically and works only for numerical input."
gsEoFverbose::usage=
"gsEoFverbose[{a,b,cp,cm}] is the verbose version of the gsEoF command. Unlike the gsEoF command it prints warnings if the state is separable 
or does not fulfill the necessary criteria to be processed by gsEoF. It is not to be used within Plot commands (due to the potentially very high number of warnings.
It returns the Entanglement of Formation of an arbitrary two-mode Gaussian state with standard form parameters a, b, cp and cm.\n
gsEoFverbose determines the value numerically and works only for numerical input."

gsLinearClusterAdjacencyMat::usage=
"gsLinearClusterAdjacencyMat[n] returns the adjacency matrix of a linear graph on n nodes; corresponding to a n-mode linear cluster state."

gsSymLinearClusterAdjacencyMat::usage=
"gsSymLinearClusterAdjacencyMat[n] returns the adjacency matrix of a \"ring\"-graph on n nodes."

gsSquareClusterAdjacencyMat::usage=
"gsSquareClusterAdjacencyMat[N] returns the adjacency matrix of a square graph on n x n nodes; 
\n corresponding to a nxn-mode two-dimensional cluster state."
gsSymSquareClusterAdjacencyMat::usage=
"gsSymSquareClusterAdjacencyMat[n] returns the adjacency matrix of a \"symmetrized\" square graph on n x n nodes; 
\n corresponding to a nxn-mode \"symmetrized\" two-dimensional cluster state."

gsCubicClusterAdjacencyMat::usage=
"gsCubicClusterAdjacencyMat[n] returns the adjacency matrix of a cubic graph on n x n x n nodes; 
\n corresponding to a nxnxn-mode three-dimensional cluster state."

gsGraphState::usage=
"gsGraphState[M,r] returns the covariance matrix of the canonical graph state corresponding to the graph specified by the adjacency matrix M
\n starting from vacuum modes finitely squeezed by an amount r (optional argument, default value r)."

gsLinearCluster::usage=
" gsLinearCluster[n, r] returns the covariance matrix of a canonical n-mode linear cluster state with initial squeezing r (optional, default r)." 
gsSquareCluster::usage=
" gsSquareCluster[n, r] returns the covariance matrix of a canonical n x n-mode square cluster state with initial squeezing r (optional, default r)." 
gsCubicCluster::usage=
" gsCubicCluster[n, r] returns the covariance matrix of a canonical n x n x n-mode cubic cluster state with initial squeezing r (optional, default r)." 

gsPhotonNumber::usage=
"gsPhotonNumber[M, n] returns the expectation value of the photon number in mode n in the Gaussian state specified by covariance matrix M."
gsMeanPhotonNumber::usage=
"gsMeanPhotonNumber[M, n] returns the average of the expectation value of the photon number per mode in the Gaussian state specified by covariance matrix M."

(* just for personal convenience: undocumented auxiliary commands *)
Sigma1
Sigma2
Sigma3
id
epspath
jpgpath

gsConventionInfo

(*  END USAGE DOCUMENTATION  *)
(*  ================================================================= *)

(*  ================================================================= *)
(* START DEFINITIONS *)

(* Begins the private environment *) 
Begin["Private`"]
(* convention for the definition of the mode operators. Default is the canonical version with [q,p]=I *)
\[Kappa]1=\[Kappa]2=\[Kappa]3=1/Sqrt[2]; (* Default convention *)
(* \[Kappa]1=\[Kappa]2=1;\[Kappa]3==1/2; Quantum-optical convention. *)
gsConventionInfo:=If[\[Kappa]1===\[Kappa]2===\[Kappa]3===1/Sqrt[2],Print["Canonical Convention:  "],Print["Unknown convention is used."]];

(* numerical constants *)
\[Epsilon]=10^(-5); (* Numerical precision to use for functions like Chop *)

(* auxiliary definitions *)
id=IdentityMatrix[2]; 
Sigma1={{0,1},{1,0}};Sigma2={{0,-I},{I,0}};Sigma3={{1,0},{0,-1}};

(* path for importing graphics in JPEG-format or and exporting graphics in EPS-format*)
epspath[filename_Symbol]:="/home/jason/Desktop/diploma_thesis/graphics/"<>ToString[filename]<>".eps"
jpgpath[filename_Symbol]:="/home/jason/Desktop/diploma_thesis/mathematica/Graphics/"<>ToString[filename]<>".jpg"
epspath[filename_String]:="/home/jason/Desktop/diploma_thesis/graphics/"<>filename<>".eps"
jpgpath[filename_String]:="/home/jason/Desktop/diploma_thesis/mathematica/Graphics/"<>filename<>".jpg"

(* function definitions *)
gsDirectSum[a_?MatrixQ, b_?MatrixQ] := Module[{lengtha, lengthb, j, k, l, m, ttt, dim},
                                                               lengtha = Length[a]; lengthb = Length[b];
                                                                dim = lengtha + lengthb; ttt = Table[0, {j, dim}, {k, dim}];
                                                               Do[ttt[[j, k]] = a[[j, k]],
                               {j, 1, lengtha}, {k, 1, lengtha}];
    Do[ttt[[l + lengtha, m + lengtha]] = b[[l, m]],
                               {l, 1, lengthb}, {m, 1, lengthb}];
                                                         ttt ]

gsVacuumState[n_Integer]:=(1/(4 \[Kappa]1^2))IdentityMatrix[2 n];
gsVac[n_Integer]:=gsVacuumState[n];

gsSqueezedState[n_?IntegerQ, v_?VectorQ]:=(1/(4 \[Kappa]1^2))DiagonalMatrix[
												Flatten[Table[{Exp[2 v[[i]]],Exp[-2 v[[i]]]},{i,1,n}]]
														]
gsTMSS[r_:Global`r]:=1/2 ArrayFlatten[{{Cosh[2 r] id ,Sinh[2r] Sigma3},{Sinh[2 r] Sigma3,Cosh[2 r] id}}];
(*  gsTMSS=gsTMSS[]; does not work as expected*)

gsStandardFormTMGS[a_:Global`a,b_:Global`b,cp_:Global`cp,cm_:Global`cm]:=ArrayFlatten[{{a id,DiagonalMatrix[{cp,cm}]},{DiagonalMatrix[{cp,cm}],b id}}];
(*  gsStandardFormTMGS=gsStandardFormTMGS[]; does not work as expected*)

gsSFParams[cm_?MatrixQ]:=Module[{x1,x2,x3,x4,\[Alpha]},
								x1=Det[Take[cm,{1,2},{1,2}]];x2=Det[Take[cm,{3,4},{3,4}]];
									x3=Det[Take[cm,{1,2},{3,4}]];x4=Det[cm];\[Alpha]=Sqrt[((Sqrt[x1 x2] + x3)^2 - x4)/Sqrt[x1 x2]] ;
										{Sqrt[x1],Sqrt[x2],1/2(\[Alpha] + Sqrt[\[Alpha]^2 -4 x3]),1/2(\[Alpha] - Sqrt[\[Alpha]^2 -4 x3])}] 
(*returns the parameters a,b,Subscript[c, +],Subscript[c, -] of the standard form of a two-mode CM; in terms of the local symplectic invariants as given in Giedke (PhD-thesis)  *)
gsStandardFormTMGSParameters[cm_?MatrixQ]:=gsSFParams[cm];

gsSymplecticForm[n_?IntegerQ]:=Module[{\[Omega],aux\[Omega]},\[Omega]={{0,1},{-1,0}};aux\[Omega]=\[Omega];
														Do[aux\[Omega]=gsDirectSum[\[Omega],aux\[Omega]],{n-1}];
																aux\[Omega]
									]

gsWignerFunction[cm_?MatrixQ]:=Module[{N,inv,vec},N=1/2 Length[cm];inv=Simplify[Inverse[cm]];vec=gsSymplecticVector[N];
									Simplify[1/(((2\[Pi])^N)(\[Kappa]2^2N)Sqrt[Det[cm]])  Exp[-1/2 vec.inv.vec],Element[#,Reals]&/@vec]
									]

(* SymplecticVector[n_?IntegerQ]:=Module[{v},v=ArrayFlatten[
											Table[
											{ToExpression["q"<>ToString[i]],ToExpression["p"<>ToString[i]]}
													,{i,1,n}] ,1];
										v=Table[Operator[v[[i]]],{i,1,2n}];
										Do[ QASet[ Commutator[ v[[i-1]],v[[i]] ]=I*Operator[1] ],{i,2,2n,2}];
									v] *)
(* The above definition of SymplecticVector as a vector of operators based on QuantumAlgebra doesn't work yet.
Problems: Transpose[v] should return just a transposed vector of the SAME operators ! (   maybe it works ;)  )*)

gsSymplecticVector[n_?IntegerQ]:=ArrayFlatten[
											Table[
											{ToExpression["q"<>ToString[i]],ToExpression["p"<>ToString[i]]}
													,{i,1,n}] ,1]
									(* <> is the StringJoin operation *)

gsSymplecticTrafo[n_?IntegerQ,hvec_?ListQ]:=Module[{hmataux},
												hmataux=SparseArray[Table[{hvec[[i]][[1]],hvec[[i]][[2]]}->hvec[[i]][[3]],{i,1,Length[hvec]}],2n];
												MatrixExp[-(1/(2*\[Kappa]1^2))gsSymplecticForm[n].(hmataux+Transpose[hmataux])]
												] (* !! Change compared to Bru\[SZ]/Leuchs Exp(-\[CapitalOmega] H) instead of Exp(H \[CapitalOmega])  *)

gsModeInterchange[mat_?MatrixQ,m1_Integer,m2_Integer]:=Module[{auxmat,id,l},
																id=IdentityMatrix[2];l=1/2 *Length[mat];
																auxmat=(ReplacePart[IdentityMatrix[l],
Join[Table[{i,i}->id,{i,Delete[Table[i,{i,1,l}],{{m1},{m2}}]}],{{m1,m1}->0,{m2,m2}->0,{m1,m2}->id,{m2,m1}->id}]])//ArrayFlatten;
																auxmat.mat.auxmat]
gsQPOrdering[cm_?MatrixQ/;And@@EvenQ[Dimensions[cm]]]:=Module[{l,permlist},l=1/2 Length[cm];
																	permlist=Riffle[Table[i,{i,1,l}],Table[i,{i,l+1,2 l}]];
																IdentityMatrix[2 l][[permlist]].cm.Transpose[IdentityMatrix[2 l][[permlist]]]]

gsBeamSplitter[n_?IntegerQ,m1_?IntegerQ,m2_?IntegerQ,t_]:=Module[{id,zero,list},id=IdentityMatrix[2];zero=0*id;
																		list=Table[Table[zero,{n}],{n}];
																		Do[list[[i,i]]=id,{i,1,n}];
																		list[[m1,m1]]=list[[m2,m2]]=Sqrt[t]*id;
																		list[[m1,m2]]=Sqrt[1-t]*id;
																		list[[m2,m1]]=-Sqrt[1-t]*id;
															ArrayFlatten[list]
																]

gsPhaseShifter[n_?IntegerQ,m_?IntegerQ,\[Phi]_]:=Module[{id,start},id=IdentityMatrix[2];start=IdentityMatrix[2*n];
									start[[2m-1,2m-1]]=Cos[\[Phi]];start[[2m-1,2m]]=Sin[\[Phi]];start[[2m,2m-1]]=-Sin[\[Phi]];start[[2m,2m]]=Cos[\[Phi]];
																start
												]

gsSqueezer[n_?IntegerQ,vec_:Automatic]:=If[vec===Automatic,
														DiagonalMatrix[Exp[Riffle[Array[Subscript[r,#]&,n],-Array[Subscript[r,#]&,n]]]],
														DiagonalMatrix[Exp[Riffle[vec,-vec]]]]/;(VectorQ[vec]\[Or]vec===Automatic)

gsPTTrafo[n_?IntegerQ,bipartitionvec_?ListQ]:=Module[{auxtrafo},auxtrafo=IdentityMatrix[2n];
						Do[auxtrafo[[2bipartitionvec[[i]],2bipartitionvec[[i]]]]=-1,{i,1,Length[bipartitionvec]}];
													auxtrafo]

(* ======== Spectra, Checks and Entanglement Measures ========================= *)

gsSymplecticSpectrum[cm_?MatrixQ]:=Module[{dim,auxmat,eval},
										dim=(1/2)Length[cm];auxmat=I gsSymplecticForm[dim].cm;eval=Eigenvalues[auxmat];
										Abs/@eval//Simplify(* //Union problematic, if some evs appear e.g. four times *)
										]

gsPTSymplecticSpectrum[cm_?MatrixQ,bipartvec_?VectorQ]:=Module[{cmpt,dim},dim=(1/2)Length[cm];
															cmpt=gsPTTrafo[dim,bipartvec].cm.gsPTTrafo[dim,bipartvec];
															gsSymplecticSpectrum[cmpt]]

(*gsHeisenbergCheck[cm_?MatrixQ]:=Module[{dim,eval,l},dim=1/2*Length[cm];eval=Abs@Eigenvalues[cm+(I/(4 \[Kappa]1^2))*gsSymplecticForm[dim]];l=Length[eval];
                                     Simplify[Apply[And,Table[eval[[i]]>=(2 \[Kappa]2)^(-2),{i,1,l}]]]
									]*)
gsHeisenbergCheck[cm_?MatrixQ]:=And[(Chop[cm-Transpose[cm],\[Epsilon]]==0*IdentityMatrix[Length[cm]]//Simplify),(* symmetric *)
									And@@((#>=0&)/@Eigenvalues[cm]),(* positive *)
										And@@((#>=(2 \[Kappa]2)^(-2)&)/@gsSymplecticSpectrum[cm])](* and congruent to a valid thermal state *)
                                     

gsPurity[cm_?MatrixQ]:=Simplify[1/((2 \[Kappa]2)^(Length[cm])Sqrt[Det[cm]])]
gsPurityCheck[cm_?MatrixQ]:=Module[{dim},dim=1/2*Length[cm];Chop[Det[cm]-(2 \[Kappa]2)^(-4 dim), \[Epsilon]]==0
									]
gsPPTSeparabilityCheck[cm_?MatrixQ,bipartvec_?VectorQ]:=Module[{cmpt,dim,eval},dim=1/2*Length[cm];cmpt=gsPTTrafo[dim,bipartvec].cm.gsPTTrafo[dim,bipartvec];
													gsHeisenbergCheck[cmpt]
															]

gsSymplecticityCheck[mat_?MatrixQ]:=And@@EvenQ/@Dimensions[mat]\[And]Equal@@Dimensions[mat]\[And]Chop[mat.gsSymplecticForm[1/2 Length[mat]].Transpose[mat]-gsSymplecticForm[1/2  Length[mat]],\[Epsilon]]==0.0*IdentityMatrix[Length[mat]]//Simplify

gsNegativity[cm_,bipartvec_?VectorQ]:=Module[{list},list=Select[gsPTSymplecticSpectrum[cm,bipartvec],Simplify[#<1]&];
											(1/2)*( Power[Times@@list,-1/2]-1) (* our symplectic spectrum contains all eigenvalues of |i\[CapitalOmega]\[Sigma]|, i.e. every ev appears at least twice *)
											]/;MatrixQ[cm,NumericQ]
(* Negativity[cm_?MatrixQ,bipartvec_?VectorQ]:=
(1/2)*( PTSymplecticSpectrum[cm,bipartvec][[1]]^(-1)-1)/;Length[bipartvec]==1 for 1xN bipartitions only the smallest eigenvalue matters*)

gsLogNegativity[cm_,bipartvec_?VectorQ]:=Module[{list},list=Select[gsPTSymplecticSpectrum[cm,bipartvec],Simplify[#<1]&];
											(-1/2)*(Plus@@(Log/@(((2Private`\[Kappa]2)^2)*list)))
(* our symplectic spectrum contains all eigenvalues of |i\[CapitalOmega]\[Sigma]|, i.e. every ev appears at least twice *)
											]/;MatrixQ[cm,NumericQ]

EoF[xm_]:=(xm+1/2)Log[xm+1/2]-(xm-1/2)Log[xm-1/2];
\[ScriptCapitalZ][{a_,b_,cp_,cm_}]:=a^2+b^2+2cp cm;
detSF[{a_,b_,cp_,cm_}]:=a^2 b^2-a b cm^2-a b cp^2+cm^2 cp^2;
\[ScriptCapitalD][{a_,b_,cp_,cm_}]:=a^2 (-(1/4)+b^2)-a b (cm^2+cp^2)+1/16 (-4 b^2+(1-4 cm cp)^2);
\[ScriptCapitalD]tilda [{a_,b_,cp_,cm_}]:=detSF[{a,b,cp,cm}] - 1/4(a^2 + b^2 +2 cp Abs[cm]) + 1/16;
(*coefficients of the quartic polynomial in p to obtain the value pm>=1*)
a0[{a_,b_,cp_,cm_}]:=(a b -cm^2)(a (a b- cm^2)-b/4)(b (a b -cm^2) - a/4);

a1[{a_,b_,cp_,cm_}]:=-(cp (a b -cm^2) + Abs[cm]/4)((a-b)^2(cp (a b -cm^2)+ Abs[cm]/4)+2a b (cp - Abs[cm])(a b -cm^2 -1/4));

a2[{a_,b_,cp_,cm_}]:=((a cp - b Abs[cm])(a Abs[cm] - b cp) + cp Abs[cm]\[ScriptCapitalZ][{a,b,cp,cm}])(detSF[{a,b,cp,cm}]+1/16)-2(a^2 b^2 -cp^2 cm^2)\[ScriptCapitalD][{a,b,cp,cm}] - cp Abs[cm]detSF[{a,b,cp,cm}];

a3[{a_,b_,cp_,cm_}]:=(-(cp/4)-(a b-cp^2) Abs[cm]) (2 a b (-(1/4)+a b-cp^2) (-cp+Abs[cm])+(a-b)^2 (cp/4+(a b-cp^2) Abs[cm]));

a4[{a_,b_,cp_,cm_}]:=(a b-cp^2) (-(b/4)+a (a b-cp^2)) (-(a/4)+b (a b-cp^2));
pEquation[{a_,b_,cp_,cm_}]:=(a0[{a,b,cp,cm}] + a1[{a,b,cp,cm}] p + a2[{a,b,cp,cm}] p^2 +a3[{a,b,cp,cm}]p^3+a4[{a,b,cp,cm}]p^4==0);
(* coefficients of the quadratic trinomial in y,which yields the optimalSubscript[y, m](as the smallest root of this trinomial),if the coefficients are evaluated for the solutionsSubscript[p, m]>=1 of the quartic polynomialSubscript[\[Sum], n]^4 Subscript[a, n]p^n=0above *)
b0[{a_,b_,cp_,cm_},p_]:=-\[ScriptCapitalD]tilda[{a,b,cp,cm}] p;
b1[{a_,b_,cp_,cm_},p_]:=-2Sqrt[p]((Abs[cm](a b -cp^2) + cp/4)p + (cp (a b - cm^2) + Abs[cm]/4));
b2[{a_,b_,cp_,cm_},p_]:=(a b - cp^2)p^2 + \[ScriptCapitalZ][{a,b,cp,cm}]p+(a b -cm^2);
yEquation[{a_,b_,cp_,cm_},p_]:=(b0[{a,b,cp,cm},p] +b1[{a,b,cp,cm},p] y + b2[{a,b,cp,cm},p] y^2 == 0);

gsEoFverbose[{a_?NumericQ,b_?NumericQ,cp_?NumericQ,cm_?NumericQ}]:=Module[{aReal,bReal,cpReal,cmReal},
If[Max[Im/@{a,b,cp,cm}]>10^(-5),
	Print["(* Standard Form parameters seem to be complex. Real numbers expected. *)"];,
	aReal=Re[a];bReal=Re[b];cpReal=Re[cp];cmReal=Re[cm];If[(* check if state is separable *)
													gsPPTSeparabilityCheck[gsStandardFormTMGS[aReal,bReal,cpReal,cmReal],{1}]==True,
													Print["(* State appears to be separable *)"];0,
													(*default action: calculate EoF *)
														Module[{params,pm,ym,xm},params={aReal,bReal,cpReal,cmReal};
																pm=Max @Select[p/.Solve[pEquation[params],p],#\[Element]Reals&];
															If[pm<1,Print["No pm bigger than one.
																Is the CM inseparable? Is it a two-mode-squeezed thermal state?"];,
													{ym=Min@Select[y/.(Solve[yEquation[params,pm],y]),#\[Element]Reals&];
														xm=Sqrt[ym^2 + 1/4];}];EoF[xm]]]]]

gsEoF[{a_?NumericQ,b_?NumericQ,cp_?NumericQ,cm_?NumericQ}]:=Module[{aReal,bReal,cpReal,cmReal},
If[Max[Im/@{a,b,cp,cm}]>10^(-5),
	Print["(* Standard Form parameters seem to be complex. Real numbers expected. *)"];,
	aReal=Re[a];bReal=Re[b];cpReal=Re[cp];cmReal=Re[cm];If[(* check if state is separable *)
													gsPPTSeparabilityCheck[gsStandardFormTMGS[aReal,bReal,cpReal,cmReal],{1}]==True,
													(*Print["(* State appears to be separable *)"];*)0,
													(*default action: calculate EoF *)
														Module[{params,pm,ym,xm},params={aReal,bReal,cpReal,cmReal};
																pm=Max @Select[p/.Solve[pEquation[params],p],#\[Element]Reals&];
															If[pm<1,(*Print["No pm bigger than one.
																Is the CM inseparable? Is it a two-mode-squeezed thermal state?"];*)0,
													{ym=Min@Select[y/.(Solve[yEquation[params,pm],y]),#\[Element]Reals&];
														xm=Sqrt[ym^2 + 1/4];}];EoF[xm]]]]]

gsLinearClusterAdjacencyMat[n_Integer]:=SparseArray[Table[{i,i+1},{i,1,n-1}]->Table[1,{n-1}],n]
gsSymLinearClusterAdjacencyMat[n_Integer]:=SparseArray[Union[Table[{i,i+1},{i,1,n-1}],{{n,1}}]->Table[1,{n}],n]

gsSquareClusterAdjacencyMat[n_Integer]:=Module[{auxmat1,auxmat2},auxmat1=auxmat2=gsLinearClusterAdjacencyMat[n];
									Do[auxmat1=gsDirectSum[auxmat1,auxmat2],{i,1,n-1}];Do[Do[auxmat1[[m*n+i,(m+1)n+i]]=1,{i,1,n}],{m,0,n-2}];
									auxmat1]
gsSymSquareClusterAdjacencyMat[n_Integer]:=Module[{auxmat1,auxmat2},auxmat1=auxmat2=gsSymLinearClusterAdjacencyMat[n];
										Do[auxmat1=gsDirectSum[auxmat1,auxmat2],{i,1,n-1}];
											Do[Do[auxmat1[[m*n+i,(m+1)n+i]]=1,{i,1,n}],{m,0,n-2}];
												Do[auxmat1[[n (n-1)+i,i]]=1,{i,1,n}];auxmat1]

gsCubicClusterAdjacencyMat[n_Integer]:=Module[{auxmat1,auxmat2},auxmat1=auxmat2=gsSquareClusterAdjacencyMat[n];
									Do[auxmat1=gsDirectSum[auxmat1,auxmat2],{i,1,n-1}];
										Do[Do[auxmat1[[m*n^2+i,(m+1)n^2+i]]=1,{i,1,n^2}],{m,0,n-2}];auxmat1];

gsGraphState[adjmat_?MatrixQ,r_:Global`r]:=Module[{auxcm,l,list,enttrafo},l=Length[adjmat];auxcm=gsSqueezedState[l,Table[r,{l}]];
										list=Join[#,{1}]&/@((2*#-1)&/@Position[Normal@adjmat,1,2]);enttrafo=gsSymplecticTrafo[l,list];
										Simplify@(enttrafo.auxcm.Transpose[enttrafo])]

gsLinearCluster[n_Integer,r_:Global`r]:=gsGraphState[gsLinearClusterAdjacencyMat[n],r]
gsSquareCluster[n_Integer,r_:Global`r]:=gsGraphState[gsSquareClusterAdjacencyMat[n],r]
gsCubicCluster[n_Integer,r_:Global`r]:=gsGraphState[gsCubicClusterAdjacencyMat[n],r]

gsPhotonNumber[cm_?MatrixQ,mode_Integer]:=Simplify[(\[Kappa]1)^2( cm[[2*mode-1,2*mode-1]] + cm[[2*mode,2*mode]] -(2(\[Kappa]1^2)))]
gsMeanPhotonNumber[cm_?MatrixQ]:=Module[{l},l=1/2*Length[cm];										
									(Sum[Simplify[(\[Kappa]1)^2(cm[[2*mode-1,2*mode-1]]+cm[[2*mode,2*mode]] - 1)],{mode,1,l}])/l]

(* END DEFINITIONS *)
(*  =================================================================  *)

End[]
(* Ends the private environment *)

EndPackage[];




