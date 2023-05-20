// $Id: computeDownAlg.h 3107 2007-12-27 12:38:05Z adist $

#ifndef ___COMPUTE_DOWN_ALG
#define ___COMPUTE_DOWN_ALG

#include "definitions.h"
#include "tree.h"
#include "suffStatComponent.h"
#include "sequenceContainer.h"
#include "computePijComponent.h"

// The down algorithm is the partial likelihood of "a side of the tree".
// Assume the following phylogenetic tree
//					     ROOT
//					 /	       \
//					 N1			  N3
//					/ \         /    \
//				   sp1 N2     sp4   sp5
//					   / \        
//					sp2   sp3
//
// First the down in the root is set to 1 for each letter (we assume, for explaining, nucleotides).
//L(ROOT,A)=L(ROOT,C)=L(ROOT,G)=L(ROOT,T)=1.
// Next, we will compute the down component for the node N1.
// Specifically, we will compute L(N1,A). That is, letter = A;
// The "fatherNode" in this case is the ROOT.
// In this case, the "fatherTerm" = 1.0.
// The brother node is N3.
// We now compute the "brotherTerm", which is
// Sigma over  x = {A,C,G,T} {pij(A,x,t3) * UP[x,n3]}
// where t3, it the branch connecing the root to N3.
// Taking it together, we have 
// DOWN[N1,A] = totalProb = fatherTerm * brotherTerm = Sigma over  x = {A,C,G,T} {pij(A,x,t3) * UP[n3,x]};
// In a sense, what we compute is the partial likelihood of the subtree
//                      ----------------------
//					    | ROOT               |
//					 /	|       \            |
//					 N1	|		  N3         |
//					/ \  \       /    \      |
//				   sp1 N2  \   sp4   sp5     |
//					   / \   \               |
//					sp2   sp3 \--------------|
//
// So down[N1,A] is the likelihood of the subtree in the above figure when the root has A.
// We next go over the computation of DOWN[N2,A]
// This will be the likelihood of the following subtree
//              |----------------------------
//				|	     ROOT                |
//				|	 /	       \             |
//				|	 N1			  N3         |
//				|   /  /\------   /  \       |
//				|  sp1/  N2    \ sp4  sp5    |
//				-----/   / \    \            |
//					    sp2 sp3  \------------|
//
// Given that N1 is assigned A.
// First, let's compute the fatherTerm.
// FatherTerm = sigma over letterInFather { Pij(A,letterInFather,t1) * down[N1,letterInFather] }
// BrotherTerm = Sigma over  x = {A,C,G,T} {pij(A,x,t2) * UP[x,SP1]}
// Here, t1 is the branch above N1 and t2 is the branch above sp1.
// So the down terms are always the likelihood of a side tree that includes
// all the tree excluding the subtree below the node.
// The down for N2 is actually the partial likelihood of the subtree that does
// not include N2, given that the father node (N1) has a character "A".




// 
class computeDownAlg {
public: 
	//FOR A SPECIFIC POSITION AND NO RATE VARIATION.
	void fillComputeDown(const tree& et,
					   const sequenceContainer& sc,
					   const size_t pos,
					   const computePijHom& pi,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup);

	void fillComputeDown(const tree& et,
					   const sequenceContainer& sc,
					   const size_t pos,
					   const stochasticProcess& sp,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup);

	void fillComputeDownSpecificRate(const tree& et,
					   const sequenceContainer& sc,
					   const size_t pos,
					   const stochasticProcess& sp,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup,
					   const MDOUBLE gRate);

/** compute the down computation for a non-reversible model:
	each down computation is conditioned on the state at the root.
	This means that the vector field is of one additional dimension (the alphabet at the root)
	and hence the use of the suffStatGlobalGamPos (=vector<suffStatGlobalHomPos>)
**/
	void fillComputeDownNonReversible(const tree& et,
		const sequenceContainer& sc,
		const size_t pos,
		const computePijHom& pi,
		suffStatGlobalGamPos& sscGivenRoot,
		const suffStatGlobalHomPos& cup);
};
#endif
