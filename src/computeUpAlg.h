#ifndef ___COMPUTE_UP_ALG
#define ___COMPUTE_UP_ALG

#include "definitions.h"
#include "tree.h"
#include "suffStatComponent.h"
#include "sequenceContainer.h"
#include "computePijComponent.h"
#include "seqContainerTreeMap.h"

// The up algorithm computes the partial likelihood of each subtree of a given rooted tree.
// It is a dynamic-programming algorithm.
// First, let's describe the algorithm for a single position, with no rate variation (the simplest case).
// We will assume a nucleotide model, such as the HKY. The model actually gives us the Pij(t) entries for each branch.
// These Pij(t) values are stored in a class called computePijComponent.
// Assume the following tree
//			ROOT
//			/	\	
//			S1   S2
// Here we assume that in the position we analyze, S1 has a "C" and S2 has a "G".
// We have a vector of double of length 4 for each node.
// The beginning of the recursion is to assign the following vector to the node represented by S1:
// L(S1) = [0.0, 1.0, 0.0, 1.0]. We write L(S1, A) = 0 and L(S1, C) = 1.
// This indicates that the conditional likelihood of the subtree that includes only the leaf S1, given that the node has "C"
// is 1.0 and otherwise it is 0.
// The vector in S2 would be [0.0, 0.0, 0.0, 1.0]
// The vector in the root is next computed, given the vector of its childs.
// For example, the conditional likelihood of "A" in the root is
// P(A -> C | along the branch going to S1) * P(A->T | along the branch going to S2). A.M.  "P(A->G | along the branch going to S2)"??
// More generatlly, it is computed by the following formula
// [Sigma over x (x being A, C, G, T) P(A -> x | along the branch leading to S1) * L(S1,x)] *  [Sigma over y (y being A, C, G, T) P(A -> y | along the branch leading to S2) * L(S2,y)]
// While it seems that this formula is complicated, it covers both the case that the child nodes are leaves or internal nodes and a mixture of the two types of nodes.
// The final likelihood of the tree is
// Sigma over x L(ROOT,X) * p(x).

// We will now describe the most basic function void fillComputeUp(const tree& et, const sequenceContainer& sc,
//																   const int pos, const computePijHom& pi, suffStatGlobalHomPos& ssc);
// It must get the tree, so that it computes the partial likelihoods for each node.
// It must get the sequences ("sequenceContainer"), so that it can intiate the recursion.
// This function computes the partial likelihoods for a specific position, "pos".
// The Pij(t) values are given in the class computePijHom.
// The way this class is used is by calling pi.getPij(mynode->id(),letter,letInSon);
// This gives the Pij(t) from letter to letInSon, along the branch LEADING to mynode.
// The partial likelihood computed are stored in a structure called "suffStatGlobalHomPos".
// It is called this way because suffStat stands for sufficient statistics. It includes all the information
// about the likelihood in the subtree.
// Global is because it is for ALL nodes. For a single node it is suffStatComponent.
// Hom is because we do not assume among site rate variation (gamma). It assumes a rate of "1".
// Pos is beacause it only stores the partial likelihood for a specific position.
// It stores the partial likelihood by calling to the function
// ssc.set(mynode->id(),letter,total_prob);
// Internally, this structure is just a two dimentional array V = vector<vector<double> >
// where V[NODE][LETTER] is the partial likelihood S(NODE, LETTER) described above.


class computeUpAlg {
public: 
	// THIS FUNCTION IS FOR A SPECIFIC POSITION, NO RATE VARIATION, ALL NODES
	// IT IS THE MOST BASIC FUNCTION, AND IS USED BY THE OTHER FUNCTIONS
	void fillComputeUp(const tree& et,
					   const sequenceContainer& sc,
					   const size_t pos,
					   const computePijHom& pi,
					   const seqContainerTreeMap &sctm,
					   suffStatGlobalHomPos& ssc);

	// THIS FUNCTION FILLS THE UP COMPUTATION FOR ALL POSITION, AND ALL RATE CATEGORIES.
	void fillComputeUp(const tree& et,
					 const sequenceContainer & sc,
					 const seqContainerTreeMap & sctm,
					 const computePijGam& pi,
					 suffStatGlobalGam& ssc);

	void fillComputeUpSpecificGlobalRate(const tree& et,
				   const sequenceContainer& sc,
				   const size_t pos,
				   const stochasticProcess& sp,
				   suffStatGlobalHomPos& ssc,
				   const MDOUBLE gRate);

// my attemp to add factors
	void fillComputeUpWithFactors(const tree& et,
				   const sequenceContainer& sc,
				   const size_t pos,
				   const computePijHom& pi,
				   suffStatGlobalHomPos& ssc,
				   vector<MDOUBLE>& factors);

	void fillComputeUpWithFactors(const tree& et,
				   const sequenceContainer& sc,
				   const size_t pos,
				   const stochasticProcess& sp,
				   suffStatGlobalHomPos& ssc,
				   vector<MDOUBLE>& factors);

	void fillComputeUpSpecificGlobalRateFactors(const tree& et,
				   const sequenceContainer& sc,
				   const size_t pos,
				   const stochasticProcess& sp,
				   suffStatGlobalHomPos& ssc,
				   const MDOUBLE gRate,
				   vector<MDOUBLE>& factors);
};
#endif


