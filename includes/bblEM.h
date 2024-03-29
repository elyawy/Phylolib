//Explanations by T.P on 4.6.2017
// BBL_EM stands for best branch lengths by the Expectation Maximization algorithm.
// By calling the construction function, all is fixed except the tree
// whose branch lengths are being optimized (hence it is not given as const).

#ifndef ___BBL_EM_H
#define ___BBL_EM_H

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "countTableComponent.h"
#include "computePijComponent.h"
#include "suffStatComponent.h"
#include "unObservableData.h"
#include "seqContainerTreeMap.h"
#include "doubleRep.h"
#include <vector>
using namespace std;

class bblEM {
public:
	explicit bblEM(tree& et,
		const sequenceContainer& sc,
		const seqContainerTreeMap& sctm,
		const stochasticProcess& sp,
		const Vdouble * weights = NULL,
		const int maxIterations=50,
		const MDOUBLE epsilon=0.05,
		const MDOUBLE tollForPairwiseDist=0.001,
		unObservableData*  unObservableData_p=NULL,
		const MDOUBLE* likelihoodLast=NULL);
	MDOUBLE getTreeLikelihood() const {return _treeLikelihood;}

private:
	MDOUBLE compute_bblEM(const int maxIterations,
					const MDOUBLE epsilon,
					const MDOUBLE tollForPairwiseDist,
					const MDOUBLE* likelihoodLast=NULL);
	void bblEM_it(const MDOUBLE tollForPairwiseDist);
	void computeDown(const int pos);
	void computeUp();
	void addCounts(const int pos);
	void addCounts(const int pos, tree::nodeP mynode, const doubleRepMantisa posProb, const MDOUBLE weig);
	void optimizeBranches(const MDOUBLE tollForPairwiseDist);
	void allocatePlace();


	MDOUBLE _treeLikelihood;
	tree& _et;
	const sequenceContainer& _sc;
	const stochasticProcess& _sp;
	const seqContainerTreeMap& _sctm;
	vector<countTableComponentGam> _computeCountsV; // for each node - a table of rate*alph*alph
	computePijGam _pij;
	suffStatGlobalGam _cup;
	suffStatGlobalGamPos _cdown;
	const Vdouble * _weights;
	vector<doubleRepMantisa> _posLike;
	unObservableData*  _unObservableData_p;
};

#endif
