
#include "../includes/bblEM.h"
#include "../includes/likelihoodComputation.h"
using namespace likelihoodComputation;
#include "../includes/computeUpAlg.h"
#include "../includes/computeDownAlg.h"
#include "../includes/computeCounts.h"
#include "../includes/treeIt.h"
#include "../includes/fromCountTableComponentToDistance.h"
//#include <ctime> // use only for debugging.

bblEM::bblEM(tree& et,
			 const sequenceContainer& sc,
			 const seqContainerTreeMap& sctm,
			 const stochasticProcess& sp,
			 const Vdouble * weights,
			 const int maxIterations,
			 const MDOUBLE epsilon,
			 const MDOUBLE tollForPairwiseDist,
			 unObservableData*  _unObservableData_p,
			 const MDOUBLE* likelihoodLast) :
_et(et),_sc(sc),_sctm(sctm),_sp(sp),_weights(weights),_unObservableData_p(_unObservableData_p) 
{
	//time_t ltime1;
	//time( &ltime1 );
	_treeLikelihood = compute_bblEM(maxIterations,epsilon,tollForPairwiseDist,likelihoodLast);
	//time_t ltime2;
	//time( &ltime2 );
	//int t = static_cast<long>(ltime2 - ltime1);
	//LOG(4,<<"Overall running time for BBL = "<<t<<" sec"<<endl);
}

MDOUBLE bblEM::compute_bblEM(
			const int maxIterations,
			const MDOUBLE epsilon,
			const MDOUBLE tollForPairwiseDist,
			const MDOUBLE* likelihoodLast){
	allocatePlace();
	MDOUBLE oldL=VERYSMALL;
	MDOUBLE currL = VERYSMALL;
	tree oldT = _et;
	for (int i=0; i < maxIterations; ++i) {
		//time_t ltime1;
		//time( &ltime1 );
		computeUp();
		currL = likelihoodComputation::getTreeLikelihoodFromUp2(_et,_sc,_sp,_cup,_posLike,_weights,_unObservableData_p);
		LOG(4,<<"--- Iter="<<i<<" logL="<<currL<<endl);
		if(oldL<=currL){	// make sure not to use tree with lower likelihood then last computed likelihood (before BBL-EM)
			if(likelihoodLast){	// likelihood from external model
				if(*likelihoodLast<=currL)
					oldT = _et;	 // L didn't go down, update old tree
				else{
					LOG(4,<<"Likelihood went down compared pre-BBL oldL="<<*likelihoodLast<<" newL="<<currL<<" Do not update old tree"<<endl);
					break;
				}
			}
			else{
				oldT = _et;	 // L didn't go down
				LOG(7,<<"LikelihoodLast was not sent to bblEM"<<endl);
			}
		}
		else
			LOG(4,<<"Likelihood went down oldL="<<oldL<<" newL="<<currL<<" Do not update tree"<<endl);

		if (currL < oldL + epsilon) { // need to break
			if (currL<oldL) {
				_et = oldT;
				if(_unObservableData_p)
					_unObservableData_p->setLforMissingData(_et,&_sp);
				return oldL; // keep the old tree, and old likelihood
			} else {
                //update the tree and likelihood and return
				return currL;
			}
		}
		bblEM_it(tollForPairwiseDist);
		oldL = currL;
		//time_t ltime2;
		//time( &ltime2 );
		//int t = static_cast<long>(ltime2 - ltime1);
		//LOG(6,<<"Time BBL iteration = "<<t<<" sec"<<endl);
	}

	// in the case were we reached max_iter, we have to recompute the likelihood of the new tree...
	computeUp();
	currL = likelihoodComputation::getTreeLikelihoodFromUp2(_et,_sc,_sp,_cup,_posLike,_weights,_unObservableData_p);
	if (currL<oldL) {
		_et = oldT;
		if(_unObservableData_p)
			_unObservableData_p->setLforMissingData(_et,&_sp);
		return oldL; // keep the old tree, and old likelihood
	} 
	else 
        return currL;
}


/********************************************************************************************
*********************************************************************************************/
void bblEM::allocatePlace() {
	_computeCountsV.resize(_et.getNodesNum()); //initiateTablesOfCounts
	for (size_t i=0; i < _computeCountsV.size(); ++i) {
		_computeCountsV[i].countTableComponentAllocatePlace(_sp.alphabetSize(),_sp.categories());
	}
	_cup.allocatePlace(_sc.seqLen(),_sp.categories(), _et.getNodesNum(), _sc.alphabetSize());
	_cdown.allocatePlace(_sp.categories(), _et.getNodesNum(), _sc.alphabetSize());
}

/********************************************************************************************
*********************************************************************************************/
void bblEM::bblEM_it(const MDOUBLE tollForPairwiseDist){
	//string costTable =  "countBBLEMTable.txt"; 
	//ofstream costTableStream(costTable.c_str());
	for (size_t i=0; i < _computeCountsV.size(); ++i) {
		_computeCountsV[i].zero();
		//_computeCountsV[i].printTable(costTableStream);
	}

	for (size_t i=0; i < _sc.seqLen(); ++i) {
		computeDown(i);
		addCounts(i); // computes the counts and adds to the table.
	}

	//for (int i=0; i < _computeCountsV.size(); ++i) { // used for Debug - check the need for 'zero()'
	//	_computeCountsV[i].printTable(costTableStream);
	//}

	optimizeBranches(tollForPairwiseDist);
	if(_unObservableData_p){
		_unObservableData_p->setLforMissingData(_et,&_sp);
	}
}
/********************************************************************************************
*********************************************************************************************/
void bblEM::optimizeBranches(const MDOUBLE tollForPairwiseDist){
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			fromCountTableComponentToDistance from1(_computeCountsV[mynode->id()],_sp,tollForPairwiseDist,mynode->dis2father(),_unObservableData_p);
			from1.computeDistance();
			mynode->setDisToFather(from1.getDistance());
			if(_unObservableData_p){	// needed only before likelihood computation
				_unObservableData_p->setLforMissingData(_et,&_sp);
			}
		}
	}
}
/********************************************************************************************
*********************************************************************************************/
void bblEM::computeUp(){
	_pij.fillPij(_et,_sp,0); // 0 is becaues we compute Pij(t) and not its derivations...
	computeUpAlg cupAlg;
	
	cupAlg.fillComputeUp(_et, _sc, _sctm, _pij, _cup);

	//for (size_t pos=0; pos < _sc.seqLen(); ++pos) {
        //for (size_t categor = 0; categor < _sp.categories(); ++categor) {
		//	cupAlg.fillComputeUp(_et,_sc,pos,_pij[categor], sctm,_cup[pos][categor]);
		//}
	//}
 }

void bblEM::computeDown(const int pos){
	computeDownAlg cdownAlg;
	for (size_t categor = 0; categor < _sp.categories(); ++categor) {
		cdownAlg.fillComputeDown(_et,_sc,pos,_pij[categor],_cdown[categor],_cup[pos][categor]);
	}
 }
/********************************************************************************************
*********************************************************************************************/
void bblEM::addCounts(const int pos){
	//MDOUBLE posProb = 
	//	likelihoodComputation::getProbOfPosWhenUpIsFilledGam(pos,_et,_sc,_sp,_cup);
						
	MDOUBLE weig = (_weights ? (*_weights)[pos] : 1.0);
	if (weig == 0) return;
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			addCounts(pos,mynode,_posLike[pos],weig);
		}
	}
}
/********************************************************************************************
*********************************************************************************************/
void bblEM::addCounts(const int pos, tree::nodeP mynode, const doubleRepMantisa posProb, const MDOUBLE weig){

	computeCounts cc;
	for (size_t categor =0; categor< _sp.categories(); ++ categor) {
			cc.computeCountsNodeFatherNodeSonHomPos(_sc,
										_pij[categor],
										_sp,
										_cup[pos][categor],
										_cdown[categor],
										weig,
										posProb,
										mynode,
										_computeCountsV[mynode->id()][categor],
										_sp.ratesProb(categor));
	}
}          

