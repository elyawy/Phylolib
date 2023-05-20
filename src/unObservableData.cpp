
#include "unObservableData.h"
#include "likelihoodComputation.h"
#include "likelihoodComputationGL.h"
#include <cmath>


using namespace std;

unObservableData::unObservableData(	const sequenceContainer& sc,
									const gainLossAlphabet& alph,
									const size_t minNumOfOnes,
									const size_t minNumOfZeros) {
	_scZero.startZeroSequenceContainerGL(sc,alph, minNumOfOnes, minNumOfZeros);
}

unObservableData::unObservableData(const unObservableData& other) { //const

	_scZero = other._scZero;
	_pi = other._pi;
	_logLforMissingData = other._logLforMissingData;
}

void unObservableData::setLforMissingData(const tree& tr,
	const stochasticProcess* sp){
	seqContainerTreeMap sctm(_scZero, tr);
	_pi.fillPij(tr,*sp);
	_logLforMissingData = 0;
	MDOUBLE prob = 0;
	for (size_t pos = 0; pos < _scZero.seqLen(); ++pos) {
		prob = convert(likelihoodComputation::getLofPos(pos, tr, _scZero, sctm, _pi, *sp));
		if (prob <= 0 || prob > 1) {
			errorMsg::reportError("error in function void unObservableData::setLforMissingData", 1);
		}
		_logLforMissingData += log(prob);
	}
}
/********************************************************************************************
*********************************************************************************************/
void unObservableData::setLforMissingData(
		const tree& tr,
		const vector<vector<stochasticProcess*> >& spVVec,
		const distribution* distGain,
		const distribution* distLoss)
{
	seqContainerTreeMap sctm(_scZero, tr);
	_logLforMissingData = 0;
	int numOfRateCategories = spVVec[0][0]->categories();
	vector<computePijGam> pi_vec(numOfRateCategories);
	vector<suffStatGlobalGam> ssc_vec(numOfRateCategories);
	vector<computeUpAlg> cup_vec(numOfRateCategories);
	likelihoodComputationGL::fillPijAndUp(tr,_scZero, sctm,spVVec,distGain,distLoss,pi_vec,ssc_vec,cup_vec);
	
	for (size_t k=0; k < _scZero.seqLen(); ++k) {
		MDOUBLE resGivenRate = 0.0;
		MDOUBLE lnL = 0;
		for(int rateIndex=0 ; rateIndex<numOfRateCategories; ++rateIndex){
			lnL = log(likelihoodComputationGL::getProbOfPosUpIsFilledSelectionGam(k,//pos,
				tr,//const tree& 
				_scZero,// sequenceContainer& sc,
				spVVec,	// only needed for sp.freq(let)
				ssc_vec[rateIndex][k],//const computePijGam& ,
				distGain, distLoss)); // distributions
			resGivenRate += lnL * spVVec[0][0]->ratesProb(rateIndex);
		}
		_logLforMissingData += exp(resGivenRate);
	}
	_logLforMissingData = log(_logLforMissingData);
	//for(int rateIndex=0 ; rateIndex<numOfRateCategories; ++rateIndex){
	//	_logLforMissingData += likelihoodComputationGL::getTreeLikelihoodFromUp2(tr,_scZero,spVVec,ssc_vec[rateIndex], distGain,distLoss,NULL)
	//		* spVVec[0][0]->ratesProb(rateIndex);
	//}
}

