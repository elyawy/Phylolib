#ifndef ___unObservableData___GL
#define ___unObservableData___GL

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "gainLossAlphabet.h"
#include "computePijComponent.h"

//*******************************************************************************************
//unObservableData
//*********************************************************************************************

//A theoretical note about un-observed data with among site rate variation:
//The question is whether one should correct for un-observed data for each rate category
//Or, for the entire position together. Short answer: for the entire position together.
//This is best explained by the following decision tree
//                  -----x
//                  |
//     rate1--------|-----observed data, not x
//       |          
//  -----|          -----column of zero
//       |          | 
//     rate2--------|-----observed data, not x
//                  |
//                  -----x
// We want to compute p(x | observed).
// According to the def of conditional probability, this is
//  p(x & observed) / p(observed) = p(x)/p(obs)
// According to the tree p(x) = p(rate1) * p(x/rate1) + p(rate2) * p(x/rate2)
// P(obs) = 1-p(column of zero) = 1-p(rate1)*p(column of zero|rate1)-p(rate2)*p(column of zero|rate2).
// The above it the RIGHT computation.
// The wrong computation would be sigma_r (p x|observed, r)*p(r).


class unObservableData{
public:
	explicit unObservableData() {};
	explicit unObservableData(const sequenceContainer& sc,
		const gainLossAlphabet& alph,
		const size_t minNumOfOnes,
		const size_t minNumOfZeros);

	unObservableData(const unObservableData& other);
	virtual ~unObservableData(){};
	virtual unObservableData* clone() const {return new unObservableData(*this);}
	MDOUBLE getlogLforMissingData() {return _logLforMissingData; };

	size_t getNumOfUnObservablePatterns() {return _scZero.seqLen();};
	void setLforMissingData(const tree& _tr, const stochasticProcess* _sp);
	void setLforMissingData(const tree& _tr, const vector<vector<stochasticProcess*> >& spVVec,	const distribution * distGain, const distribution* distLoss);

private:
//members
	sequenceContainer _scZero;
	MDOUBLE _logLforMissingData;
	computePijGam _pi;
};


#endif
