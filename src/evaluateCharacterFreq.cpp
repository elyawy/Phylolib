// $Id: evaluateCharacterFreq.cpp 10474 2012-03-18 07:54:07Z itaymay $

#include "../includes/evaluateCharacterFreq.h"
#include "../includes/someUtil.h"
#include <cassert>

vector<MDOUBLE> sumAlphabetCounts(const sequenceContainer & sc) {
	vector<MDOUBLE> charFreq(sc.alphabetSize(),0.0);
	sequenceContainer::constTaxaIterator tIt;
	sequenceContainer::constTaxaIterator tItEnd;
	tIt.begin(sc);
	tItEnd.end(sc);
	while (tIt!= tItEnd) {
		sequence::constIterator sIt;
		sequence::constIterator sItEnd;
		sIt.begin(*tIt);
		sItEnd.end(*tIt);
		while (sIt != sItEnd) {
			if ((*sIt >= 0) && (*sIt <charFreq.size())) ++charFreq[(*sIt)];
			++sIt;
		}
		++tIt;
	}
	return charFreq;
}

void changeCountsToFreqs(vector<MDOUBLE>& charFreq){
	MDOUBLE sumA = 0;
	size_t i=0;
	for (i=0; i < charFreq.size(); ++i) {
        sumA+=charFreq[i] ;
	}
	for (i=0; i < charFreq.size(); ++i) {
		charFreq[i] /= sumA;
	}
}

void makeSureNoZeroFreqs(vector<MDOUBLE> & charFreq){
	// CORRECT SO THAT THERE ARE NO ZERO FREQUENCIES.
	// ALL FREQS THAT WERE ZERO ARE CHANGED
	MDOUBLE ZERO_FREQ = 0.0000000001;
	MDOUBLE sumB=0;
	int charWithZeroFreq = 0;
	size_t i=0;
	for (i=0; i < charFreq.size(); ++i) {
		if (DSMALL_EQUAL(charFreq[i], ZERO_FREQ)) {
			charFreq[i] = ZERO_FREQ;
			++charWithZeroFreq;
		}
		else sumB +=charFreq[i];
	}
	if (!DEQUAL(sumB, 1.0, 0.01))
	{
		cerr.precision(10);	
		cerr<<"sumFreq = "<<sumB<<endl;
		//errorMsg::reportError("error in makeSureNoZeroFreqs(). Input frequencies must sum to 1.0");
	}
	scaleVec(charFreq, 1.0/charFreq.size());
	//MDOUBLE scaleFactor = sumB - (charWithZeroFreq * ZERO_FREQ);
	//for (i=0; i < charFreq.size(); ++i) {
	//	if (charFreq[i] != ZERO_FREQ)
	//		charFreq[i] *= scaleFactor;
	//}
}


vector<MDOUBLE> evaluateCharacterFreq(const sequenceContainer & sc) {
	vector<MDOUBLE> charFreq=sumAlphabetCounts(sc);
	changeCountsToFreqs(charFreq);
	makeSureNoZeroFreqs(charFreq);
	return charFreq;
}

VVdouble evaluateCharacterFreqOneForEachGene(const vector<sequenceContainer> & scVec){
	VVdouble charFreq;
	for (int k=0; k < scVec.size(); ++k) {
		charFreq.push_back(evaluateCharacterFreq(scVec[k]));
	}
	return charFreq;
}

		


vector<MDOUBLE> evaluateCharacterFreqBasedOnManyGenes(const vector<sequenceContainer> & scVec) {
	// note: all alphabets have to be the same!
	vector<MDOUBLE> charFreq(scVec[0].alphabetSize(),0.0);
	for (size_t i=0; i < scVec.size();++i) {
		assert(scVec[0].getAlphabet()->size()==scVec[i].getAlphabet()->size());
        vector<MDOUBLE> charFreqTmp=sumAlphabetCounts(scVec[i]);
		for (size_t z=0; z < charFreq.size();++z) charFreq[z]+=charFreqTmp[z];
	}
	changeCountsToFreqs(charFreq);
	makeSureNoZeroFreqs(charFreq);
	return charFreq;
}

//returns the number of each character in each position. 
//NOTE: returns also the number of unknown charecters in the last place in each vector, so that the actual vector size for each position is alphabetSize()+1
void getCharacterCounts(const sequenceContainer & sc, VVint& counts4pos)
{
	const alphabet* pAlph = sc.getAlphabet();
	int alphSize = sc.alphabetSize();
	size_t pos;
	counts4pos.resize(sc.seqLen());
	for (pos = 0; pos < sc.seqLen(); ++pos)
		counts4pos[pos].resize(alphSize + 1, 0);

	for (int seq = 0; seq < sc.numberOfSeqs();++seq) 
	{
		int id = sc.placeToId(seq);
		for (pos = 0; pos < sc.seqLen(); ++pos)
		{
			int charType = sc[id][pos];
			if (pAlph->isSpecific(charType))
			{
				++counts4pos[pos][charType];
			}
			else
				++counts4pos[pos][alphSize];
		}
	}
}

//returns the number of different character types in each position
void getCharacterType4pos(const sequenceContainer & sc, Vint& charactersType4pos)
{
	VVint counts4Pos;
	getCharacterCounts(sc, counts4Pos);
	charactersType4pos.resize(sc.seqLen(), 0);
	for (size_t pos = 0; pos < sc.seqLen(); ++pos)
	{
		for (int c = 0; c < counts4Pos[pos].size()-1; ++c)
		{
			if (counts4Pos[pos][c] > 0)
				++charactersType4pos[pos];
		}
	}
}

//returns the distribution of the different character types in each position along the whole alignment
void getCharacterTypeDistribution(const sequenceContainer & sc, Vint& charactersTypeDist)
{
	Vint charactersType4pos;
	getCharacterType4pos(sc, charactersType4pos);
	charactersTypeDist.resize(sc.numberOfSeqs()+1, 0);
	for (size_t pos = 0; pos < sc.seqLen(); ++pos)
	{
		int count = charactersType4pos[pos];
		++charactersTypeDist[count];
	}

}
