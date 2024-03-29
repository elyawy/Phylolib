// $Id: sequenceContainer.cpp 11751 2013-09-12 21:52:03Z cohenofi $

#include "../includes/sequenceContainer.h"
#include "../includes/logFile.h"
#include "../includes/someUtil.h"
#include "../includes/fastaFormat.h"

sequenceContainer::sequenceContainer(const sequenceContainer& other,const alphabet *inAlph) :
_generalRemarks(other._generalRemarks),
_id2place(other._id2place)
{
	for (size_t i=0; i < other._seqDataVec.size(); ++i)
		_seqDataVec.push_back(sequence(other._seqDataVec[i],inAlph));
}


//if bAugumentShorterSeqs=true then add gap characters at the end of short seqeunces
const int sequenceContainer::makeSureAllSeqAreSameLengthAndGetLen(bool bAugumentShorterSeqs) {
	if (_seqDataVec.size() == 0) return 0;
	const int len = _seqDataVec[0].seqLen();
	for (size_t i=1; i < _seqDataVec.size(); ++i) {
		if (_seqDataVec[i].seqLen()!=len) {
			if (bAugumentShorterSeqs) {
				for (int pos = _seqDataVec[i].seqLen(); pos < len; ++pos) 
					_seqDataVec[i].push_back(getAlphabet()->gap());
			}
			else {
                cerr<<_seqDataVec[i].name()<<" length = "<<_seqDataVec[i].seqLen()<<" "<<_seqDataVec[0].name()<<" length = "" "<<len<<endl;
                errorMsg::reportError("not all sequences are of the same lengths");
			}
		}
	}

	return len;
}

//void sequenceContainer::addFromsequenceContainer(sequenceContainer& seqToAdd){
//	if (_seqDataVec.empty()) { // first sequence to add
//		sequenceContainer::taxaIterator tit;
//		sequenceContainer::taxaIterator titEND;
//		tit.begin(seqToAdd);
//		titEND.end(seqToAdd);
//		while (tit!=titEND) {
//			_seqDataVec.push_back(*tit);
//
//		}
//	}
//	else {// now we are adding sequences to sequences that are already there.
//		sequenceContainer::taxaIterator tit;
//		sequenceContainer::taxaIterator titEND;
//		tit.begin(seqToAdd);
//		titEND.end(seqToAdd);
//		while (tit!=titEND) {
//			for (int i=0; i < _seqDataVec.size(); ++i) {
//				if (tit->name() == _seqDataVec[i].name()) {
//					_seqDataVec[i]+=(*tit);
//					break;
//				}
//			}
//			++tit;
//		}
//	}
//}

void sequenceContainer::changeGaps2MissingData() {

	for (size_t i = 0; i < seqLen();++i) {//going over al positions
		for (size_t j = 0; j < _seqDataVec.size();++j) {
			if (_seqDataVec[j][i] == getAlphabet()->gap()){
				 _seqDataVec[j][i]=getAlphabet()->unknown(); // missing data
			}
		}
	}
}

const int sequenceContainer::getId(const string &seqName, bool issueWarningIfNotFound) const {
	size_t k;
	for (k=0 ; k < _seqDataVec.size() ; ++k) {
		if (_seqDataVec[k].name() == seqName) return (_seqDataVec[k].id());
	}
	if (k == _seqDataVec.size() && issueWarningIfNotFound) {
		// debuggin
		LOG(5,<<"seqName = "<<seqName<<endl);
		for (k=0 ; k < _seqDataVec.size() ; ++k) {
			LOG(5,<<"_seqDataVec["<<k<<"].name() ="<<_seqDataVec[k].name()<<endl);
		}
		//end dubug
		LOG(0,<<seqName<<endl);
		vector<string> err;
		err.push_back("Could not find a sequence that matches the sequence name  ");
		err.push_back(seqName);
		err.push_back("in function sequenceContainer::getSeqPtr ");
		err.push_back(" make sure that names in tree file match name in sequence file ");
		errorMsg::reportError(err); // also quit the program
	}
	return -1;
}

const Vstring sequenceContainer::names() const {
	vector<string> res;
	for (size_t i=0; i < _seqDataVec.size(); ++i) {
		res.push_back(_seqDataVec[i].name());
	}
	return res;
}

sequenceContainer::sequenceContainer() {
	_id2place.resize(100,-1);
}

sequenceContainer::~sequenceContainer(){}

void sequenceContainer::add(const sequence& inSeq) {
	_seqDataVec.push_back(inSeq);
	if (_id2place.size() < inSeq.id()+1) {
		_id2place.resize(inSeq.id()+100,-1);
	}
	if (_id2place[inSeq.id()] != -1) {
		string err = "Two sequences with the same id - error in function sequenceContainer::add";
		err+= "\nThe id of the sequence you are trying to add = ";
		err += int2string(inSeq.id());
		errorMsg::reportError(err);
	}
	_id2place[inSeq.id()] = _seqDataVec.size()-1;
}


//given a sequence id the sequence is removed from the sequence container
//and the vector _id2place is updated.
void sequenceContainer::remove(const size_t idSeq)  {
	if (idSeq+1 > _id2place.size()) 	
		errorMsg::reportError("the id of sequence is not mapped by id2place in function sequenceContainer::remove");
	int place = _id2place[idSeq];
	
	if (place < 0) 
		errorMsg::reportError("cannot find place of the id in the sequence container in function sequenceContainer::remove");
	_seqDataVec.erase(_seqDataVec.begin()+place);

	_id2place[idSeq] = -1;
	for (size_t i=place;i<_seqDataVec.size();i++) {
		int id = _seqDataVec[i].id();
		_id2place[id]--;
	}
}
// remove all sequences from the sequence container
void sequenceContainer::removeAll(){
	Vint  ids2remove(numberOfSeqs());
	for(size_t i= 0; i<numberOfSeqs() ;i++){
		ids2remove[i] =placeToId(i);
	}
	for(size_t i= 0; i<ids2remove.size() ;i++){
		remove(ids2remove[i]);
	}
}


 
//removes identical sequences in the sequence container.
void sequenceContainer::removeIdenticalSequences(){
	bool exist;
	for (size_t i=1;i<_seqDataVec.size();i++){
		sequence sq1 = _seqDataVec[i];
		for (size_t j=0;j<i;j++){
			sequence sq2 = _seqDataVec[j];
			exist = true;
			if (sq1.seqLen() != sq2.seqLen()) continue;
			for (size_t pos=0;pos<sq1.seqLen();pos++){
				if (sq1[pos] != sq2[pos]){
					exist = false;
					break;
				}
			}
			if (exist) { 
				remove(sq1.id());
				i--;
				break;
				
			}

		}
	
	}

}

void sequenceContainer::removeGapPositions(){
	vector<int> posToRemove(seqLen(),0);
	bool gapCol;
	size_t i,j;
	for (i = 0; i < seqLen();++i) {//going over al positions
		gapCol = false;
		for (j = 0; j < _seqDataVec.size();++j) {
			if (_seqDataVec[j][i] == -1) posToRemove[i] = 1;
		}
	}
	removePositions(posToRemove);
}
void sequenceContainer::removeGapPositionsAllSeqs(){
	vector<int> posToRemove(seqLen(),1);
	bool gapCol;
	size_t i,j;
	for (i = 0; i < seqLen();++i) {//going over al positions
		gapCol = false;
		for (j = 0; j < _seqDataVec.size();++j) {
			if (_seqDataVec[j][i] != -1) posToRemove[i] = 0;
		}
	}
	removePositions(posToRemove);
}
void sequenceContainer::removeGapPositionsAccordingToAReferenceSeq(const string & seqName){
	int idOfRefSeq = getId(seqName,true);
	vector<int> posToRemove(seqLen(),0);
	size_t i;
	for (i = 0; i < seqLen();++i) {//going over al positions
		if (_seqDataVec[idOfRefSeq][i] == -1) posToRemove[i] = 1;
	}
	removePositions(posToRemove);
}

void sequenceContainer::removeUnknownPositionsAccordingToAReferenceSeq(const string & seqName){
	int idOfRefSeq = getId(seqName,true);
	vector<int> posToRemove(seqLen(),0);
	size_t i;
	for (i = 0; i < seqLen();++i) {//going over al positions
		if (_seqDataVec[idOfRefSeq][i] == getAlphabet()->unknown()) posToRemove[i] = 1;
	}
	removePositions(posToRemove);
}

//removePositions: the positions to be removed are marked as '1' in posToRemoveVec
//all othehr positions are '0' 	
void sequenceContainer::removePositions(const Vint & posToRemoveVec) {
	for (size_t z = 0; z < _seqDataVec.size();++z) {
		_seqDataVec[z].removePositions(posToRemoveVec);
	}
}


sequenceContainer sequenceContainer::getSubSeq(const int startPos, const int endPos) {
	sequenceContainer subSeq(*this);

	vector<int> posToRemove(seqLen(),true);
	for (int i = startPos; i <= endPos;++i) {//going over al positions
		posToRemove[i] = false;
	}
	subSeq.removePositions(posToRemove);

	return subSeq;
}


void sequenceContainer::changeDotsToGoodCharacters() {
	for (size_t i = 0; i < seqLen();++i) {//going over al positions
		int charInFirstSeq = _seqDataVec[0][i];
		if (charInFirstSeq == -3) {
			LOG(5,<<" position is "<<i<<endl);
			errorMsg::reportError(" the first line contains dots ");
		}
		for (size_t j = 1; j < _seqDataVec.size();++j) {
			if ((_seqDataVec[j][i] == -3)) {
				_seqDataVec[j][i] = charInFirstSeq; // missing data
			}
		}
	}
}

int sequenceContainer::numberOfSequencesWithoutGaps (const int pos) const {
	int numOfNonCharPos = numberOfSeqs();
	int gap = getAlphabet()->gap();
	int unknown = getAlphabet()->unknown();
	for (size_t i=0; i < numberOfSeqs(); ++i) {
		if ((*this)[i][pos] >= gap || (*this)[i][pos] >= unknown) --numOfNonCharPos;
	}
	return numOfNonCharPos;
}

int sequenceContainer::numberOfSequencesWithoutUnknowns (const int pos) const {
	int numOfNonCharPos = numberOfSeqs();
	int unknown = getAlphabet()->unknown();
	for (size_t i=0; i < numberOfSeqs(); ++i) {
		if ((*this)[i][pos] == unknown ) 
			--numOfNonCharPos;
	}
	return numOfNonCharPos;
}

bool sequenceContainer::isInvariable(const int pos) const {
	int charFound = getAlphabet()->unknown(); 
	for (size_t i=0; i < numberOfSeqs(); ++i) {
		if ((*this)[i][pos] >= 0) {
			if (charFound == getAlphabet()->unknown())
				charFound = (*this)[i][pos];
			else if (charFound != (*this)[i][pos])
				return false;
		}
	}
	return true;
}

int sequenceContainer::getInvariablePosNum() const {
	int sum = 0;
	for (size_t pos = 0; pos < seqLen(); ++pos) {
		if (isInvariable(pos))
			++sum;
	}
	return sum;
}

// new func for gainLoss project
void sequenceContainer::startZeroSequenceContainerGL(	const sequenceContainer &sc,
														const gainLossAlphabet& alph,
														const int minNumOfOnes,
														const int minNumOfZeros)
{
	//if(minNumOfOnes==0 && minNumOfZeros==0)
	//	return;

	string str0 = "0";
	string str1 = "1";
	vector<string> strV;
	strV.resize(sc.numberOfSeqs());
	string remark ="";
	switch (minNumOfOnes) {
			case (1) :
				for(size_t i=0; i<sc.numberOfSeqs();i++){
					// add patterns of 0 ones
					strV[i] = str0;
				}
				break;
			case (2) :
				for(size_t i=0; i<sc.numberOfSeqs();i++){
					// add patterns of 0 ones
					strV[i] = str0;
				}
				for(size_t i=0; i<sc.numberOfSeqs();i++){
					// add patterns of only 1 ones
					for(size_t j=0; j<sc.numberOfSeqs(); j++){
						if(j==i){
							strV[i]+=str1;
						}
						else{
							strV[i]+=str0;
						}
					}
				}
				break;
			case (3) :
				for(size_t i=0; i<sc.numberOfSeqs();i++){
					// add patterns of 0 ones
					strV[i] = str0;
				}
				for(size_t i=0; i<sc.numberOfSeqs();i++){
					// add patterns of only 1 ones
					for(size_t j=0; j<sc.numberOfSeqs(); j++){
						if(j==i){
							strV[i]+=str1;
						}
						else{
							strV[i]+=str0;
						}
					}
				}
				// add patterns of only 2 ones
				for(size_t onePosition1=0; onePosition1<sc.numberOfSeqs(); onePosition1++){
					for(size_t onePosition2=0; onePosition2<sc.numberOfSeqs(); onePosition2++){
						if(onePosition2<=onePosition1)
							continue;
						for(size_t i=0; i<sc.numberOfSeqs();i++){
							if(i==onePosition1 || i==onePosition2){
								strV[i]+=str1;
							}
							else{
								strV[i]+=str0;
							}			
						}
					}
				}				
				break;
	}
	switch (minNumOfZeros) {
			case (0) :
				break;			
			case (1) :
				for(size_t i=0; i<sc.numberOfSeqs();i++){
					// add patterns of 0 zeroes (only '1')
					strV[i] += str1;
				}
				break;
	}
	//////////////////////////////////////////////////////////////////////////
	
	for(size_t i=0; i<sc.numberOfSeqs();i++){
		//cout<<strV[i]<<endl;
		this->add(sequence(strV[i],sc.name(i),remark,i,&alph));
	}
}


//concatenate two sequecneContainers. 
//The sequence names must be identical in the two containers.
//returns false if: (1) A sequence_name in one of the containers does not match any sequence_name in the other container.
void sequenceContainer::concatenate(sequenceContainer& other) {
	if (other.numberOfSeqs() != numberOfSeqs()){
		string msg = "Not the same number of taxa, can't concatenate: other="+ int2string(other.numberOfSeqs()) + " this=" + int2string( numberOfSeqs()) +"\n";
		errorMsg::reportError(msg);
		return;
	}
	for (sequenceContainer::taxaIterator itThis=(*this).taxaBegin();itThis!=(*this).taxaEnd();++itThis) {
	//for(int i = 0; i < numberOfSeqs(); ++i)	{
		bool bFound = false;
		//out << (*this)[i].name()<<endl;

		for (sequenceContainer::taxaIterator itOther=other.taxaBegin();itOther!=other.taxaEnd();++itOther) {
        //for (int j = 0; j < other.numberOfSeqs(); ++j) {			
			//if((*this)[i].name().compare(other[j].name()) == 0)
			if(itThis->name().compare(itOther->name()) == 0)
			{
				//(*this)[i] += other[j]; // was i ?????
				*(itThis) += *(itOther);
				bFound = true;
				break;
			}
		}
		if (bFound == false) 
		{
			string msg = "Can't find sequence name in the second MSA: " +itThis->name();
            errorMsg::reportError(msg);			
		}
	}	
}
//////////////////////////////////////////////////////////////////////////
const bool sequenceContainer::operator==(const sequenceContainer& sq) const {
	if (_seqDataVec.size() != sq._seqDataVec.size())	// not the same number of sequences in  sequenceContainer
		return false;		
	const int numberOfSeqs = _seqDataVec.size();
	const int len = _seqDataVec[0].seqLen();
	for (int i=0; i < numberOfSeqs; ++i) {
		string nameI = name(i);
		int idI = getId(nameI);
		int idSq = sq.getId(nameI);
		if (_seqDataVec[idI].seqLen()!=sq._seqDataVec[idSq].seqLen())
			return false;
		for (int pos = 0; pos < len; ++pos)
		{
			if (_seqDataVec[idI][pos]!=sq._seqDataVec[idSq][pos])
				return false;
		}
	}
	return true;
}



//////////////////////////////////////////////////////////////////////////
int sequenceContainer::getNumOfOccurancesPerPos(const int pos, const char charId){
	int numOfOccurancesPerPos = 0;
	const int numberOfSeqs = _seqDataVec.size();
	const int len = _seqDataVec[0].seqLen();
	
	for (int i=0; i < numberOfSeqs; ++i) {
		string nameI = name(i);
		int idI = getId(nameI);
		if (_seqDataVec[idI][pos]==charId)
			numOfOccurancesPerPos++;
	}
	return numOfOccurancesPerPos;
}

//////////////////////////////////////////////////////////////////////////
vector<string> sequenceContainer::getSeqNamesThatMatchPos(const int pos, const char charId){
	vector<string> SeqNamesThatMatchPos;
	const int numberOfSeqs = _seqDataVec.size();
	const int len = _seqDataVec[0].seqLen();

	for (int i=0; i < numberOfSeqs; ++i) {
		string nameI = name(i);
		int idI = getId(nameI);
		if (_seqDataVec[idI][pos]==charId)
			SeqNamesThatMatchPos.push_back(nameI);
	}
	return SeqNamesThatMatchPos;
}

//////////////////////////////////////////////////////////////////////////
// added counts for unKnown data
const vector<size_t> sequenceContainer::getAlphabetDistribution() const {
	vector<size_t> alphabetVec;
	size_t alphSizePlusOne = alphabetSize()+1; //unKnown
	alphabetVec.resize(alphSizePlusOne);
	const size_t len = _seqDataVec[0].seqLen();
	const size_t numberOfSeqs = _seqDataVec.size();
	for (size_t pos = 0; pos < len; ++pos) {
		vector<size_t> alphabetVecPos;
		alphabetVecPos.resize(alphSizePlusOne);
		alphabetVecPos = getAlphabetDistribution(pos);
		for (size_t k=0; k < alphabetVec.size(); ++k) {
			alphabetVec[k] += alphabetVecPos[k];
		}
	}
	return alphabetVec;
}

//////////////////////////////////////////////////////////////////////////
const vector<size_t> sequenceContainer::getAlphabetDistribution(size_t pos) const {
	vector<size_t> alphabetVec;
	size_t alphSize = alphabetSize() + 1; //unKnown
	size_t UnknownVal = getAlphabet()->unknown();
	alphabetVec.resize(alphSize);
	const size_t numberOfSeqs = _seqDataVec.size();
	for (size_t i = 0; i < numberOfSeqs; ++i) {
		for (size_t alph = 0; alph < alphSize; ++alph) {
			if (_seqDataVec[i][pos] == alph) ++alphabetVec[alph];
			else if (_seqDataVec[i][pos] == UnknownVal)
				++alphabetVec[alph];
			
		}
	}
	return alphabetVec;
}



