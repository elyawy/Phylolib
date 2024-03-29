
#include "../includes/seqContainerTreeMap.h"
#include "../includes/errorMsg.h"
#include "../includes/treeUtil.h"
#include <cstdlib>
using namespace std;

// TAL PUPKO WENT OVER THIS FUNCTION ON 2.6.2017.
seqContainerTreeMap::seqContainerTreeMap(const sequenceContainer& sc, const tree& et) {
	fillSeqContainerTreeMap(sc, et);
}
void seqContainerTreeMap::fillSeqContainerTreeMap(const sequenceContainer& sc, const tree& et) {
	checkThatNamesInTreeAreSameAsNamesInSequenceContainer(et, sc);
	_V.resize(et.getNodesNum());
	treeIterTopDownConst tit(et);
	for (tree::nodeP myN = tit.first(); myN != tit.end(); myN = tit.next()) {
		if (myN->isInternal()) {
			_V[myN->id()] = -1;
		}
		else {
			_V[myN->id()] = sc.getId(myN->name(), false); // false -> do not issue a warning of not found. We know is it found because we checked it above.
		}
	}
}

// TAL PUPKO WENT OVER THIS FUNCTION ON 2.6.2017. Made it static.
//if bLeavesOnly == true then checks only leaves, otherwise the sequence container includes also internal nodes (as may be the result of simlations
void seqContainerTreeMap::checkThatNamesInTreeAreSameAsNamesInSequenceContainer(const tree& et, 
																				const sequenceContainer & sc,
																				bool bLeavesOnly) {
	treeIterDownTopConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		bool bFound = false;
		if (bLeavesOnly) {
			if (mynode->isInternal())
				continue;
		}
		sequenceContainer::constTaxaIterator it = sc.constTaxaBegin();
		for (; it != sc.constTaxaEnd(); ++it)
		{
			string scName = it->name();
			string treeNodeName = mynode->name();

			if (it->name() == mynode->name())
			{
				bFound = true;
				break;
			}
		}
		if (bFound == false)
		{
			string errMsg = "The sequence name: ";
			errMsg += mynode->name();
			errMsg += " was found in the tree file but not found in the sequence file.\n";
			errMsg += " Please, Re-run program with _intersectTreeAndSeq to produce new MSA and Tree.\n";
			errorMsg::reportError(errMsg);
		}
	}
}



void seqContainerTreeMap::intersectNamesInTreeAndSequenceContainer(tree& et, sequenceContainer & sc, bool bLeavesOnly){
	//LOG(4,<<"\n intersectNames Tree vs Sequence. Before intersect numOfSeq= "<<sc.numberOfSeqs()<<" nunOfTaxa= "<<et.getLeavesNum()<<" Remove "<<abs(et.getLeavesNum() -sc.numberOfSeqs())<<" taxa"<<endl);
	treeIterDownTopConst tIt(et);
	vector<tree::nodeP> nodes2remove;
	vector<int> seqIDs2remove;

	//cout<<"tree names:"<<endl;

	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		bool bFound = false;
		bool bFound_more = false;

		if (bLeavesOnly) {
			if (mynode->isInternal()) 
				continue;
		}
		sequenceContainer::constTaxaIterator it=sc.constTaxaBegin();
		for (;it != sc.constTaxaEnd(); ++it) 
		{
			string scName = it->name();
			string treeNodeName = mynode->name();

			if (it->name() == mynode->name()) 
			{
				if(bFound)
					bFound_more = true;
				bFound = true;
				//break;
			}
			if (bFound_more == true) 
			{
				string errMsg = "The taxID:\t";
				errMsg += mynode->name();
				errMsg += "\twas found again in the sequence file. Removed from sequence.";
				LOG(4,<<errMsg<<endl);
				seqIDs2remove.push_back(it->id());
				bFound_more = false;
			}
		}
		if (bFound == false) 
		{
			string errMsg = "The taxID:\t";
			errMsg += mynode->name();
			errMsg += "\twas found in the tree file but not found in the sequence file. Removed from tree.";
			LOG(4,<<errMsg<<endl);
			nodes2remove.push_back(mynode);			
		}

	}
	for(size_t i=0; i<nodes2remove.size(); ++i){
		et.removeLeaf(nodes2remove[i]);
	}
	sequenceContainer::constTaxaIterator myseq=sc.constTaxaBegin();
	for (;myseq != sc.constTaxaEnd(); ++myseq){
		bool bFound = false;
		bool bFound_more = false;
		for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {			
			if (bLeavesOnly)
			{
				if (mynode->isInternal()) 
					continue;
			}
			if (myseq->name() == mynode->name()) 
			{
				if(bFound)
					bFound_more = true;
				bFound = true;
				//break;
			}
			if (bFound_more == true) 
			{
				string errMsg = "The taxID name:\t";
				errMsg += myseq->name();
				errMsg += "\twas found again in the tree file. Removed.";
				LOG(4,<<errMsg<<endl);
				nodes2remove.push_back(mynode);
				bFound_more = false;
			}
		}
		if (bFound == false) 
		{
			string errMsg = "The taxID name:\t";
			errMsg += myseq->name();
			errMsg += "\twas found in the sequence file but not found in the tree file. Removed.";
			LOG(4,<<errMsg<<endl);
			seqIDs2remove.push_back(myseq->id());			
		}
	}
	for(size_t i=0; i<seqIDs2remove.size(); ++i){
		sc.remove(seqIDs2remove[i]);
	}
}



/********************************************************************************************
*********************************************************************************************/

/********************************************************************************************
// input: a tree and a sequence-container containing all of the leaves sequences.
// output: fills sc_leaves with the sequences of the leaves only.
*********************************************************************************************/
void seqContainerTreeMap::getLeavesSequences(const sequenceContainer& sc,
						const tree& tr, sequenceContainer& sc_leaves) {
	vector<string> leavesNames = getSequencesNames(tr);
	vector<string>::iterator itr_leaves;
	for (itr_leaves=leavesNames.begin();itr_leaves!=leavesNames.end();++itr_leaves) {
		sequenceContainer::constTaxaIterator it_sc=sc.constTaxaBegin();
		for (;it_sc != sc.constTaxaEnd(); ++it_sc) {
			if (it_sc->name() == *(itr_leaves)) {
				sc_leaves.add(*it_sc);
				break;
			}
		}
	}
	if (tr.getLeavesNum() != sc_leaves.numberOfSeqs()) {
		string errMsg = "getLeavesSequencese: the number of leaves is not equal to the number of leaves' sequences";
		errorMsg::reportError(errMsg);
	}
}
