
#include "../includes/definitions.h"
#include "../includes/computeDownAlg.h"
#include "../includes/treeIt.h"

// THIS IS THE MOST BASIC FUNCTION. IT FILLS THE DOWN COMPONENT FOR A SPECIFIC POSITION
// AND ASSUMES NO RATE VARIATION (I.E., A SINGLE RATE CATEGORY).

void computeDownAlg::fillComputeDown(const tree& et,
					   const sequenceContainer& sc,
					   const size_t pos,
					   const computePijHom& pi,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup){
	ssc.allocatePlace(et.getNodesNum(), pi.alphabetSize()); // for each node and alphabet we will assign a value
	treeIterTopDownConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		size_t letter,letterInFather,bro,letterInSon;
		if (mynode->father()==NULL) {// if root
			for(letter=0; letter<pi.alphabetSize();letter++) {
				ssc.set(mynode->id(),letter,1.0);
			}
			mynode = tIt.next(); //continue
		}
		tree::nodeP fatherNode=mynode->father();
		const size_t n_bro=fatherNode->getNumberOfSons();
		for(letter=0; letter<pi.alphabetSize();letter++) {//alpha
			doubleRepMantisa totalProb=1.0;
			doubleRepMantisa fatherTerm=0;
			if (fatherNode->father()!=NULL) {
				for(letterInFather=0; letterInFather<pi.alphabetSize();letterInFather++)
					fatherTerm += pi.getPij(fatherNode->id(),letter,letterInFather)*
					ssc.get(fatherNode->id(),letterInFather);
			}
			else {
				fatherTerm=1.0;
			}
				doubleRepMantisa brotherTerm=1.0;
			for(bro = 0; bro < n_bro; bro++) {
				tree::nodeP brother = fatherNode->getSon(bro);
				if (brother != mynode) {
					doubleRepMantisa tmp_bro=0.0;
					for(letterInSon=0; letterInSon<pi.alphabetSize();letterInSon++) {
						tmp_bro+=pi.getPij(brother->id(),letter,letterInSon)*
						cup.get(brother->id(),letterInSon);
					}
					brotherTerm *=tmp_bro;
				}
			}
			totalProb = fatherTerm * brotherTerm;
			ssc.set(mynode->id(),letter,totalProb);
		}
	}
}


//use Pij(t) from the stochastic process instead of precomputed probabilities (via the computePijHom class)
void computeDownAlg::fillComputeDown(const tree& et,
					   const sequenceContainer& sc,
					   const size_t pos,
					   const stochasticProcess& sp,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup){
	ssc.allocatePlace(et.getNodesNum(), sp.alphabetSize());
	treeIterTopDownConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		size_t letter, letterInFather, bro, letterInSon;
		if (mynode->isRoot()) {// if root: set all values to 1.0
			for(letter = 0; letter < sp.alphabetSize(); letter++) {
				ssc.set(mynode->id(), letter, 1.0);
			}
			mynode = tIt.next(); //continue
		}
		tree::nodeP fatherNode = mynode->father();
		const size_t n_bro = fatherNode->getNumberOfSons();
		for(letter = 0; letter < sp.alphabetSize(); letter++) {
			doubleRepMantisa totalProb=1.0;
			doubleRepMantisa fatherTerm=0;
			if (fatherNode->isRoot()) 
			{
				fatherTerm = 1.0;
			}
			else
			{
				for(letterInFather = 0; letterInFather < sp.alphabetSize(); letterInFather++)
				{
					MDOUBLE dist = fatherNode->dis2father() * sp.getGlobalRate(); 
					fatherTerm += sp.Pij_t(letter, letterInFather, dist)
					* ssc.get(fatherNode->id(), letterInFather);
				}
			}
			doubleRepMantisa brotherTerm = 1.0;
			for(bro = 0; bro < n_bro; bro++) {
				tree::nodeP brother = fatherNode->getSon(bro);
				if (brother != mynode) {
					doubleRepMantisa tmp_bro=0.0;
					for(letterInSon = 0; letterInSon < sp.alphabetSize(); letterInSon++) 
					{
						MDOUBLE dist = brother->dis2father() * sp.getGlobalRate();
						tmp_bro += sp.Pij_t(letter, letterInSon, dist)
						* cup.get(brother->id(), letterInSon);
					}
					brotherTerm *= tmp_bro;
				}
			}
			totalProb = fatherTerm * brotherTerm;
			ssc.set(mynode->id(), letter, totalProb);
		}
	}
}


//compute probabilities with a site-specific rate
void computeDownAlg::fillComputeDownSpecificRate(const tree& et,
					   const sequenceContainer& sc,
					   const size_t pos,
					   const stochasticProcess& sp,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup,
					   const MDOUBLE gRate){
	ssc.allocatePlace(et.getNodesNum(), sp.alphabetSize());
	treeIterTopDownConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		size_t letter, letterInFather, bro, letterInSon;
		if (mynode->isRoot()) {// if root: set all values to 1.0
			for(letter = 0; letter < sp.alphabetSize(); letter++) {
				ssc.set(mynode->id(), letter, 1.0);
			}
			mynode = tIt.next(); //continue
		}
		tree::nodeP fatherNode = mynode->father();
		const size_t n_bro = fatherNode->getNumberOfSons();
		for(letter = 0; letter < sp.alphabetSize(); letter++) {
			doubleRepMantisa totalProb=1.0;
			doubleRepMantisa fatherTerm=0;
			if (fatherNode->isRoot()) 
			{
				fatherTerm = 1.0;
			}
			else
			{
				for(letterInFather = 0; letterInFather < sp.alphabetSize(); letterInFather++)
				{
					MDOUBLE dist = fatherNode->dis2father() * gRate * sp.getGlobalRate(); 
					fatherTerm += sp.Pij_t(letter, letterInFather, dist)
					* ssc.get(fatherNode->id(), letterInFather);
				}
			}
			doubleRepMantisa brotherTerm = 1.0;
			for(bro = 0; bro < n_bro; bro++) {
				tree::nodeP brother = fatherNode->getSon(bro);
				if (brother != mynode) {
					doubleRepMantisa tmp_bro=0.0;
					for(letterInSon = 0; letterInSon < sp.alphabetSize(); letterInSon++) 
					{
						MDOUBLE dist = brother->dis2father() * gRate * sp.getGlobalRate();
						tmp_bro += sp.Pij_t(letter, letterInSon, dist)
						* cup.get(brother->id(), letterInSon);
					}
					brotherTerm *= tmp_bro;
				}
			}
			totalProb = fatherTerm * brotherTerm;
			ssc.set(mynode->id(), letter, totalProb);
		}
	}
}

// The filled sscGivenRoot is using the "Gam" class (over all rate categories) for placing letter@root hidden state
void computeDownAlg::fillComputeDownNonReversible(const tree& et,
		const sequenceContainer& sc,
		const size_t pos,
		const computePijHom& pi,
		suffStatGlobalGamPos& sscGivenRoot,
		const suffStatGlobalHomPos& cup)
{
			sscGivenRoot.allocatePlace(pi.alphabetSize(),et.getNodesNum(), pi.alphabetSize());
			treeIterTopDownConst tIt(et);
			for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
				size_t letter,letterInFather,bro,letterInSon;
				if (mynode->father()==NULL) {//root
					for (size_t letterAtRoot=0; letterAtRoot<pi.alphabetSize();letterAtRoot++){
						for(letter=0; letter<pi.alphabetSize();letter++) {
							MDOUBLE ind = (letterAtRoot==letter?1.0:0.0);
							sscGivenRoot.set(letterAtRoot,mynode->id(),letter,ind);
						}
					}
					mynode = tIt.next(); //continue
				}
				tree::nodeP fatherNode=mynode->father();
				const size_t n_bro=fatherNode->getNumberOfSons();
				for(size_t letterAtRoot=0; letterAtRoot<pi.alphabetSize();letterAtRoot++) {//root state
					for(letter=0; letter<pi.alphabetSize();letter++) {//letter for current down calc (at father of node)
						doubleRepMantisa totalProb=1.0;
						doubleRepMantisa fatherTerm=0;
						//down of father
						if (fatherNode->father()!=NULL) { // not son of root
							for(letterInFather=0; letterInFather<pi.alphabetSize();letterInFather++)//father of father
								fatherTerm += pi.getPij(fatherNode->id(),letterInFather,letter)*
								sscGivenRoot.get(letterAtRoot,fatherNode->id(),letterInFather);
						}
						else {//son of root
							fatherTerm=(letterAtRoot==letter?1.0:0.0);
						}
						doubleRepMantisa brotherTerm=1.0;
						for(bro = 0; bro < n_bro; bro++) {
							tree::nodeP brother = fatherNode->getSon(bro);
							if (brother != mynode) {
								doubleRepMantisa tmp_bro=0.0;
								for(letterInSon=0; letterInSon<pi.alphabetSize();letterInSon++) {
									tmp_bro+=pi.getPij(fatherNode->getSon(bro)->id(),letter,letterInSon)*
										cup.get(brother->id(),letterInSon);
								}
								brotherTerm *=tmp_bro;
							}
						}
						totalProb = fatherTerm * brotherTerm;
						sscGivenRoot.set(letterAtRoot,mynode->id(),letter,totalProb);
					}
				}
			}			
		}