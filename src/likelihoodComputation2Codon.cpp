
#include "../includes/likelihoodComputation2Codon.h"

#include "../includes/wYangModel.h"
#include "../includes/definitions.h"
#include "../includes/tree.h"
#include "../includes/computeUpAlg.h"
#include "../includes/likelihoodComputation.h"

#include <cmath>
#include <cassert>

using namespace likelihoodComputation2Codon;



MDOUBLE likelihoodComputation2Codon::getTreeLikelihoodAllPosAlphTheSame(const tree& et,
							const sequenceContainer& sc,
							const seqContainerTreeMap &sctm,
							const vector<stochasticProcess>& spVec,const distribution * distr){
	computePijGam pi;
	pi._V.resize(distr->categories());
	for (size_t i=0; i < spVec.size(); ++i) {
		pi._V[i].fillPij(et,spVec[i]);
	}
	
	suffStatGlobalGam ssc;
	computeUpAlg cup;
	cup.fillComputeUp(et,sc,sctm,pi,ssc);
	
	MDOUBLE res = 0.0;
	for (size_t k=0; k < sc.seqLen(); ++k) {
		MDOUBLE lnL = log(likelihoodComputation2Codon::getProbOfPosUpIsFilledSelectionGam(k,//pos,
			  et,//const tree& 
			  sc,// sequenceContainer& sc,
			  spVec[0],
			  ssc[k],//const computePijGam& ,
			  distr)); //W distribution ,
		LOG(20,<<"pos= "<<k<<" lnL= "<<lnL<<endl);
		res += lnL;
		//if (k==5) exit(0);
			  
	}
	return res;
	


}


MDOUBLE likelihoodComputation2Codon::getProbOfPosUpIsFilledSelectionGam(const int pos,const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGamPos& cup,const distribution * distr){

	doubleRepMantisa tmp=0.0;
	for (size_t categor = 0; categor < distr->categories(); ++categor) {
		doubleRepMantisa veryTmp =0;
		for (size_t let =0; let < sc.alphabetSize(); ++let) {
			veryTmp+=cup.get(categor,et.getRoot()->id(),let) * sp.freq(let);
			
		}
		//cout<<"category= "<<categor<<" fh= "<<veryTmp<<" freqCategor= "<<distr->ratesProb(categor)<<endl;
		tmp += veryTmp*distr->ratesProb(categor);
	}
	assert(tmp>0.0);
	return convert(tmp);
}

MDOUBLE likelihoodComputation2Codon::getTreeLikelihoodFromUp2(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGam& cup,
						Vdouble& posLike, // fill this vector with each position likelihood but without the weights.
						const distribution * distr,
						const Vdouble * weights) {
	posLike.clear();
	MDOUBLE like = 0;
	//computing the likelihood from up:
	for (size_t pos = 0; pos < sc.seqLen(); ++pos) {
		doubleRepMantisa tmp=0;
		for (size_t categor = 0; categor < distr->categories(); ++categor) {
			doubleRepMantisa veryTmp =0;
			for (size_t let =0; let < sc.alphabetSize(); ++let) {
				veryTmp+=cup.get(pos,categor,et.getRoot()->id(),let) * sp.freq(let);
			}
			tmp += veryTmp*distr->ratesProb(categor);
		}
		assert(tmp>0.0);
		like += log(tmp) * (weights?(*weights)[pos]:1);
		posLike.push_back(convert(tmp));
	}
	return like;

}
