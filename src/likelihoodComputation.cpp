
#include "../includes/definitions.h"
#include "../includes/tree.h"
#include "../includes/computeUpAlg.h"
#include "../includes/likelihoodComputation.h"
#include "../includes/gammaUtilities.h"
#include <cmath>
#include <cassert>
using namespace likelihoodComputation;
using namespace std;

/********************************************************************************************
likelihood computation - full data (1)
*********************************************************************************************/
MDOUBLE likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(const tree& et,
																  const sequenceContainer& sc,
																  const stochasticProcess& sp,
																  const Vdouble * const weights,
																  unObservableData *unObservableData_p)
{
	computePijGam pi;
	pi.fillPij(et,sp);
	seqContainerTreeMap sctm(sc, et);
	MDOUBLE logLforMissingData;
	MDOUBLE LforMissingData;
	if(unObservableData_p){
		logLforMissingData = unObservableData_p->getlogLforMissingData();
		LforMissingData = exp(logLforMissingData);
	}	
	MDOUBLE res =0;
	doubleRepMantisa LofPos;
	for (size_t k=0; k < sc.seqLen(); ++k) {
		LofPos = likelihoodComputation::getLofPos(k,//pos,
			et,		//const tree& 
			sc,		// sequenceContainer& sc,
			sctm,
			pi,		//const computePijGam& ,
			sp,
			NULL);
		if(unObservableData_p){		// conditioning on observability for all rateCat.
			LofPos = LofPos / (1- LforMissingData);
		}
		res += log(LofPos) * (weights?(*weights)[k]:1);//const stochasticProcess& );
	}
	//if(unObservableData_p){		// conditioning on observability for allPos & allRateCat
	//	res = res - sc.seqLen()*log(1- exp(unObservableData_p->getlogLforMissingData()));
	//}
	return res;
}

/********************************************************************************************
likelihood computation - per pos (1.1)
*********************************************************************************************/
doubleRepMantisa likelihoodComputation::getLofPos(const int pos,
										   const tree& et,
										   const sequenceContainer& sc,
											const seqContainerTreeMap & sctm,
										   const computePijGam& pi,
										   const stochasticProcess& sp,
										   unObservableData *unObservableData_p)
{
	//	with the pi already computed.
	doubleRepMantisa tmp=0;
	size_t numOfCat = sp.categories();
	vector<doubleRepMantisa> tmpPerCat;
	tmpPerCat.resize(numOfCat);	

	for (size_t i=0; i < sp.categories();++i) {
		tmpPerCat[i] = getLofPos(pos,et,sc, sctm,pi[i],sp);
		tmp += tmpPerCat[i]*sp.ratesProb(i);
	}
	if(unObservableData_p){
		tmp = tmp / (1- exp(unObservableData_p->getlogLforMissingData()));
	}
	return tmp;
}

/********************************************************************************************
likelihood computation - per pos, per cat (1.1.1)
*********************************************************************************************/
doubleRepMantisa likelihoodComputation::getLofPos(const int pos,
					  const tree& et,
					  const sequenceContainer& sc,
					  const seqContainerTreeMap & sctm,
					  const computePijHom& pi,
					  const stochasticProcess& sp,
					  unObservableData *unObservableData_p)
{
	computeUpAlg cup;
	suffStatGlobalHomPos ssc;
	cup.fillComputeUp(et,sc,pos,pi, sctm,ssc);

	doubleRepMantisa tmp = 0.0;
	for (size_t let = 0; let < sp.alphabetSize(); ++let) {
		doubleRepMantisa tmpLcat=
				ssc.get(et.getRoot()->id(),let)*
				sp.freq(let);
		if (!DBIG_EQUAL(convert(tmpLcat), 0.0))
		{
			cerr<<"tmpLcat = "<<tmpLcat<<endl;
			errorMsg::reportError("error in likelihoodComputation::getLofPos. likelihood is smaller than zero");
		}		
		//assert(tmpLcat>=0.0);
		tmp+=tmpLcat;
	}
//	cout<<"likelihoodComputation::getLofPos: tmp = "; tmp.outputn(cout);	// DEBUG EP
	if (!DBIG_EQUAL(convert(tmp), 0.0)){
		LOG(5,<<"likelihoodComputation::getLofPos: "<< tmp<<endl;);
		LOG(5,<<"pos = "<< pos <<endl;);
		tmp = EPSILON;
		//errorMsg::reportError("likelihoodComputation::getLofPos: likelihood of pos was zero!",1);
	}

	if(unObservableData_p){		// conditioning on observability
		tmp = tmp / (1- exp(unObservableData_p->getlogLforMissingData()));
	}
	return tmp;
}

//r4s_proportional
/********************************************************************************************
likelihood computation - full data (1)
*********************************************************************************************/
Vdouble likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(const tree& et,
																  const vector<sequenceContainer>& sc,
																  multipleStochasticProcess* msp,
																  const gammaDistribution* pProportionDist,
																  const Vdouble * const weights)
{
	Vdouble geneLikelihoodVec;
	//geneRateLikelihoodVec[geneN][globalRateCateg] will hold the LL of the gene given the global rate
	VVdouble geneRateLikelihoodVec;
	geneLikelihoodVec.resize(sc.size(),0.0);
	geneRateLikelihoodVec.resize(sc.size());
	for(size_t geneN = 0;geneN < sc.size();++geneN){
		seqContainerTreeMap sctm(sc[geneN], et);
		geneRateLikelihoodVec[geneN].resize(pProportionDist->categories(),0.0);
		for(size_t globalRateCateg = 0;globalRateCateg < pProportionDist->categories();++globalRateCateg){
			msp->getSp(geneN)->setGlobalRate(pProportionDist->rates(globalRateCateg));
			computePijGam pi;
			pi.fillPij(et,*msp->getSp(geneN));
			doubleRepMantisa LofPos;
			for (size_t k=0; k < sc[geneN].seqLen(); ++k) {
				//LofPos is sum LofPos_LocalRateCat_i*p(LocalRateCat_i)
				LofPos = likelihoodComputation::getLofPosProportional(k,//pos,
					et,		//const tree& 
					sc[geneN],		// sequenceContainer& sc,
					sctm,
					pi,		//const computePijGam& ,
					*msp->getSp(geneN)); //removed the prior of the globar rate categ cause it is multiplied below
				geneRateLikelihoodVec[geneN][globalRateCateg] += log(LofPos)*(weights?(*weights)[k]:1);
			}
		}
		//Once we are finished iterating over all globalRateCategs we need to sum the log likelihood for this gene
		//which is: log(prior(globalRateCateg_i)*exp(geneRateLikelihoodVec[geneN][globalRateCateg_i]+prior(globalRateCateg_j)*exp(geneRateLikelihoodVec[geneN][globalRateCateg_j]..)
		//assuming a flat prior this equals: log(prior(globalRateCateg))+log(exp(geneRateLikelihoodVec[geneN][globalRateCateg_i]+exp(geneRateLikelihoodVec[geneN][globalRateCateg_j]..)
		//which can be written as:log(prior(globalRateCateg))+log(exp(geneRateLikelihoodVec[geneN][globalRateCateg_i]))(1+exp(geneRateLikelihoodVec[geneN][globalRateCateg_j]-geneRateLikelihoodVec[geneN][globalRateCateg_i]..)
        geneLikelihoodVec[geneN] = log(pProportionDist->ratesProb(0))+exponentResolver(geneRateLikelihoodVec[geneN]);//Strictly assumes a flat prior distribution
	}
	return geneLikelihoodVec;
}

/********************************************************************************************
likelihood computation - per pos (1.1)
*********************************************************************************************/
//Old - remove when QA is done
doubleRepMantisa likelihoodComputation::getLofPosProportional(const int pos,
										   const tree& et,
										   const sequenceContainer& sc,
										   const seqContainerTreeMap & sctm,
										   const computePijGam& pi,
										   const stochasticProcess& sp,
										   const MDOUBLE globalRateProb)
{
	//	with the pi already computed.
	doubleRepMantisa tmp=0;
	int numOfCat = sp.categories();
	vector<doubleRepMantisa> tmpPerCat;
	tmpPerCat.resize(numOfCat);	

	for (size_t i=0; i < sp.categories();++i) {
		tmpPerCat[i] = getLofPos(pos,et,sc, sctm,pi[i],sp);
		tmp += tmpPerCat[i]*sp.ratesProb(i)*globalRateProb; //old - now globalRateProb is multipled outside
	}
	return tmp;
}

/********************************************************************************************
likelihood computation - per pos (1.1)
*********************************************************************************************/
doubleRepMantisa likelihoodComputation::getLofPosProportional(const int pos,
										   const tree& et,
										   const sequenceContainer& sc,
										   const seqContainerTreeMap& sctm,
										   const computePijGam& pi,
										   const stochasticProcess& sp)
{
	//	with the pi already computed.
	doubleRepMantisa tmp=0;
	int numOfCat = sp.categories();
	vector<doubleRepMantisa> tmpPerCat;
	tmpPerCat.resize(numOfCat);	

	for (size_t i=0; i < sp.categories();++i) {
		tmpPerCat[i] = getLofPos(pos,et,sc, sctm,pi[i],sp);
		tmp += tmpPerCat[i]*sp.ratesProb(i);
	}
	return tmp;
}

//r4s_proportional


/********************************************************************************************
*********************************************************************************************/
doubleRepMantisa likelihoodComputation::getProbOfPosWhenUpIsFilledHom(const int pos,
					  const tree& et,
					  const sequenceContainer& sc,
					  const stochasticProcess& sp,
					  const suffStatGlobalHomPos& ssc){
// using the pij of stochastic process rather than pre computed pij's...
	if (ssc.size()==0) {errorMsg::reportError("error in function likelihoodComputation::getLofPosWhenUpIsFilled");}
	doubleRepMantisa tmp = 0.0;
	for (size_t let = 0; let < sp.alphabetSize(); ++let) {
		doubleRepMantisa tmpLcat=
				ssc.get(et.getRoot()->id(),let)*
				sp.freq(let);
		tmp+=tmpLcat;
	}
	return tmp;
}

/********************************************************************************************
*********************************************************************************************/
doubleRepMantisa likelihoodComputation::getLofPosHomModelEachSiteDifferentRate(const int pos,
					  const tree& et,
					  const sequenceContainer& sc,
					  const seqContainerTreeMap & sctm,
					  const stochasticProcess& sp){
// using the pij of stochastic process rather than pre computed pij's...
	if (sp.categories()!=1) {
		  errorMsg::reportError("num of categories in function getLofPosHomModel must be one");
	}
	computeUpAlg cup;
	suffStatGlobalHomPos ssc;
	computePijHom cpij;
	cpij.fillPij(et,sp);
	cup.fillComputeUp(et,sc,pos,cpij, sctm,ssc);
	return getProbOfPosWhenUpIsFilledHom(pos,et,sc,sp,ssc);
}
/********************************************************************************************
*********************************************************************************************/
doubleRepMantisa likelihoodComputation::getLofPosGamModelEachSiteDifferentRate(const int pos,
					  const tree& et,
					  const sequenceContainer& sc,
					  const seqContainerTreeMap & sctm,
					  const stochasticProcess& sp){
	computePijGam pi;
	pi.fillPij(et,sp);
	return getLofPos(pos,et,sc, sctm,pi,sp);
}
/********************************************************************************************
*********************************************************************************************/
doubleRepMantisa likelihoodComputation::getLofPos(const int pos,
					  const tree& et,
					  const sequenceContainer& sc,
					  const stochasticProcess& sp,
					  const MDOUBLE gRate){ // when there is a global rate for this position
// using the pij of stochastic process rather than pre computed pij's...
	computeUpAlg cup;
	suffStatGlobalHomPos ssc;
	cup.fillComputeUpSpecificGlobalRate(et,sc,pos,sp,ssc,gRate);

	doubleRepMantisa tmp = 0.0;
	for (size_t let = 0; let < sp.alphabetSize(); ++let) {
		doubleRepMantisa tmpLcat=
			ssc.get(et.getRoot()->id(),let)*
			sp.freq(let);;
		assert(tmpLcat>=0.0);
		tmp+=tmpLcat;
	}
	return tmp;
}

/********************************************************************************************
*********************************************************************************************/
doubleRepMantisa likelihoodComputation::getLofPosAndPosteriorOfRates(const int pos,
															  const tree& et,
															  const sequenceContainer& sc,
															  const seqContainerTreeMap& sctm,
															  const computePijGam& pi,
															  const stochasticProcess& sp,
															  vector<doubleRepMantisa>& postrior){
//	with the pi already computed.
	doubleRepMantisa tmp=0;
	for (size_t i=0; i < sp.categories();++i) {
	  postrior[i]=getLofPos(pos,et,sc, sctm,pi[i],sp)*sp.ratesProb(i);
	  tmp += postrior[i]; 
	}
	for (size_t i=0; i < sp.categories();++i)
	  postrior[i] /= tmp;
	return tmp;
}
/********************************************************************************************
*********************************************************************************************/
MDOUBLE likelihoodComputation::getTreeLikelihoodFromUp(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGam& cup,
						const Vdouble * weights) {
	MDOUBLE like = 0;
	//computing the likelihood from up:
	for (size_t pos = 0; pos < sc.seqLen(); ++pos) {
		doubleRepMantisa tmp=0;
		for (size_t categor = 0; categor < sp.categories(); ++categor) {
			doubleRepMantisa veryTmp =0;
			for (size_t let =0; let < sc.getAlphabet()->size(); ++let) {
				veryTmp+=cup.get(pos,categor,et.getRoot()->id(),let) * sp.freq(let);
			}
			tmp += veryTmp*sp.ratesProb(categor);
		}
		like += log(tmp) * (weights?(*weights)[pos]:1);
	}
	return like;
}
/********************************************************************************************
*********************************************************************************************/
MDOUBLE likelihoodComputation::getTreeLikelihoodFromUp2(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGam& cup,
						vector<doubleRepMantisa>& posLike, // fill this vector with each position likelihood but without the weights.
						const Vdouble * weights,
						unObservableData* unObservableData_p) {
	posLike.clear();
	MDOUBLE like = 0;
	//computing the likelihood from up:
	for (size_t pos = 0; pos < sc.seqLen(); ++pos) {
		doubleRepMantisa tmp=0;
		for (size_t categor = 0; categor < sp.categories(); ++categor) {
			doubleRepMantisa veryTmp =0;
			for (size_t let =0; let < sc.alphabetSize(); ++let) {
				veryTmp+=cup.get(pos,categor,et.getRoot()->id(),let) * sp.freq(let);
			}
			tmp += veryTmp*sp.ratesProb(categor);
		}
		assert(tmp>0.0);
		if(unObservableData_p){
			tmp = tmp/(1- exp(unObservableData_p->getlogLforMissingData()));
		}
		like += log(tmp) * (weights?(*weights)[pos]:1);
		posLike.push_back(tmp);
	}
	return like;
}
/********************************************************************************************
*********************************************************************************************/
//old
MDOUBLE likelihoodComputation::getTreeLikelihoodFromUp2(const tree& et,
						const sequenceContainer& sc,
						stochasticProcess& sp,
						const suffStatGlobalGamProportional& cup,
						const gammaDistribution* pProportionDist,
						vector<doubleRepMantisa>& posLike, // fill this vector with each position likelihood but without the weights.
						const Vdouble * weights) {
	posLike.clear();
	MDOUBLE like = 0.0;
	//computing the likelihood from up:
	for (size_t pos = 0; pos < sc.seqLen(); ++pos) {
		doubleRepMantisa tmp(0.0);
		for(size_t globalRateCategor = 0;globalRateCategor < pProportionDist->categories();++globalRateCategor){
			for (size_t localRateCategor = 0; localRateCategor < sp.categories(); ++localRateCategor) {
				doubleRepMantisa veryTmp =0;
				for (size_t let =0; let < sc.alphabetSize(); ++let) {
					veryTmp+=cup.get(pos,globalRateCategor,localRateCategor,et.getRoot()->id(),let) * sp.freq(let);
				}
				tmp += veryTmp*pProportionDist->ratesProb(globalRateCategor)*sp.ratesProb(localRateCategor);
			}
		}
		assert(tmp>0.0);
		like += log(tmp) * (weights?(*weights)[pos]:1);		
		posLike.push_back(tmp);
	}
	return like;
}

//new
MDOUBLE likelihoodComputation::getTreeLikelihoodFromUp2(const tree& et,
						const sequenceContainer& sc,
						stochasticProcess& sp,
						const suffStatGlobalGamProportional& cup,
						const gammaDistribution* pProportionDist,
						vector<vector<doubleRepMantisa> >& posLike,
						const Vdouble * weights) {
	for(size_t pos = 0;pos < sc.seqLen();++pos){
		posLike[pos].resize(pProportionDist->categories(),0.0);
	}
	Vdouble geneRateLikelihoodVec;
	geneRateLikelihoodVec.resize(pProportionDist->categories(),0.0);
	MDOUBLE like = 0.0;
	//computing the likelihood from up:
	for (size_t pos = 0; pos < sc.seqLen(); ++pos) {
		vector<doubleRepMantisa> tmpVec; //hold the LofPos for each global rate category
		tmpVec.resize(pProportionDist->categories(),0.0);//This would sum for every global rate category
		for(size_t globalRateCategor = 0;globalRateCategor < pProportionDist->categories();++globalRateCategor){
			doubleRepMantisa tmp1(0.0);
			doubleRepMantisa tmp2(0.0);
			for (size_t localRateCategor = 0; localRateCategor < sp.categories(); ++localRateCategor) {
				doubleRepMantisa veryTmp(0.0);
				for (size_t let =0; let < sc.alphabetSize(); ++let) {
					veryTmp+=cup.get(pos,globalRateCategor,localRateCategor,et.getRoot()->id(),let) * sp.freq(let);
				}
				tmp1 += veryTmp;
				tmp2 += veryTmp*sp.ratesProb(localRateCategor);
			}
			tmpVec[globalRateCategor] += tmp2;
			posLike[pos][globalRateCategor] = tmp1;
		}
		for(size_t globalRateCategor = 0;globalRateCategor < pProportionDist->categories();++globalRateCategor){
			assert(tmpVec[globalRateCategor]>0.0);
			geneRateLikelihoodVec[globalRateCategor] += log(tmpVec[globalRateCategor])*(weights?(*weights)[pos]:1);
		}
	}
	like = log(pProportionDist->ratesProb(0))+exponentResolver(geneRateLikelihoodVec);
	return like;
}

/********************************************************************************************
 fill the posteriorLike matrix with each position posterior rate (p(r|D))
 but without the weights.
*********************************************************************************************/
MDOUBLE likelihoodComputation::getPosteriorOfRates(const tree& et,
						const sequenceContainer& sc,
						const seqContainerTreeMap &sctm,
						const stochasticProcess& sp,
						vector<vector<doubleRepMantisa> >& posteriorLike, 
						const Vdouble * weights) {
	suffStatGlobalGam cup;
	computeUpAlg cupAlg;
	computePijGam cpGam;
	cpGam.fillPij(et,sp);
	cupAlg.fillComputeUp(et,sc,sctm,cpGam,cup);
	return getPosteriorOfRates(et,sc,sctm,sp,cup,posteriorLike,weights);
}

// fill the posteriorLike matrix with each position posterior rate (p(r|D))
// but without the weights.
MDOUBLE likelihoodComputation::getPosteriorOfRates(const tree& et,
						const sequenceContainer& sc,
						const seqContainerTreeMap &sctm,
						const stochasticProcess& sp,
						const suffStatGlobalGam& cup,
						vector<vector<doubleRepMantisa> >& posteriorLike, 
						const Vdouble * weights) {
	posteriorLike.clear();
	posteriorLike.resize(sc.seqLen());
	for (size_t z=0; z < posteriorLike.size(); ++z) posteriorLike[z].resize(sp.categories());
	MDOUBLE like = 0;
	//computing the likelihood from up:
	for (size_t pos = 0; pos < sc.seqLen(); ++pos) {
		doubleRepMantisa posProb=0;
		for (size_t categor = 0; categor < sp.categories(); ++categor) {
			doubleRepMantisa veryTmp =0;
			for (size_t let =0; let < sc.getAlphabet()->size(); ++let) {
				veryTmp+=cup.get(pos,categor,et.getRoot()->id(),let) * sp.freq(let);
			}
			posProb += veryTmp*sp.ratesProb(categor);
			posteriorLike[pos][categor] += veryTmp*sp.ratesProb(categor);
		}
		like += log(posProb) * (weights?(*weights)[pos]:1);
		for (size_t categor1 = 0; categor1 < sp.categories(); ++categor1) {
			posteriorLike[pos][categor1] /= posProb;
		}
	}

	return like;
}


// fill the posteriorLike matrix with each position posterior rate (p(r|D))
// and the LLPP,  but without the weights.
MDOUBLE likelihoodComputation::getPosteriorOfRatesAndLLPP(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGam& cup,
						vector<vector<doubleRepMantisa> >& posteriorLike, 
						vector<doubleRepMantisa>& LLPerPos, 
						const Vdouble * weights) {
	posteriorLike.clear();
	posteriorLike.resize(sc.seqLen());
	for (size_t z=0; z < posteriorLike.size(); ++z) posteriorLike[z].resize(sp.categories());
	MDOUBLE like = 0;
	//computing the likelihood from up:
	for (size_t pos = 0; pos < sc.seqLen(); ++pos) {
	  LLPerPos[pos] = 0.0;
		for (size_t categor = 0; categor < sp.categories(); ++categor) {
			doubleRepMantisa veryTmp =0;
			for (size_t let =0; let < sc.getAlphabet()->size(); ++let) {
				veryTmp+=cup.get(pos,categor,et.getRoot()->id(),let) * sp.freq(let);
			}
			LLPerPos[pos] += veryTmp*sp.ratesProb(categor);
			posteriorLike[pos][categor] += veryTmp*sp.ratesProb(categor);
		}
		like += log(LLPerPos[pos]) * (weights?(*weights)[pos]:1);
		for (size_t categor1 = 0; categor1 < sp.categories(); ++categor1) {
			posteriorLike[pos][categor1] /= LLPerPos[pos];
		}
	}

	return like;
}

// this function forces non gamma computation of likelihoods from up.
// i.e., even if the stochastic process is really gamma - the likelihood is computed as if there's no gamma.
MDOUBLE likelihoodComputation::getTreeLikelihoodFromUpSpecifcRates(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalHom& cup,
						vector<doubleRepMantisa>& posLike, // fill this vector with each position likelihood but without the weights.
						const Vdouble * weights)
{
	posLike.clear();
	MDOUBLE like = 0;
	//computing the likelihood from up:
	for (size_t pos = 0; pos < sc.seqLen(); ++pos)
	{
		doubleRepMantisa tmp=0;
		for (size_t let =0; let < sc.getAlphabet()->size(); ++let) {
			tmp += cup.get(pos, et.getRoot()->id(), let) * sp.freq(let);
		}
		
		assert(tmp > 0);
		like += log(tmp) * (weights?(*weights)[pos]:1);
		posLike.push_back(tmp);
	}
	return like;	
}
/********************************************************************************************
*********************************************************************************************/
doubleRepMantisa likelihoodComputation::getProbOfPosWhenUpIsFilledGam(const int pos,
						const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGamPos& cup) {
	doubleRepMantisa tmp=0;
	for (size_t categor = 0; categor < sp.categories(); ++categor) {
		doubleRepMantisa veryTmp =0;
		for (size_t let =0; let < sc.alphabetSize(); ++let) {
			veryTmp+=cup.get(categor,et.getRoot()->id(),let) * sp.freq(let);
		}
		tmp += veryTmp*sp.ratesProb(categor);
	}
	assert(tmp>0.0);
	return tmp;
}
/********************************************************************************************
*********************************************************************************************/
MDOUBLE likelihoodComputation::computeLikelihoodAndLikelihoodPerPosition(	const sequenceContainer &sc,
																			const tree &et, 
																			const stochasticProcess &sp, 
																			Vdouble &LLPerPos) {
	MDOUBLE treeLogLikelihood = 0.0;
	computePijGam cpij;
	cpij.fillPij(et, sp);
	seqContainerTreeMap sctm(sc, et);
	LLPerPos.resize(sc.seqLen());
	doubleRepMantisa LofPos;
	for (size_t pos=0; pos < sc.seqLen() ;++pos) {
		LofPos = likelihoodComputation::getLofPos(pos, et, sc, sctm, cpij, sp);
		MDOUBLE tmpLL = log(LofPos);
		treeLogLikelihood += tmpLL;
		LLPerPos[pos] = tmpLL;
	}
	return treeLogLikelihood;
}
/********************************************************************************************
likelihood for each category - used for unObservableData
*********************************************************************************************/
Vdouble likelihoodComputation::getLofPosPerCat(const int pos,
									const tree& et,
									const sequenceContainer& sc,
									const seqContainerTreeMap& sctm,
									const computePijGam& pi,
									const stochasticProcess& sp)
{
//	with the pi already computed.
    int numOfCat = sp.categories();
	Vdouble tmp;
	tmp.resize(numOfCat);
	for (int i=0; i < numOfCat;++i) {
		tmp[i] = convert(getLofPos(pos,et,sc, sctm,pi[i],sp))*sp.ratesProb(i);
	}
	return tmp;
}

//doubleRepMantisa likelihoodComputation::getLofPos(const int pos,
//										   const tree& et,
//										   const sequenceContainer& sc,
//										   const computePijGam& pi,
//										   const stochasticProcess& sp){
////	with the pi already computed.
//	doubleRepMantisa tmp=0;
//	for (int i=0; i < sp.categories();++i) {
//		tmp += getLofPos(pos,et,sc,pi[i],sp)*sp.ratesProb(i);
//	}
//	return tmp;
//}

// MDOUBLE likelihoodComputation::getTreeLikelihoodFromPosteriorAndAlpha(const MDOUBLE alpha,
// 																	  const Vdouble originalBounderi,
// 																	  const VVdouble& posteriorLike,
// 																	  const vector<doubleRepMantisa>& LLPP,
// 																	  const Vdouble* weights)
// {
//   int nCategories = originalBounderi.size()-1;
//   Vdouble rateWeights; rateWeights.resize(nCategories);
//   for (int i=0; i<n; ++i) 
// 	rateWeights[i]=(gammp(alpha, originalBounderi[i+1]*alpha)-gammp(alpha, originalBounderi[i]*alpha))*nCategories;

// }
