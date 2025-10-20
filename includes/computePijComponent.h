// $Id: computePijComponent.h 9253 2011-01-31 01:37:21Z rubi $

#ifndef ___COMPUTE_PIJ_COMPONENT
#define ___COMPUTE_PIJ_COMPONENT

#include <memory>

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "multipleStochasticProcess.h"
#include "gammaDistribution.h"
#include "DiscreteDistribution.h"

class computePijHomSpec {//specific node, no rate variation
public:
	virtual ~computePijHomSpec(){};
	void fillPij(const MDOUBLE dis, const stochasticProcess& sp, int derivationOrder = 0, bool isReversible =true);
	void fillDistributions();

	void resize(const int alphabetSize) {
		_V.resize(alphabetSize);
		for (int z=0;z<alphabetSize;++z) _V[z].resize(alphabetSize);
	}

	size_t alphabetSize() const {return _V.size();}
	MDOUBLE getPij(const int let1,const int let2)const{
		return _V[let1][let2];
	}

	template<typename RngType = std::mt19937_64>
	ALPHACHAR getRandomChar(const ALPHACHAR originLetter, RngType &rng) const {
		return _D[originLetter]->drawSample(rng) - 1;
	}
	VVdouble _V; // let, let
	std::vector<std::shared_ptr<DiscreteDistribution>> _D;
};

class computePijHom {//all nodes, no rate variation
public:
	virtual ~computePijHom(){};
	void fillPij(const tree& et, const stochasticProcess& sp, int derivationOrder = 0, bool isReversible =true);
	size_t alphabetSize() const {return _V[0].alphabetSize();}
	size_t getNodesNum() const {return _V.size();}
	MDOUBLE getPij(const int nodeId,const int let1,const int let2)const{
		return _V[nodeId].getPij(let1,let2);
	}

template<typename RngType = std::mt19937_64>
	ALPHACHAR getRandomChar(const int nodeId, const ALPHACHAR originLetter, RngType &rng) const {
		return _V[nodeId].getRandomChar(originLetter, rng);
	}

	vector<computePijHomSpec> _V; // let, let
};

class computePijGam {//
public:
	virtual ~computePijGam(){};
	void fillPij(const tree& et, const stochasticProcess& sp, int derivationOrder = 0, bool isReversible =true);
	size_t categories() const {return _V.size();}
	size_t alphabetSize() const {return _V[0].alphabetSize();}
	size_t getNodesNum() const {return _V[0].getNodesNum();}

	MDOUBLE getPij(const int rateCategor,const int nodeId,const int let1,const int let2)const{
		return _V[rateCategor].getPij(nodeId,let1,let2);
	}

	template<typename RngType = std::mt19937_64>
	ALPHACHAR getRandomChar(const int rateCategor,const int nodeId, const ALPHACHAR originLetter, RngType &rng) const {
		return _V[rateCategor].getRandomChar(nodeId,originLetter, rng);
	}
	computePijHom& operator[] (int i) {return _V[i];}
	const computePijHom& operator[] (int i) const {return _V[i];}
	vector<computePijHom> _V; // each rate category
};


class PijGamMultipleStochasticProcesses
{
public:
	virtual ~PijGamMultipleStochasticProcesses() {};
	void fillPij(const tree& et, const vector<stochasticProcess>& V_sp);
	computePijGam& operator[] (size_t i) { return _V[i]; }
	const computePijGam& operator[] (int i) const { return _V[i]; }
	vector<computePijGam> _V; // each stochastic process

private:

};




#endif
