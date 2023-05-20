#ifndef _MULTIPLE_STOCHASTIC_PROCESS
#define _MULTIPLE_STOCHASTIC_PROCESS

#include "stochasticProcess.h"


class multipleStochasticProcess {
public:
	multipleStochasticProcess();
	virtual ~multipleStochasticProcess();
	virtual MDOUBLE getProb(size_t spPlace) const;
	virtual stochasticProcess* getSp(size_t spPlace);
	virtual size_t getSPVecSize() const {return _spVec.size();}
	virtual void setSpVec(vector<stochasticProcess>& spVec);

    
protected:
	virtual void copy(const multipleStochasticProcess * pOther);
protected:
	vector<stochasticProcess> _spVec;
	Vdouble _spProb;
};
#endif
