// $Id: aaJC.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___AA_JC
#define ___AA_JC

#include "replacementModel.h"
#include <cmath>
using namespace std;

namespace aaDef {
	const MDOUBLE Alp = 20.0;
	const MDOUBLE odAl = 1.0/Alp; // one divided by alphabet
	const MDOUBLE om_odAl = 1.0-odAl; // one minus odAl;
	const MDOUBLE alDiv_omalp = Alp/(Alp-1.0);
	const MDOUBLE m_alDiv_omalp = -alDiv_omalp;
}

class aaJC : public replacementModel {
private:
	mutable MDOUBLE cached_exponent;
	mutable MDOUBLE cached_d;
public:

	virtual replacementModel* clone() const { return new aaJC(*this); }// see note down:
//	virtual aaJC* clone() const { return new aaJC(*this); } 
	const size_t alphabetSize() const {return 20;}

	explicit aaJC(){};
	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const {
		if (d != cached_d) {//could be a bug if d==0 in the first call
			cached_exponent = exp(aaDef::m_alDiv_omalp*d);
			cached_d = d;
		}
//(wrong!)		return ((i==j) ?  0.05+0.95*exp(-20.0*d): 0.05-0.05*exp(-20.0*d));
		return ((i==j) ?  aaDef::odAl*(1+cached_exponent): aaDef::odAl*(1-cached_exponent));

	}

	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const{
		//(worng!)return ((i==j) ?  -19.0*exp(-20.0*d): exp(-20.0*d));
		if (d != cached_d) {
			cached_exponent = exp(aaDef::m_alDiv_omalp*d);
			cached_d = d;
		}
		return ((i==j) ?  -cached_exponent: cached_exponent/(aaDef::Alp-1));
	}
	const MDOUBLE freq(const size_t i) const {return aaDef::odAl;};

	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const{
	//(wrong!)	return ((i==j) ?  19.0*20.0*exp(-20.0*d): 0.0-20.0*exp(-20.0*d));
		if (d != cached_d) {
			cached_exponent = exp(aaDef::m_alDiv_omalp*d);
			cached_d = d;
		}
		return ((i==j) ?  aaDef::alDiv_omalp*cached_exponent: aaDef::m_alDiv_omalp*cached_exponent);
	}

};

#endif 

// note: according to the new C++ rules, the clone function should be like this:
//	virtual aaJC* clone() const { return new aaJC(*this); }
// however, not all compiler support it yet. look at More Effective C++ page 126.



