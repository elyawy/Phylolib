// customDistribution.h - Elya 12.11.2025

#ifndef ___CUSTOM_DIST
#define ___CUSTOM_DIST

/************************************************************ 
This represents a custom discrete distribution with user-specified rates and probabilities.

The distribution allows users to provide arbitrary rate categories and their associated
probabilities, useful for non-parametric rate heterogeneity models.

_rates is a vector with the rate value for each category.
_ratesProb represents the probability of each category.
_globalRate represents a scaling factor applied to all rates.
************************************************************/

#include "definitions.h"
#include "distribution.h"

class customDistribution : public distribution {

public:
	explicit customDistribution(const Vdouble& rates, const Vdouble& ratesProb);
	explicit customDistribution();
	explicit customDistribution(const customDistribution& other);
        
	virtual ~customDistribution() {};

	const size_t categories() const {return _rates.size();}
	// virtual void change_number_of_categories(int in_number_of_categories); // Does nothing for custom distribution
	virtual const MDOUBLE rates(const size_t i) const {return _rates[i]*_globalRate;}
	virtual const MDOUBLE ratesProb(const size_t i) const {return _ratesProb[i];}
	virtual distribution* clone() const { return new customDistribution(*this); }
	virtual void setGlobalRate(const MDOUBLE x) {_globalRate = x;}
 	virtual MDOUBLE getGlobalRate() const {return _globalRate;}

	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const;

	void setCustomParameters(const Vdouble& rates, const Vdouble& ratesProb);

private:	
	Vdouble _rates;
	Vdouble _ratesProb;
	MDOUBLE _globalRate;
};

#endif