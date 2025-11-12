// customDistribution.cpp - Elya 12.11.2025

#include "../includes/customDistribution.h"
#include "../includes/errorMsg.h"
#include <algorithm>

customDistribution::customDistribution(const Vdouble& rates, const Vdouble& ratesProb) 
	: distribution() {
	_globalRate = 1.0;
	setCustomParameters(rates, ratesProb);
}

customDistribution::customDistribution() : distribution() {
	_globalRate = 1.0;
}

// Copy constructor
customDistribution::customDistribution(const customDistribution& other) : 
	_rates(other._rates),
	_ratesProb(other._ratesProb),
	_globalRate(other._globalRate)
{
}

void customDistribution::setCustomParameters(const Vdouble& rates, const Vdouble& ratesProb) {
	if (rates.size() != ratesProb.size()) {
		errorMsg::reportError("customDistribution: rates and probabilities must have the same size");
	}
	
	if (rates.empty()) {
		errorMsg::reportError("customDistribution: rates and probabilities cannot be empty");
	}
	
	// Validate probabilities sum to 1.0 (within tolerance)
	MDOUBLE sum = 0.0;
	for (size_t i = 0; i < ratesProb.size(); ++i) {
		if (ratesProb[i] < 0.0) {
			errorMsg::reportError("customDistribution: probabilities must be non-negative");
		}
		sum += ratesProb[i];
	}
	
	if (std::abs(sum - 1.0) > 1e-6) {
		errorMsg::reportError("customDistribution: probabilities must sum to 1.0");
	}
	
	_rates = rates;
	_ratesProb = ratesProb;
}

const MDOUBLE customDistribution::getCumulativeProb(const MDOUBLE x) const {
	// Return sum of probabilities for all categories with rate <= x
	MDOUBLE cumulativeProb = 0.0;
	
	for (size_t i = 0; i < _rates.size(); ++i) {
		if (_rates[i] * _globalRate <= x) {
			cumulativeProb += _ratesProb[i];
		}
	}
	
	return cumulativeProb;
}

// void customDistribution::change_number_of_categories(int in_number_of_categories) {
// 	// Do nothing - custom distribution has fixed user-provided categories
// 	// Silently ignore this call
// 	return;
// }