// $Id: pijAccelerator.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___PIJ_ACCELERATOR
#define ___PIJ_ACCELERATOR

#include "definitions.h"
#include "replacementModel.h"

class pijAccelerator {
public:
	virtual pijAccelerator* clone() const = 0;
	virtual ~pijAccelerator() = 0;

	virtual const MDOUBLE Qij(const int i, const int j) const = 0;

	virtual const MDOUBLE Pij_t(const int i, const int j, const MDOUBLE t) const = 0;
	virtual void fillPijRow(const int i, const double t, double* outRow) const {
		// Default fallback loop so other accelerators compile out-of-the-box
		// Alphabet size can be fetched via alphabetSize() or the derived class properties
		const size_t alphSize = alphabetSize();
		for (size_t j = 0; j < alphSize; ++j) {
			outRow[j] = Pij_t(i, j, t);
		}
	}
	virtual const MDOUBLE freq(const size_t i) const = 0;	// P(i)
	virtual const MDOUBLE dPij_dt(const int i, const int j, const MDOUBLE t) const =0;
	virtual const MDOUBLE d2Pij_dt2(const int i, const int j, const MDOUBLE t) const =0;
	virtual replacementModel* getReplacementModel() const =0; // @@@@ this const is a lie !!!
	virtual const size_t alphabetSize() const =0;
};





#endif

