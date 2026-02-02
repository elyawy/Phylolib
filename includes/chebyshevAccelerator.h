// $Id: chebyshevAccelerator.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___CHEBYSHEV_ACCELERATOR
#define ___CHEBYSHEV_ACCELERATOR

#include "pijAccelerator.h"
#include "replacementModel.h"

class chebyshevAccelerator : public pijAccelerator {
public:

	explicit chebyshevAccelerator(	 replacementModel* pb,	
									const int alphanetSize=20,
									const int totalNumOfCoef=60,
									const int usingNumberOfCoef=13,
	const MDOUBLE rightRange=0,const MDOUBLE leftRange=2);

	// Constructor with precomputed coefficients for optimization
	chebyshevAccelerator(
		replacementModel* pb,
		const VVVdouble& precomputed_coff,
		const VVVdouble& precomputed_derv_coff,
		const VVVdouble& precomputed_sec_derv_coff,
		const int alphanetSize=20,
		const int totalNumOfCoef=60,
		const int usingNumberOfCoef=13,
		const MDOUBLE rightRange=0,
		const MDOUBLE leftRange=2
	);
  	chebyshevAccelerator(const chebyshevAccelerator& other);

	const MDOUBLE Qij(const int i,const int j) const {return _pb->Qij(i,j);}

	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE freq(const size_t i) const {return _pb->freq(i);}
	virtual pijAccelerator* clone() const { return new chebyshevAccelerator(*this); }
	virtual ~chebyshevAccelerator() {delete _pb;}
	virtual replacementModel* getReplacementModel() const {return (_pb);}
	virtual const size_t alphabetSize() const {return _pb->alphabetSize();}

	const VVVdouble& getChebiCoff() const { return chebi_coff; }
	const VVVdouble& getChebiDervCoff() const { return chebi_dervation_coff; }
	const VVVdouble& getChebiSecDervCoff() const { return chebi_sec_dervation_coff; }

private:
	VVVdouble chebi_coff;//[N_ABC][N_ABC][NUMBER_OF_TOTAL_COFF+1];
	VVVdouble chebi_dervation_coff;//[N_ABC][N_ABC][NUMBER_OF_TOTAL_COFF+1];
	VVVdouble chebi_sec_dervation_coff;//[N_ABC][N_ABC][NUMBER_OF_TOTAL_COFF+1];

	const int _alphabetSize;
	const int _totalNumOfCoef;
	const int _usingNumberOfCoef;
	
    replacementModel* _pb;

	void chebft(Vdouble& c, int n, int from_aa, int to_aa);
	void chder(Vdouble &c, Vdouble &cder, int n);

	const MDOUBLE _rightRange;
	const MDOUBLE _leftRange;

};

// This is an accelerator of Pij(t) calculation, using a proximity to polynomial.
#endif

