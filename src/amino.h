// $Id: amino.h 1901 2007-03-15 13:21:06Z nimrodru $

#ifndef ____AMINO
#define ____AMINO

#include "definitions.h"
#include "errorMsg.h"
#include "alphabet.h"
#include "geneticCodeHolder.h"
#include "codon.h"


//utility of amino acid
class aminoUtility {
public:

	static vector<int> codonOf(const int a, codon &cod); //returns vector of codons that code to a under a specific genetic code.
	
};

//based on the amino-acid list found in http://www.dur.ac.uk/~dbl0www/Bioinformatics/aminoacids.htm
class amino : public alphabet {
public:
	explicit amino();
	virtual ~amino() {}
	virtual alphabet* clone() const { return new amino(*this); }
	size_t unknown() const  {return 22;}
	size_t gap() const  {return 23;}
	size_t size() const {return 20;}
	size_t stringSize() const {return 1;} // one letter code.
	size_t relations(const size_t charInSeq, const size_t charToCheck) const;
	size_t fromChar(const string& str, const size_t pos) const;
	size_t fromChar(const char s) const;
	string fromInt(const size_t in_id) const;
	vector<size_t> fromString(const string& str) const;
	// "specific" here is not unknown, nor ambiguity, nor gap (for example, for nucleotides it will true for A,C,G, or T).
	bool isSpecific(const size_t id) const {return (id>=0 && id < size());}

private:
	int relations_internal(const int charInSeq, const int charToCheck) const;
	VVint _relation;
};//end of class

#endif


