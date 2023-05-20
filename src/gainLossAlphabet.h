#ifndef ___GAIN_LOSS_ALPH
#define ___GAIN_LOSS_ALPH

#include "alphabet.h"
#include "errorMsg.h"

class gainLossAlphabet : public alphabet {
public:
	explicit gainLossAlphabet(); 
	virtual ~gainLossAlphabet() {}
	virtual alphabet* clone() const { return new gainLossAlphabet(*this); }
	size_t unknown() const  {return 2;}
	size_t gap() const  {return 3;} // 
	size_t size() const {return 2;} // presence or absence only
	size_t stringSize() const {return 1;} // one letter code.
	size_t relations(const size_t charInSeq, const size_t charToCheck) const;
	size_t fromChar(const string& str, const size_t pos) const;
	size_t fromChar(const char s) const;
	string fromInt(const size_t in_id) const;
	vector<size_t> fromString(const string& str) const;
	bool isSpecific(const size_t id) const {return (id>=0 && id < static_cast<int>(size()));}

};

#endif
