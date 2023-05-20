// $Id: mulAlphabet.h 6420 2009-06-25 11:17:08Z adist $

// version 1.01
// last modified 1 Jan 2004

#ifndef ___MUL_ALPHABET_H
#define ___MUL_ALPHABET_H

#include "definitions.h"
#include "alphabet.h"
#include "someUtil.h"
//A.M. What is this thing? no record whatsoever! how can I understand from "mulAlphabet" what the purpose is???
//Sad, very sad
class mulAlphabet : public alphabet {

public:
	mulAlphabet(const alphabet* baseAlphabet, int mulFactor);
	mulAlphabet(const mulAlphabet& other);
	virtual ~mulAlphabet();
    virtual alphabet* clone() const { return new mulAlphabet(*this); }
	mulAlphabet& operator=(const mulAlphabet &other);

	size_t unknown() const ;
	size_t gap() const;
	
	size_t size() const {return _size;}
	size_t stringSize() const ;
	bool isSpecific(const size_t id) const ;

	size_t fromChar(const string& str, const size_t pos) const;
	vector<size_t> fromString(const string& str) const;

	string fromInt(const size_t id) const;

	size_t relations(const size_t charInSeq, const size_t charToCheck) const;
	int compareCategories(int charA, int charB) const;
	enum rateShiftType {noRateShift=0, acceleration, deceleration};
	static rateShiftType compareCategories(int charA, int charB, int baseAlphabetSize, int multiplicationFactor) ;
	const alphabet* getBaseAlphabet() const {return _baseAlphabet;}
	
public:
	int convertFromBasedAlphaInt(int id) const;
	int convertToBasedAlphaInt(int id) const;	
	
private:
	alphabet* _baseAlphabet; // This alphabet must use single characters, i.e. - not codon. (or we will have to add to every alphabet a member which holds its character's size)
	int _mulFactor ; // number of times that the alphabet is multiplied by = Number of categories (g in Galtier paper)
	size_t _size ; // this is simply the _baseAlphabet->size() * _mulFactor


};

#endif

