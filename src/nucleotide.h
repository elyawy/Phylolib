#ifndef ___NUCLEOTIDE_H
#define ___NUCLEOTIDE_H

#include "definitions.h" 
#include "alphabet.h"

#include <cassert>
using namespace std;

/* =======================================================================
This is the nucleotide class. It is derived from the class alphabet.
All alphabets are internally stored as integers. So what has to implement
is a way to translate from strings to array (vector) of integers and back.

Starting with the easiest functions to explain: 
size() gives the size of the alphabet: 4 in this case.
stringSize() say if it is a one letter code (unlike codon which is 3 letters code).

clone() is a general machanism in C++. The idea is that if you have a derived class,
and a pointer to the base class, and you want to self-copy the derived class.
In such case you use the clone() machanism. Ref: More effective C++ page. 126.

int unknown(): sometimes one doesn't know if it is A, C, G, or T. In such case we use
the int that represents unknown. In this class it is set to 15. This is used for example
when gap characters are converted to unknown characters.


int fromChar(const string& str, const int pos) and int fromChar(const char s)
give the same answer: there is a map from integers to characters.
For example, A is zero, C is 1, etc. However, the function fromChar(const char s)
is specific to nucleotide and to amino - because these are one letter alphabet.
For codon - this function won't work. This is why the general function is 
in the form int fromChar(const string& str, const int pos);
In the case of codon - it will read 3 letters each time.
=========================================================================*/



class nucleotide : public alphabet {
public:
	explicit nucleotide();
	virtual ~nucleotide() {}
    virtual alphabet* clone() const { return new nucleotide(*this); }
	size_t unknown() const  {return 15;}
	size_t gap() const  {return 16;}
	size_t size() const {return 4;}
	size_t stringSize() const {return 1;} // one letter code.

	size_t fromChar(const string& str, const size_t pos) const;
	size_t fromChar(const char s) const;
	vector<size_t> fromString(const string& str) const;

	string fromInt(const size_t id) const;
	size_t relations(const size_t charInSeq, const size_t charToCheck) const{ // see explanation below
		assert (charInSeq != gap());//gaps in the sequences
	  	return _relation[charToCheck][charInSeq];
	}
	
	// "specific" here is not unknown, nor ambiguity, nor gap (for example, for nucleotides it will true for A,C,G, or T).
	// in this speical case, in fact it will be true also for U which is coded by 4.
	// this is why it is <= size.
	bool isSpecific(const size_t id) const {return (id>=0 && id <= size());}

private:
	VVint _relation;
	char fromIntInternal(const int in_id)  const;
	int relationsInternal(const int ctc,const int charInSeq) const;

};

	
	
#endif



	

  



