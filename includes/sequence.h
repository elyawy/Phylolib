#ifndef ___SEQUENCE
#define ___SEQUENCE
#include "definitions.h"
#include "errorMsg.h"
#include "alphabet.h"
#include "mulAlphabet.h"
#include <iostream>
using namespace std;

class sequence {


public:
	class Iterator;
	friend class Iterator;
	class constIterator;
	friend class constIterator;

	// constructors
	explicit sequence(const string& str,
					const string& name,
					const string& remark,
					const size_t id,
					const alphabet* inAlph);

	sequence(const sequence& other);
	sequence(const sequence& other,const alphabet* inAlph); // convert the other sequence to the alphabet inAlph. 
	explicit sequence(const alphabet* inAlph) {
		if (inAlph == NULL) {
			errorMsg::reportError("must give a non Null alphabet when constructing sequences");
		}
		_alphabet = inAlph->clone();
		if (_alphabet->size() == 4) _lostBits = 0;
		if (_alphabet->size() == 20) _lostBits = 4;
		if (_alphabet->size() == 61) _lostBits = 6;
		_length = 0;
	}
	virtual ~sequence();

	void push_back(const ALPHACHAR &charId);

	size_t seqLen() const {return _vec.size();}
	size_t seqLenSpecific() const; //return the number of sites that are isSpecific()
	const string& name() const {return _name;}
	void setName(const string & inName)  { _name =inName ;}
	const size_t id() const {return _id;}
	void setID(const int inID)  { _id =inID ;}
	const string& remark() const {return _remark;}
	void setRemarks(const string & inRemarks)  { _remark =inRemarks ;}
	string toString() const;
	string toString(const int pos) const;

	void addFromString(const string& str);
	//push_back: add a single characer to the sequence
	void push_back(ALPHACHAR p) {_vec.push_back(p);}
	void resize(const int k, const int* val = NULL);
	void removePositions(const vector<int> & parCol);

	void setAlphabet(const alphabet* inA) {if (_alphabet) delete _alphabet;
		_alphabet=inA->clone();
	} 
	const alphabet* getAlphabet() const {return _alphabet;}

	inline sequence& operator=(const sequence& other);
	inline sequence& operator+=(const sequence& other);
	ALPHACHAR& operator[](const int i) {
		_vec[i]

		return _vec[i];
	}
	const ALPHACHAR& operator[](const int pos) const {return _vec[pos];}

	bool isUnknown(const int pos) const {return _vec[pos] == _alphabet->unknown();}
	
	// "specific" here is not unknown, nor ambiguity, nor gap (for example, for nucleotides it will true for A,C,G, or T).
	bool isSpecific(const size_t pos) const {return _alphabet->isSpecific(_vec[pos]);}

private: 
	vector<uint64_t> _vec;
	size_t _length;
	size_t _charsPerInt64;
	size_t _bitOffset;

	size_t _lostBits;
	const alphabet* _alphabet;
	string _remark;
	string _name;
	size_t _id;


public:
	class Iterator {
	public:
		explicit Iterator(){};
		~Iterator(){};
		void begin(sequence& seq){_pointer = seq._vec.begin();}
		void end(sequence& seq){_pointer = seq._vec.end();}
		ALPHACHAR& operator* (){return *_pointer;}
		ALPHACHAR const &operator* () const {return *_pointer;}
		void operator ++() {++_pointer;}
		void operator --() { --_pointer; }
		bool operator != (const Iterator& rhs){return (_pointer != rhs._pointer);}
		bool operator == (const Iterator& rhs){return (_pointer == rhs._pointer);}
	private:
		vector<uint64_t>::iterator _pointer;
  };

	class constIterator {
	public:
		explicit constIterator(){};
		~constIterator(){};
		void begin(const sequence& seq){_pointer = seq._vec.begin();}
		void end(const sequence& seq){_pointer = seq._vec.end();}
		ALPHACHAR const &operator* () const {return *_pointer;}
		void operator ++(){++_pointer;}
		void operator --(){--_pointer;}
		bool operator != (const constIterator& rhs) {
		  return (_pointer != rhs._pointer);
		}
		bool operator == (const constIterator& rhs) {
		  return (_pointer == rhs._pointer);
		}
	private:
		vector<ALPHACHAR>::const_iterator _pointer;
	};


} ;

inline sequence& sequence::operator=(const sequence& other) {
	_vec = other._vec;
	_alphabet = other._alphabet->clone();
	_name=other.name();
	_id=other.id();
	_remark=other.remark();

	return *this;
}

inline sequence& sequence::operator+=(const sequence& other) {
	for (size_t i=0; i <other._vec.size();++i) { 
		_vec.push_back(other._vec[i]);
	}
	return *this;
}


inline ostream & operator<<(ostream & out, const sequence &Seq){
    out<< Seq.toString();
    return out;
}


#endif