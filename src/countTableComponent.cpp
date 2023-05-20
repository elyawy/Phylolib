// $Id: countTableComponent.cpp 9595 2011-06-30 18:56:40Z rubi $

// version 1.00
// last modified 3 Nov 2002

#include "countTableComponent.h"
#include "logFile.h"

void countTableComponentHom::zero() {
	for (size_t alphabetChar1=0; alphabetChar1 < _countValues.size() ;++alphabetChar1) {
		for (size_t alphabetChar2=0; alphabetChar2 < _countValues[alphabetChar1].size() ;++alphabetChar2) {
			_countValues[alphabetChar1][alphabetChar2] = 0;
		}
	}
}

void countTableComponentHom::countTableComponentAllocatePlace(
		const int alphabetSize) {
	_countValues.resize(alphabetSize);
	for (int alphabetChar=0; alphabetChar < alphabetSize;++alphabetChar) _countValues[alphabetChar].resize(alphabetSize);
}

void countTableComponentHom::printTable(ostream& out) const {
	MDOUBLE sumCheck = 0.0;
	for (size_t i=0; i < _countValues.size();++i) {
		for (size_t k=0; k <  _countValues.size();++k) {
			out<<"counts["<<i<<"]["<<k<<"]"<<_countValues[i][k];
			sumCheck += _countValues[i][k];
			out<<endl;
		}
	}
	out<<"sum is: "<<sumCheck<<endl;
}

