// $Id: amino.cpp 2414 2007-10-08 14:34:42Z adist $

#include "../includes/amino.h"

//VVint amino::_relation;

amino::amino() {
	_relation.resize(24); 	// relation should realy be an allocted, two dimentional array, not a vector. 
	for (size_t i=0; i < _relation.size(); ++i) { // this implementation would be much faster.  with some c-tricks, this checkup could be done with one access only.
		_relation[i].resize(20);
	}

	for (size_t k=0;k<24;++k){
		for (size_t j=0;j<20;++j){
			_relation[k][j]=relations_internal(k,j);
		}
	}
}

ALPHACHAR amino::fromChar(const char s) const{
	switch (s) {
	case 'A' : case'a' : return 0 ; break;
	case 'R' : case'r' : return 1 ; break;
	case 'N' : case'n' : return 2 ; break;
	case 'D' : case'd' : return 3 ; break;
	case 'C' : case'c' : return 4 ; break;
	case 'Q' : case'q' : return 5 ; break;
	case 'E' : case'e' : return 6 ; break;
	case 'G' : case'g' : return 7 ; break;
	case 'H' : case'h' : return 8 ; break;
	case 'I' : case'i' : return 9 ; break;
	case 'L' : case'l' : return 10; break;
	case 'K' : case'k' : return 11; break;
	case 'M' : case'm' : return 12; break;
	case 'F' : case'f' : return 13; break;
	case 'P' : case'p' : return 14; break;
	case 'S' : case's' : return 15; break;
	case 'T' : case't' : return 16; break;
	case 'W' : case'w' : return 17; break;
	case 'Y' : case'y' : return 18; break;
	case 'V' : case'v' : return 19; break;
	case 'B' : case'b' : return 20 ; break; // aspartate(D) or asparagine(N)
	case 'Z' : case'z' : return 21 ; break; // glutamate (E) or glutamine(Q)
	case '-' : case'_' : return 23; break;
	case '?' : case'*' : return 22; break;
	case 'x' : case'X' : return 23; break;
	case '.' : return 24; break;
	default:
	  vector<string> err;
	  err.push_back(" The amino-acid sequences contained the character: ");
	  err[0]+=s;
	  err.push_back(" Amino acid was not one of the following: ");
	  err.push_back(" A, B, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, X, Z, -, ?");
	  err.push_back(" a, b, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v, x, z, _, *");
	  errorMsg::reportError(err);
	}// end of switch
	return -99; // never suppose to be here.	
}// end of function

vector<ALPHACHAR> amino::fromString(const string &str) const {
	vector<ALPHACHAR> vec;
	for (size_t i=0;i<str.size();i++)
	  vec.push_back(fromChar(str[i]));
	return vec;
}

string amino::fromInt(const ALPHACHAR in_id) const{
  char res = 0;
	switch (in_id) {
		case 0 : res = 'A'  ; break;
		case 1 : res = 'R'  ; break;
		case 2 : res = 'N'  ; break;
		case 3 : res = 'D'  ; break;
		case 4 : res = 'C'  ; break;
		case 5 : res = 'Q'  ; break;
		case 6 : res = 'E'  ; break;
		case 7 : res = 'G'  ; break;
		case 8 : res = 'H'  ; break;
		case 9 : res = 'I'  ; break;
		case 10: res = 'L'  ; break;
		case 11: res = 'K'  ; break;
		case 12: res = 'M'  ; break;
		case 13: res = 'F'  ; break;
		case 14: res = 'P'  ; break;
		case 15: res = 'S'  ; break;
		case 16: res = 'T'  ; break;
		case 17: res = 'W'  ; break;
		case 18: res = 'Y'  ; break;
		case 19: res = 'V'  ; break;
		case 20: res = 'B'  ; break;
		case 21: res = 'Z'  ; break;
		case 22: res = '?'  ; break;
		case 23: res = '-'  ; break;
		
		default:
		vector<string> err;
		err.push_back(" unable to print amino ac_id. amino ac_id was not one of the following: ");
		err.push_back("A, B, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, Z, -, ?");
		err.push_back("a, b, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v, z, _, *");
		errorMsg::reportError(err);
		}//end of switch
	string vRes;
	vRes.append(1,res);
	return vRes;
}// end of function

ALPHACHAR amino::relations(const ALPHACHAR charInSeq, const ALPHACHAR charToCheck) const{
	if (charInSeq == gap()) {
		errorMsg::reportError("gaps in the sequences. Either change gaps to ? or remove gap positions");
	}
	return _relation[charInSeq][charToCheck];
}

ALPHACHAR amino::fromChar(const string& str, const size_t pos) const{
	return fromChar(str[pos]);
}

ALPHACHAR amino::relations_internal(const ALPHACHAR charInSeq, const ALPHACHAR charToCheck) const{
	if (charInSeq == charToCheck) return 1;
	else if (charInSeq == fromChar('?')) return 1;
	else if ((charInSeq == fromChar('B')) && 
		 ((charToCheck == fromChar('N')) || 
		  (charToCheck == fromChar('D')))) return 1; // B is either N or D
	else if ((charInSeq == fromChar('Z')) && 
		 ((charToCheck == fromChar('Q')) ||
		  (charToCheck == fromChar('E')))) return 1; // Z is either E or Q
	return 0;
}


vector<ALPHACHAR> aminoUtility::codonOf(const ALPHACHAR a, codon &cod){
	vector<ALPHACHAR> codons;
	amino amin;
	string strAmino=amin.fromInt(a);
	map <string, string> genCode=cod.geneticCode();
	map <string, string>::iterator it=genCode.begin();
	int tmp2=genCode.size();
	while (it!=genCode.end()){
		string tmp=(*it).second;
		if ((*it).second==strAmino){
			string strCodon=(*it).first;
			int c=cod.fromChar(strCodon,0);
			codons.push_back(c);		
		}
		it++;
	}
	if (codons.empty()){
		cout<<tmp2<<" amino is  = "<<a<<endl;
		errorMsg::reportError("error in function aminoUtility::codonOf: no codon found for amino acid");
	}
	return codons;
}
