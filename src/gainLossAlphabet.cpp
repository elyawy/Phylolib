
#include "gainLossAlphabet.h"

gainLossAlphabet::gainLossAlphabet() {}

size_t gainLossAlphabet::fromChar(const char s) const{
	switch (s) {
		case '0': return 0; break;	
		case '1': return 1; break;
		//case '2': return 1; break;	// added to read seq with paralogs
		//case '3': return 1; break;	// added to read seq with paralogs
		//case '4': return 1; break;	// added to read seq with paralogs
		//case '5': return 1; break;	// added to read seq with paralogs
		//case '6': return 1; break;	// added to read seq with paralogs
		//case '7': return 1; break;	// added to read seq with paralogs
		//case '8': return 1; break;	// added to read seq with paralogs
		//case '9': return 1; break;	// added to read seq with paralogs
		case '-' : case'_' : return 3; break; // gap
		case '?' : case'*' : return 2; break; // unknown
		case 'x' : case'X' : return 2; break; // unknown
		default:
			vector<string> err;
			err.push_back(" The gainLoss sequences contained the character: ");
			err[0]+=s;
			err.push_back(" gainLoss was not one of the following: ");
			err.push_back(" 0, 1, or for unknown '?'/'-'");
			errorMsg::reportError(err);
	}// end of switch
	return -99; // never suppose to be here.	
}// end of function

vector<size_t> gainLossAlphabet::fromString(const string &str) const {
	vector<size_t> vec;
	for (unsigned int i=0;i<str.size();i++)
		vec.push_back(fromChar(str[i]));
	return vec;
}

string gainLossAlphabet::fromInt(const size_t in_id) const{
	char res = 0;
	switch (in_id) {
		case 0 : res = '0'  ; break;
		case 1 : res = '1'  ; break;
		case 3 : res = '-'; break;
		default:
			vector<string> err;
			err.push_back("unable to print gainLoss_id. gainLossl_id was not one of the following: ");
			err.push_back("0,1");
			errorMsg::reportError(err);
	}//end of switch
	string vRes;
	vRes.append(1,res);
	return vRes;
}// end of function

// There are no relations here.
size_t gainLossAlphabet::relations(const size_t charInSeq, const size_t charToCheck) const{
	if (charInSeq == charToCheck)
		return 1;
	if(charInSeq == -1 || charInSeq == -2) 
		return 1 ;//	missing data
	return 0;
}

size_t gainLossAlphabet::fromChar(const string& str, const size_t pos) const{
	return fromChar(str[pos]);
}




