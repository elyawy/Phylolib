// $Id: recognizeFormat.cpp 6780 2009-09-15 00:55:05Z itaymay $


#include "../includes/recognizeFormat.h"
#include "../includes/maseFormat.h"
#include "../includes/sequenceContainer.h"
#include "../includes/molphyFormat.h"
#include "../includes/phylipFormat.h"
#include "../includes/nexusFormat.h"
#include "../includes/fastaFormat.h"
#include "../includes/clustalFormat.h"
#include "../includes/nexusFormat.h"
#include "../includes/phylipSequentialFormat.h"


sequenceContainer recognizeFormat::read(istream &infile, const alphabet* alph) {
	sequenceContainer mySeqData = readUnAligned(infile, alph);
	mySeqData.makeSureAllSeqAreSameLengthAndGetLen();
	return mySeqData;
}

sequenceContainer recognizeFormat::readUnAligned(istream &infile, const alphabet* alph) {
	// recognize a format and returns the sequence container of it.
	sequenceContainer sc;
	if (!infile){
		string tmp = "error unable to open sequence input file ";
		errorMsg::reportError(tmp);
	}

	// this part eats spaces, tabs and such.
	char check = infile.peek();
	while ((check==' ') || (check == '\n') || (check == '\t')) {
		infile.get();
		check = infile.peek();
	}

	switch (check){
	case '#':
		sc=nexusFormat::readUnAligned(infile,alph);
		break;
	case '>':
		sc=fastaFormat::readUnAligned(infile,alph);
		break;
	case 'C':
		sc=clustalFormat::readUnAligned(infile,alph);
		break;
	case ';':
		sc=maseFormat::readUnAligned(infile,alph);
		break;	

	default:
		if (isdigit(check)){ 
			// here it can be either MOLPHY format or one of the PHYLIP type formats (interleaved, sequential)
			// in PHYLIP format there are lines that are not empty, but the first 10 characters
			// are space.
			string s;
			getline(infile,s, '\n' ); // read the first line which are numbers in both formats
			getline(infile,s, '\n' ); // read the second line
			bool phylipFormat = false;
			int r = s.find_first_of(' '); // if there is a space somewhere - this is phylip format
			if ((r==(s.size()-1)) || (r==-1)) phylipFormat = false;
			else phylipFormat = true;
			
			
			if (phylipFormat == false) {
				infile.seekg(0, ios::beg); // file return to the beginning
				sc=molphyFormat::readUnAligned(infile,alph);
			} else {
				getline(infile,s, '\n' ); // read the third line: interleaved will begin with a space, sequential not
				infile.seekg(0, ios::beg); // file return to the beginning
				if (s[0] == ' ') 
					sc = phylipSequentialFormat::readUnAligned(infile, alph);
				else 
					sc = phylipFormat::readUnAligned(infile,alph);
			}
		}
		else{
			string line;
			getline(infile, line, '\n');
			string tmp2 = "The program can't recognise your format!";
			tmp2+="\nThis is the first line in your format:\n";
			tmp2+=line;
			errorMsg::reportError(tmp2);
		}
		break;
	}
	return sc;
}
