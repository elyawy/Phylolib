// $Id: phylipFormat.h 1812 2007-03-01 09:29:12Z adist $

#ifndef ___PHYLIP_FORMAT
#define ___PHYLIP_FORMAT

#include "definitions.h"
#include "sequenceContainer.h"
#include <ostream>
using namespace std;

class phylipFormat {
public:
	static sequenceContainer read(istream &infile, const alphabet* alph);
	static void write(ostream &out, const sequenceContainer& sd,
		const size_t numOfPositionInLine = 50,
		const size_t spaceEvery = 10);
		
	//readUnAligned: the input sequences do not need to be aligned (not all sequences are the same length).
	static sequenceContainer readUnAligned(istream &infile, const alphabet* alph);
};

#endif

/* EXAMPLE OF PHYLIP FORMAT (interleaved):

6   128
Langur     KIFERCELAR TLKKLGLDGY KGVSLANWVC LAKWESGYNT EATNYNPGDE
Baboon     KIFERCELAR TLKRLGLDGY RGISLANWVC LAKWESDYNT QATNYNPGDQ
Human      KVFERCELAR TLKRLGMDGY RGISLANWMC LAKWESGYNT RATNYNAGDR
Rat        KTYERCEFAR TLKRNGMSGY YGVSLADWVC LAQHESNYNT QARNYDPGDQ
Cow        KVFERCELAR TLKKLGLDGY KGVSLANWLC LTKWESSYNT KATNYNPSSE
Horse      KVFSKCELAH KLKAQEMDGF GGYSLANWVC MAEYESNFNT RAFNGKNANG

           STDYGIFQIN SRYWCNNGKP GAVDACHISC SALLQNNIAD AVACAKRVVS
           STDYGIFQIN SHYWCNDGKP GAVNACHISC NALLQDNITD AVACAKRVVS
           STDYGIFQIN SRYWCNDGKP GAVNACHLSC SALLQDNIAD AVACAKRVVR
           STDYGIFQIN SRYWCNDGKP RAKNACGIPC SALLQDDITQ AIQCAKRVVR
           STDYGIFQIN SKWWCNDGKP NAVDGCHVSC SELMENDIAK AVACAKKIVS
           SSDYGLFQLN NKWWCKDNKR SSSNACNIMC SKLLDENIDD DISCAKRVVR

           DQGIRAWVAW RNHCQNKDVS QYVKGCGV
           DQGIRAWVAW RNHCQNRDVS QYVQGCGV
           DQGIRAWVAW RNRCQNRDVR QYVQGCGV
           DQGIRAWVAW QRHCKNRDLS GYIRNCGV
           EQGITAWVAW KSHCRDHDVS SYVEGCTL
           DKGMSAWKAW VKHCKDKDLS EYLASCNL


*/

