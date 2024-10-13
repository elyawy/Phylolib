#include "../includes/sequence.h"
#include "../includes/nucleotide.h"

int main() {

    nucleotide a;
    std::string strSeq = "AAATTCCATG";

    sequence seq(strSeq, "Seq1","",1,&a);

    


    // seq.resize(100);

    // seq[50] = 0;

    std::cout << seq << "\n";


    return 0;
}