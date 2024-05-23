#include <iostream>
#include "../includes/tree.h"



int main() {

    tree t1("(A:0.5,B:0.5);", false);
    std::cout << (t1.getRoot())->getSon(0)->dis2father() << "\n";
    std::cout << t1.getNodesNum() << "\n";
    std::cout << t1.getAllBranchesLengthSum() << "\n";

    tree t2("(A:0.5);", false);
    std::cout << (t2.getRoot())->getSon(0)->dis2father() << "\n";
    std::cout << t2.getNodesNum() << "\n";
    std::cout << t2.getAllBranchesLengthSum() << "\n";

    tree t3("(((A:0.5):0.4,B:0.2):0.3);", false);
    std::cout << (t3.getRoot())->getSon(0)->dis2father() << "\n";
    std::cout << t3.getNodesNum() << "\n";
    std::cout << t3.getAllBranchesLengthSum() << "\n";


    return 0;
}