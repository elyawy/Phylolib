#include <iostream>
#include <cassert>
#include <iomanip>
#include "../includes/readDatMatrix.h"
#include "../includes/gtrModel.h"



int main() {
    auto model = pupAll(datMatrixHolder::wag, 20); 
    std::cout << model.Qij(0,0) << "\n";


    double diagonal_sum = 0.0;
    
    for (size_t i = 0; i < model.alphabetSize(); i++)
    {
        diagonal_sum += (model.Qij(i,i)*model.freq(i));
        MDOUBLE sum = 0;
        for (size_t j = 0; j < model.alphabetSize(); j++)
        {
            sum += model.Qij(i,j);
            // std::cout  << model.Qij(i,j) << "  ";

            assert(model.Qij(i,j)*model.freq(i) - model.Qij(j,i)*model.freq(j) < 0.00001);
            std::cout << model.Qij(i,j)*model.freq(i) << " = " << model.Qij(j,i)*model.freq(j) << "\n";
        }
        std::cout << "\n" << sum << "\n";
        
    }
    std::cout << diagonal_sum << "\n";


}
