#include <iostream>
#include <cassert>
#include "../includes/readDatMatrix.h"




int main() {

    hky model(0.2,0.1,0.4,0.3, 0.5);

    
    for (size_t i = 0; i < model.alphabetSize(); i++)
    {
        for (size_t j = 0; j < model.alphabetSize(); j++)
        {
            assert(model.dPij_dt(i,j,0) - model.Qij(i,j) < 0.00001);
            std::cout << model.dPij_dt(i,j,0) << " = " << model.Qij(i,j) << "\n";
        }
        
    }
    

}