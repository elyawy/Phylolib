// Standalone program to generate precomputed Chebyshev coefficients
// Compile and run this once to generate PrecomputedChebyshevData.cpp

#include <iostream>
#include <fstream>
#include <iomanip>
#include "../includes/chebyshevAccelerator.h"
#include "../includes/datMatrixHolder.h"
#include "../includes/readDatMatrix.h"


void writeVVVdoubleToFile(const std::string& filename, const VVVdouble& data) {
    std::ofstream out(filename);
    
    // Write just the array data - no variable name, will be included directly
    out << "{\n";
    
    for (size_t i = 0; i < data.size(); ++i) {
        out << "    {\n";
        for (size_t j = 0; j < data[i].size(); ++j) {
            out << "        {";
            for (size_t k = 0; k < data[i][j].size(); ++k) {
                out << std::setprecision(17) << data[i][j][k];
                if (k < data[i][j].size() - 1) out << ", ";
            }
            out << "}";
            if (j < data[i].size() - 1) out << ",";
            out << "\n";
        }
        out << "    }";
        if (i < data.size() - 1) out << ",";
        out << "\n";
    }
    
    out << "}\n";
    out.close();
}

void generateModelData(const std::string& modelName, const datMatrixString& matrixData) {
    std::cout << "Generating " << modelName << " coefficients..." << std::endl;
    
    // Create the model
    pupAll model(matrixData);
    
    // Create chebyshev accelerator (this does the expensive computation)
    chebyshevAccelerator accel(&model);
    
    // Create directory if it doesn't exist (you may need to do this manually)
    std::string dir = "chebyshevData/";
    
    // Write three separate .dat files for this model
    writeVVVdoubleToFile(dir + modelName + "_coff.cheby.dat", accel.getChebiCoff());
    writeVVVdoubleToFile(dir + modelName + "_derv_coff.cheby.dat", accel.getChebiDervCoff());
    writeVVVdoubleToFile(dir + modelName + "_sec_derv_coff.cheby.dat", accel.getChebiSecDervCoff());
    
    std::cout << modelName << " - 3 .dat files created successfully!" << std::endl;
}

int main() {
    std::cout << "Generating precomputed Chebyshev coefficients for all models..." << std::endl;
    std::cout << "Creating chebyshevData/ directory (create manually if this fails)..." << std::endl;
    
    // Create output directory (may fail on some systems - create manually if needed)
    system("mkdir -p chebyshevData");
    
    // Generate all amino acid models
    generateModelData("cpREV45", datMatrixHolder::cpREV45);
    generateModelData("dayhoff", datMatrixHolder::dayhoff);
    generateModelData("jones", datMatrixHolder::jones);      // JTT
    generateModelData("mtREV24", datMatrixHolder::mtREV24);
    generateModelData("wag", datMatrixHolder::wag);
    generateModelData("HIVb", datMatrixHolder::HIVb);
    generateModelData("HIVw", datMatrixHolder::HIVw);
    generateModelData("lg", datMatrixHolder::lg);

    generateModelData("EX_BURIED", datMatrixHolder::EX_BURIED);
    generateModelData("EX_EXPOSED", datMatrixHolder::EX_EXPOSED);
    generateModelData("EHO_EXTENDED", datMatrixHolder::EHO_EXTENDED);
    generateModelData("EHO_HELIX", datMatrixHolder::EHO_HELIX);
    generateModelData("EHO_OTHER", datMatrixHolder::EHO_OTHER);
    generateModelData("EX_EHO_BUR_EXT", datMatrixHolder::EX_EHO_BUR_EXT);
    generateModelData("EX_EHO_BUR_HEL", datMatrixHolder::EX_EHO_BUR_HEL);
    generateModelData("EX_EHO_BUR_OTH", datMatrixHolder::EX_EHO_BUR_OTH);
    generateModelData("EX_EHO_EXP_EXT", datMatrixHolder::EX_EHO_EXP_EXT);
    generateModelData("EX_EHO_EXP_HEL", datMatrixHolder::EX_EHO_EXP_HEL);
    generateModelData("EX_EHO_EXP_OTH", datMatrixHolder::EX_EHO_EXP_OTH);
    
    std::cout << "\n===========================================\n";
    std::cout << "SUCCESS! All .dat files generated in chebyshevData/\n";
    std::cout << "Total: " << (20 * 3) << " .dat files (3 per model)\n";
    std::cout << "===========================================\n";
    
    return 0;
}