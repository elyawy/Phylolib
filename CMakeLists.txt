cmake_minimum_required(VERSION 3.15)

project(Phylolib VERSION 24.2)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_CXX_STANDARD 14)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
endif(CMAKE_COMPILER_IS_GNUCXX)


set(CMAKE_CXX_FLAGS_RELEASE "-O2")
############################################################
# Create a library
############################################################

#Generate the static library from the library sources
add_library(Phylolib STATIC 
    src/aaJC.cpp
    src/alphabet.cpp
    src/amino.cpp
    src/bblEM.cpp
    src/bblEM2codon.cpp
    src/bblEMProportionalEB.cpp
    src/bblLSProportionalEB.cpp
    src/bestAlpha.cpp
    src/bestAlphaAndK.cpp
    src/bestGtrModelParams.cpp
    src/bestHKYparam.cpp
    src/bestTamura92param.cpp
    src/betaDistribution.cpp
    src/betaOmegaDistribution.cpp
    src/betaUtilities.cpp
    src/chebyshevAccelerator.cpp
    src/clustalFormat.cpp
    src/codon.cpp
    src/codonUtils.cpp
    src/computeCounts.cpp
    src/computeDownAlg.cpp
    src/computeMarginalAlg.cpp
    src/computePijComponent.cpp
    src/computeUpAlg.cpp
    src/ConversionUtils.cpp
    src/countTableComponent.cpp
    src/customDistribution.cpp
    src/datMatrixHolder.cpp
    src/distanceTable.cpp
    src/distribution.cpp
    src/doubleRep.cpp
    src/errorMsg.cpp
    src/evaluateCharacterFreq.cpp
    src/fastaFormat.cpp
    src/fastStartTree.cpp
    src/fromCountTableComponentToDistance.cpp
    src/fromCountTableComponentToDistance2Codon.cpp
    src/fromCountTableComponentToDistancePropEB.cpp
    src/fromQtoPt.cpp
    src/gainLossAlphabet.cpp
    src/gammaDistribution.cpp
    src/gammaUtilities.cpp
    src/generalGammaDistribution.cpp
    src/geneticCodeHolder.cpp
    src/getRandomWeights.cpp
    src/gtrModel.cpp
    src/hky.cpp
    src/likeDist.cpp
    src/likeDist2Codon.cpp
    src/likeDistPropEB.cpp
    src/likelihoodComputation.cpp
    src/likelihoodComputation2Codon.cpp
    src/likelihoodComputationFactors.cpp
    src/likelihoodComputationGL.cpp
    src/logFile.cpp
    src/maseFormat.cpp
    src/matrixUtils.cpp
    src/molphyFormat.cpp
    src/mulAlphabet.cpp
    src/multipleStochasticProcess.cpp
    src/nexusFormat.cpp
    src/nj.cpp
    src/njConstrain.cpp
    src/nucJC.cpp
    src/nucleotide.cpp
    src/numRec.cpp
    src/Parameters.cpp
    src/phylipFormat.cpp
    src/phylipSequentialFormat.cpp
    src/pijAccelerator.cpp
    src/readDatMatrix.cpp
    src/readTree.cpp
    src/recognizeFormat.cpp
    src/replacementModel.cpp
    src/seqContainerTreeMap.cpp
    src/sequence.cpp
    src/sequenceContainer.cpp
    src/siteSpecificRate.cpp
    src/someUtil.cpp
    src/stochasticProcess.cpp
    src/suffStatComponent.cpp
    src/talRandom.cpp
    src/tamura92.cpp
    src/tree.cpp
    src/treeIt.cpp
    src/treeUtil.cpp
    src/uniDistribution.cpp
    src/uniformDistribution.cpp
    src/unObservableData.cpp
    src/wYangModel.cpp
)

target_include_directories(Phylolib
    PUBLIC 
        ${PROJECT_SOURCE_DIR}/includes
)

set_target_properties(Phylolib PROPERTIES PUBLIC_HEADER "includes/aaJC.h;includes/allModels.h;includes/alphabet.h;includes/amino.h;includes/bblEM.h;includes/bblEM2codon.h;includes/bblEMProportionalEB.h;includes/bblLSProportionalEB.h;includes/bestAlpha.h;includes/bestAlphaAndK.h;includes/bestGtrModelParams.h;includes/bestHKYparam.h;includes/bestTamura92param.h;includes/betaDistribution.h;includes/betaOmegaDistribution.h;includes/betaUtilities.h;includes/chebyshevAccelerator.h;includes/checkcovFanctors.h;includes/checkcovFanctorsWithFactors.h;includes/clustalFormat.h;includes/codon.h;includes/codonUtils.h;includes/computeCounts.h;includes/computeDownAlg.h;includes/computeMarginalAlg.h;includes/computePijComponent.h;includes/computeUpAlg.h;includes/ConversionUtils.h;includes/countTableComponent.h;includes/src/customDistribution.h;includes/datMatrixHolder.h;includes/definitions.h;includes/distanceMethod.h;includes/distances2Tree.h;includes/distanceTable.h;includes/distribution.h;includes/doubleRep.h;includes/errorMsg.h;includes/evaluateCharacterFreq.h;includes/fastaFormat.h;includes/fastStartTree.h;includes/fromCountTableComponentToDistance.h;includes/fromCountTableComponentToDistance2Codon.h;includes/fromCountTableComponentToDistancePropEB.h;includes/fromQtoPt.h;includes/gainLossAlphabet.h;includes/gammaDistribution.h;includes/gammaUtilities.h;includes/generalGammaDistribution.h;includes/geneticCodeHolder.h;includes/getopt.h;includes/getRandomWeights.h;includes/gtrModel.h;includes/hky.h;includes/jcDistance.h;includes/likeDist.h;includes/likeDist2Codon.h;includes/likeDistPropEB.h;includes/likelihoodComputation.h;includes/likelihoodComputation2Codon.h;includes/likelihoodComputationFactors.h;includes/likelihoodComputationGL.h;includes/logFile.h;includes/maseFormat.h;includes/matrixUtils.h;includes/modelFactory.h;includes/molphyFormat.h;includes/mulAlphabet.h;includes/multipleStochasticProcess.h;includes/nexusFormat.h;includes/nj.h;includes/njConstrain.h;includes/nucJC.h;includes/nucleotide.h;includes/numRec.h;includes/Parameters.h;includes/phylipFormat.h;includes/phylipSequentialFormat.h;includes/pijAccelerator.h;includes/rateMatrixSim.h;includes/readDatMatrix.h;includes/readTree.h;includes/recognizeFormat.h;includes/replacementModel.h;includes/seqContainerTreeMap.h;includes/sequence.h;includes/sequenceContainer.h;includes/siteSpecificRate.h;includes/someUtil.h;includes/stochasticProcess.h;includes/suffStatComponent.h;includes/talRandom.h;includes/tamura92.h;includes/targetver.h;includes/tree.h;includes/treeIt.h;includes/treeUtil.h;includes/trivialAccelerator.h;includes/uniDistribution.h;includes/uniformDistribution.h;includes/unObservableData.h;includes/wYangModel.h")


install (TARGETS Phylolib
    ARCHIVE DESTINATION build
    PUBLIC_HEADER DESTINATION build/includes)

include(CPack)