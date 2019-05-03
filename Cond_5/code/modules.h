//  MABE is a product of The Hintza Lab @ MSU
//     for general research information:
//         http://hintzelab.msu.edu/
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//          github.com/Hintzelab/MABE/wiki

//  This file was auto-generated with MBuilder.py

#ifndef __AutoBuild__Modules__
#define __AutoBuild__Modules__
#include "World/TestWorld/TestWorld.h"
#include "World/EmptyWorld/EmptyWorld.h"
#include "Genome/CircularGenome/CircularGenome.h"
#include "Brain/ANNBrain/ANNBrain.h"
#include "Optimizer/TwoSexesOptimizer/TwoSexesOptimizer.h"

#include "Archivist/DefaultArchivist.h"
#include "Archivist/SSwDArchivist/SSwDArchivist.h"


//create a world
std::shared_ptr<AbstractWorld> makeWorld(std::shared_ptr<ParametersTable> PT){
  std::shared_ptr<AbstractWorld> newWorld;
  bool found = false;
  std::string worldType = AbstractWorld::worldTypePL->get(PT);
  if (worldType == "Test") {
    newWorld = std::make_shared<TestWorld>(PT);
    found = true;
    }
  if (worldType == "Empty") {
    newWorld = std::make_shared<EmptyWorld>(PT);
    found = true;
    }
  if (!found){
    std::cout << "  ERROR! could not find WORLD-worldType \"" << worldType << "\".\n  Exiting." << std::endl;
    exit(1);
    }
  return newWorld;
}


//create an optimizer
std::shared_ptr<AbstractOptimizer> makeOptimizer(std::shared_ptr<ParametersTable> PT){
  std::shared_ptr<AbstractOptimizer> newOptimizer;
  bool found = false;
  std::string optimizerType = AbstractOptimizer::Optimizer_MethodStrPL->get(PT);
  if (optimizerType == "TwoSexes") {
    newOptimizer = std::make_shared<TwoSexesOptimizer>(PT);
    found = true;
    }
  if (!found){
    std::cout << "  ERROR! could not find OPTIMIZER-optimizer \"" << optimizerType << "\".\n  Exiting." << std::endl;
    exit(1);
    }
  return newOptimizer;
}


//create an archivist
std::shared_ptr<DefaultArchivist> makeArchivist(std::vector<std::string> popFileColumns, std::shared_ptr<Abstract_MTree> _maxFormula, std::shared_ptr<ParametersTable> PT, std::string groupPrefix = ""){
  std::shared_ptr<DefaultArchivist> newArchivist;
  bool found = false;
  std::string archivistType = DefaultArchivist::Arch_outputMethodStrPL->get(PT);
  if (archivistType == "SSwD") {
    newArchivist = std::make_shared<SSwDArchivist>(popFileColumns, _maxFormula, PT, groupPrefix);
    found = true;
    }
  if (archivistType == "Default") {
    newArchivist = std::make_shared<DefaultArchivist>(popFileColumns, _maxFormula, PT, groupPrefix);
    found = true;
    }
  if (!found){
    std::cout << "  ERROR! could not find ARCHIVIST-outputMethod \"" << archivistType << "\".\n  Exiting." << std::endl;
    exit(1);
    }
  return newArchivist;
}


//create a template genome
std::shared_ptr<AbstractGenome> makeTemplateGenome(std::shared_ptr<ParametersTable> PT){
  std::shared_ptr<AbstractGenome> newGenome;
  bool found = false;
  std::string genomeType = AbstractGenome::genomeTypeStrPL->get(PT);
  if (genomeType == "Circular") {
    newGenome = CircularGenome_genomeFactory(PT);
    found = true;
    }
  if (found == false){
    std::cout << "  ERROR! could not find GENOME-genomeType \"" << genomeType << "\".\n  Exiting." << std::endl;
    exit(1);
    }
  return newGenome;
}


//create a template brain
std::shared_ptr<AbstractBrain> makeTemplateBrain(int inputs, int outputs, std::shared_ptr<ParametersTable> PT){
  std::shared_ptr<AbstractBrain> newBrain;
  bool found = false;
  std::string brainType = AbstractBrain::brainTypeStrPL->get(PT);
  if (brainType == "ANN") {
    newBrain = ANNBrain_brainFactory(inputs, outputs, PT);
    found = true;
    }
  if (found == false){
    std::cout << "  ERROR! could not find BRAIN-brainType \"" << brainType << "\".\n  Exiting." << std::endl;
    exit(1);
    }
  return newBrain;
}


//configure Defaults and Documentation
void configureDefaultsAndDocumentation(){
  Parameters::root->setParameter("BRAIN-brainType", (std::string)"ANN");
  Parameters::root->setDocumentation("BRAIN-brainType", "brain to be used, [ANN]");

  Parameters::root->setParameter("GENOME-genomeType", (std::string)"Circular");
  Parameters::root->setDocumentation("GENOME-genomeType", "genome to be used, [Circular]");

  Parameters::root->setParameter("ARCHIVIST-outputMethod", (std::string)"SSwD");
  Parameters::root->setDocumentation("ARCHIVIST-outputMethod", "output method, [SSwD, Default]");

  Parameters::root->setParameter("OPTIMIZER-optimizer", (std::string)"TwoSexes");
  Parameters::root->setDocumentation("OPTIMIZER-optimizer", "optimizer to be used, [TwoSexes]");

  Parameters::root->setParameter("WORLD-worldType", (std::string)"Test");
  Parameters::root->setDocumentation("WORLD-worldType","world to be used, [Test, Empty]");
}


#endif /* __AutoBuild__Modules__ */
