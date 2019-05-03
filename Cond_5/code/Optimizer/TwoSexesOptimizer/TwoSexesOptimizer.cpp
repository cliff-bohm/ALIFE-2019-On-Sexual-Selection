//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#include "TwoSexesOptimizer.h"

//using namespace std;




std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::additiveNamesPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_GEN_DISPLAYS-additiveNames", (std::string)"", "list of names for additive markers. example: marker_A,marker_B");
std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::additiveSourceGenomesPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_GEN_DISPLAYS-additiveSourceGenomes", (std::string)"", "genome (nameSpace) which should be seached for each additive marker. "
	"If organisms do not already have a named genome, it will be added. example: root::,root::");
std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::additiveStartCodonsPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_GEN_DISPLAYS-additiveStartCodons", (std::string)"", "list of start codons for additive markers."
	"\nstartcodons are lists of 1 and 8 integer values between [0 and 255] seperated by ':'. example 11:12,15:16");
std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::assignedNamesPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_GEN_DISPLAYS-assignedNames", (std::string)"Marker_Sex,Marker_Preference", "list of names for assigned markers. example Marker_Sex,Marker_Preference");
std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::assignedSourceGenomesPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_GEN_DISPLAYS-assignedSourceGenomes", (std::string)"root::,root::", "genome (nameSpace) which should be seached for each assigned marker. "
	"If organisms do not already have a named genome, it will be added. example: genome_A::,root::");
std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::assignedStartCodonsPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_GEN_DISPLAYS-assignedStartCodons", (std::string)"21:34,22:33", "list of start codons for assigned markers"
	"\nstartcodons are lists of 1 and 8 integer values between [0 and 255] seperated by ':'. example: 21:34,22:33");
std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::assignedRangesPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_GEN_DISPLAYS-assignedRanges", (std::string)"1,1", "list of ranges (ints) for assigned markers, genetic values will be between 0 and the value"
	"\nexample [1,3] would mean that the first assigned trait could be 0 or 1 and the second could be 0,1,2 or, 3");
std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::assignedDuplicatesRulesPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_GEN_DISPLAYS-assignedDuplicatesRules", (std::string)"Unique,Unique",
	"list of rules for assigned markers defining behavior in case of duplicates. example: Most,Unique"
	"\nUnique - only allow one codon in genome (more then one will return -1 i.e. none found)"
	"\nMost - locate all condons and return the highest frequency (ties decided by last in genome)"
	"\nAverage - locate all condons and return the average of all values found (rounded down)"
	"\nGreatest - greatest value found"
	"\nLeast - least value found");
std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::seedGenomesPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_GEN_DISPLAYS-seedGenomes", (std::string)"root::,root::", "list of names of genomes to seed with codons");
std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::seedCodonsPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_GEN_DISPLAYS-seedCodons", (std::string)"100<21:34,1x22:33", "list of codons to seed. seedGenomes will determine which genome will be seeded by each. format: [copies]x[start codon list],[location]@[start codon list]\n"
	"If x is used, then the number before the start codon list will be the number of copies, placed randomly\n"
	"If < is used, then the number before the start codon list will be the location were one copy of the start codon will be placed.");
std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::seedRandomizeGenomesPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_GEN_DISPLAYS-seedRandomizeGenomes", (std::string)"", "genomes in this list will be randomized on update 0 before seeding codons");

std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::sexFormulaPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES-sexFormula", 
	(std::string) "DM_AVE[Marker_Sex]", "MTree used to determine sex of organism");

std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::preferenceFormulaPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES-preferenceFormula",
	(std::string) "DM_AVE[Marker_Preference]", "MTree used to determine sexual preference of organism");

std::shared_ptr<ParameterLink<int>> TwoSexesOptimizer::tournamentSizeMomPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_PARENTS_SELECTION-tournamentSizeMom", 5, "if using Tournament Method, how many orgs are considered when selecting a mom?");
std::shared_ptr<ParameterLink<int>> TwoSexesOptimizer::tournamentSizeDadPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_PARENTS_SELECTION-tournamentSizeDad", 5, "if using Tournament Method, how many orgs are considered for each dad (i.e. selection group member)?");

std::shared_ptr<ParameterLink<int>> TwoSexesOptimizer::selectionGroupSizePL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES-selectionGroupSize", 5, "size of lek if selectionSystem is Lek, number of males female will see if selectionSystem is Threshold"
	"\ni.e. the number of dads from which each mom will need to choose."
	"\nNOTE: this should not be confused with the tournamentSizeMom and tournamentSizeDad which are used when selecting moms and dads (for a selection group)");
std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::selectionSystemsPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES-selectionSystems", (std::string)"Lek", "List of methods used with each sex trait (if only one system is listed, then all traits will use that method)"
	"\noptions are Lek or Threshold");

std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::sexTraitsPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES-sexTraits", (std::string) "0,DM_AVE[score]%0.0%0.0", "list of sexual traits"
	"\na constant value trait will result in all organisms having the same value for the trait - this value can be effected by optimizeValueDadMateEffect/optimizeValueMomMateEffect."
	"\n%p%p may optionally follow a preference:"
	"\nfirst %p sets detection error; probability that selection will be proportional (in range of values in lek) as apposed to perfectly selecting for best"
	"\nsecond %p is probability that selection will simply fail (i.e. female selects random)"
	"\ntraits are defined with MTree with access to organisms DM, and this optimizers PT");

std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::sexTraitsThresholdsPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES-sexTraitsThresholds",
	(std::string) "", "list of sexual traits thresholds, to be used if selectionSystem is Threshold, for a given trait.(MTree)"
	"\nIf a mix of Threshold and Lex methods are being used, the Lek methods must be assigned a threshold as a place holder (which will be ignored, and can be 0)");

std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::optimizeValueMomPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_PARENTS_SELECTION-optimizeValueMom", (std::string) "DM_AVE[score]", "value used to pick moms (MTree)");
std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::optimizeValueDadPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_PARENTS_SELECTION-optimizeValueDad", (std::string) "DM_AVE[score]", "value to choose dads when picking members for selection group (MTree)");

std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::traitConditionEffectsPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES-traitConditionEffects", (std::string) "[0]", "how is each trait correlated with reproductve success? only applies to males."
	"\nif [0], then there are no condtion effects."
	"\nif there is more then one value in this list, then there must be one value for each sexual trait."
	"\neach value in list will, relative to the male population scores, additivly determin the likelihood"
	"\nthat a selected male will produce a viable offspring (this is experimental, biologically questionable!)");

std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::selectionMethodPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_PARENTS_SELECTION-selectionMethod", (std::string) "Roulette", "Tournament (will use tournament sizes) or Roulette will be used to select moms and dads for selection groups");

std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::surviveRatePL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_PARENTS_SELECTION-surviveRate", (std::string) "0", "value between 0 and 1, likelyhood that an organism survive (MTree)");

std::shared_ptr<ParameterLink<bool>> TwoSexesOptimizer::elitismPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_PARENTS_SELECTION-elitism", false, "if true, dad with best dad score and mom with best mom score will survive");

std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::reMapFormulaMomPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_PARENTS_SELECTION-reMapFormulaMom", (std::string) "REMAP[$score$,$scoreMomAve$,$scoreMomMax$]",
	"function used during roulette selection to remap scores genereated by optimizeValueMom."
	"\nUsed to control strength of selection and allows for negitive scores."
	"\nin addtion of MTree operators, the $score$, $scoreMomMin$, $scoreMomAve$, $scoreMomMax$, $scoreDadMin$, $scoreDadAve$ and, $scoreDaddMax$ can be used");
std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::reMapFormulaDadPL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES_PARENTS_SELECTION-reMapFormulaDad", (std::string) "REMAP[$score$,$scoreDadAve$,$scoreDadMax$]",
	"see optimizeFormulaMom");

std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::nextPopSizePL = Parameters::register_parameter("OPTIMIZER_TWO_SEXES-nextPopSize", (std::string)"-1", "size of population after optimization(MTree). -1 will maintian the current population size.");

std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::optimizeValueMateEffectDadPL = 
	Parameters::register_parameter("OPTIMIZER_TWO_SEXES_PARENTS_SELECTION-optimizeValueAfterDadMates",
	(std::string)"NONE", "when dad mates (i.e. is chosen ) even if no ofspring are produced, score will be temporarily updated using this MTree"
	" (NONE indicates no change) $score$ (current temporary score) can be used (e.g. $score$*.5 = 1/2 score after each mating)");
std::shared_ptr<ParameterLink<std::string>> TwoSexesOptimizer::optimizeValueMateEffectMomPL = 
	Parameters::register_parameter("OPTIMIZER_TWO_SEXES_PARENTS_SELECTION-optimizeValueAfterMomMates", 
	(std::string)"NONE", "when mom mates (even if no offspring are produced) or mom does not accept a mate (if system is Threshold), score will be temporarily updated using this MTree"
	" (NONE indicates no change) $score$ (current temporary score) can be used (e.g. $score$*.5 = 1/2 score after each mating)");


std::shared_ptr<ParameterLink<int>> TwoSexesOptimizer::litterSizePL =
	Parameters::register_parameter("OPTIMIZER_TWO_SEXES-litterSize", 
	1, "number of offspring produced from each mating event."
	"\nif condition effects are being used, then the % chance to succeed will be applied to each potential offspring.");


void stringReplace(std::string &s, const std::string &search, const std::string &replace) {
	for (size_t pos = 0; ; pos += replace.length()) {
		// Locate the substd::string to replace
		pos = s.find(search, pos);
		if (pos == std::string::npos) break;
		// Replace by erasing and inserting
		s.erase(pos, search.length());
		s.insert(pos, replace);
	}
}

std::string showUnsignedlonglongArray(unsigned long long arr[]) {
	std::string returnStr;
	for (size_t i = 0; i < 8; i++) {
		returnStr += std::to_string(arr[i]) + " ";
	}
	return returnStr;
}
std::string showIntDeque(std::deque<int> deq) {
	std::string returnStr;
	for (auto v : deq) {
		returnStr += std::to_string(v) + " ";
	}
	return returnStr;
}

void TwoSexesOptimizer::getGeneticMarkers(std::shared_ptr<Organism> org, std::shared_ptr<ParametersTable> PT_) {
	int alphabetSize = 256; // should this be parameter?
	int hashOffset = 8; // if so, this will need to be set to number of bits in alphabetSize

	// will contain counts for each additive maker
	std::vector<int> additiveCounts(additiveNames.size()); 

	// will contain lists of site values (in the range of the associated codon) for each assigned marker
	std::vector<std::vector<int>> assignedCollections(assignedNames.size()); 

	// will hold 9 values (8 to compare by hash and 1 extra for the assigned value with a 8 site assigned codon)
	std::deque<int> sitesBuffer;
	sitesBuffer.resize(9);
	// hash values for the 8 possible start codon sizes (size 1 to size 8)
	unsigned long long hashValues[8];

	// for each genome with atleast one associated genetic marker
	for (auto g : searchGenomes) {
		if (Global::update == 0) { // we may need to initalize genomes and/or seed codons
			if (std::find(seedRandomizeGenomes.begin(), seedRandomizeGenomes.end(), g) != seedRandomizeGenomes.end()) {
				// if this genome is in seedRandomizeGenomes then fillRandom
				//std::cout << "seeing genome " << g << std::endl;
				org->genomes[g]->fillRandom();
			}
			for (size_t i = 0; i < seedGenomes.size(); i++) {
				// for each genome in seedGenome...
				if (seedGenomes[i] == g) {
					// is this the genome that needs to be seeded?
					int copies = seedCodons[i][0]; // get the number of codons to seed
					int location = -1; // if -1, place codons at random locations - this is used if '@' and not 'x' is used in codon seeding
					if (copies <= 0) { // this is a location, not a count!
						location = -1*copies;
						copies = 1;
					}
					size_t codonSize = seedCodons[i].size() - 1; // get size (first value was number of copies)
					int genomeSize = org->genomes[g]->countSites(); // get size of this genome
					int offsetIntoGenome;
					for (int c = 0; c < copies; c++) {
						// go to random location within genome between 0 and genomeSize - (codonSize*2)
						auto handler = org->genomes[g]->newHandler(org->genomes[g]);
						if (location < 0) {
							offsetIntoGenome = Random::getIndex(genomeSize - (codonSize * 2)); // find a location between start and (end - (2 * size of this codon))
						}
						else {
							offsetIntoGenome = location;
						}
						handler->advanceIndex(offsetIntoGenome);
						// write values [1..end] from seedCodons[i] into genome
						for (size_t p = 0; p < codonSize; p++) {
							handler->writeInt(seedCodons[i][p + 1],0, alphabetSize-1);
						}
					}
				}
			}
		} // done with initialization / seeding

		auto genomeHandler = org->genomes[g]->newHandler(org->genomes[g]);
		auto genomeBufferHandler = genomeHandler->makeCopy(); // points to first site for first element in buffer.
		
		sitesBuffer.clear();
		

		// read 9 sites from genome, 8 for codons + 1 value incase of assigned codon with length 8
		// and load into sitesBuffer
		for (int i = 0; i < 8; i++) {
			sitesBuffer.push_back(genomeHandler->readInt(0, alphabetSize-1));
		}
		
		// we are now set up to be able to read the genome

		// fill the hash table
		hashValues[0] = sitesBuffer[0];
		for (int i = 1; i < 8; i++) {
			hashValues[i] = (hashValues[i - 1] << hashOffset) | sitesBuffer[i];
		}

		// check if the first site in genome defines an additive marker
		for (auto m : additiveMarkersSourceGenomes[g]) {
			if (additiveMarkersHashValues[m] == hashValues[additiveStartCodons[m].size() - 1]) {
				additiveCounts[m]++;
			}
		}

		// check if the first site in genome defines an assigned marker
		for (auto m : assignedMarkersSourceGenomes[g]) {
			if (assignedMarkersHashValues[m] == hashValues[assignedStartCodons[m].size() - 1]) {
				auto tempHandler = genomeBufferHandler->makeCopy();
				// from the first site in the buffer, read one value for each value in the start codon
				for (size_t s = 0; s < assignedStartCodons[m].size(); s++) {
					tempHandler->readInt(0, alphabetSize - 1);
				}				
				// now that the handler is positioned on the first site after the start codon
				// read a value from in the correct range.
				assignedCollections[m].push_back(tempHandler->readInt(0, assignedRanges[m]));
			}
		}

		// now we have addressed if there is a codon starting at the first site and set up the buffer.
		// read sites from genome one at a time. For each site recalculate the hash table and compair
		// markers (in this genome) start codons to the hash table entry with the same size.
		while (!genomeHandler->EOG) { // note: conversion will stop 9 sites before end of genome!
			// delete first value in buffer
			sitesBuffer.pop_front();
			// get next site from genome and add to the back of buffer
			sitesBuffer.push_back(genomeHandler->readInt(0, alphabetSize-1));
			genomeBufferHandler->readInt(0, alphabetSize - 1); // advance buffer handler

			// recalculate the hash table
			hashValues[0] = sitesBuffer[0];
			for (int i = 1; i < 8; i++) {
				hashValues[i] = (hashValues[i - 1] << hashOffset) | sitesBuffer[i];
			}
			
			//std::cout << "sites: " << showIntDeque(sitesBuffer) << " -> " << showUnsignedlonglongArray(hashValues) << std::endl;
			
			// check for additive marker
			for (auto m : additiveMarkersSourceGenomes[g]) {
				if (additiveMarkersHashValues[m] == hashValues[additiveStartCodons[m].size()-1]) {
					//std::cout << "found match with marker: " << additiveNames[m] << std::endl;
					additiveCounts[m]++;
				}
			}

			// check for assigned marker
			for (auto m : assignedMarkersSourceGenomes[g]) {
				if (assignedMarkersHashValues[m] == hashValues[assignedStartCodons[m].size() - 1]) {
					auto tempHandler = genomeBufferHandler->makeCopy();
					// from the first site in the buffer, read one value for each value in the start codon
					for (size_t s = 0; s < assignedStartCodons[m].size();s++) {
						tempHandler->readInt(0, alphabetSize - 1);
					}
					// now that the handler is positioned on the first site after the start codon
					// read a value from in the correct range.
					assignedCollections[m].push_back(tempHandler->readInt(0, assignedRanges[m]));
				}
			}
		}
	} // go onto the next genome

	bool debug = false;
	if (!debug) {
		// write markers information into orgs data map
		for (size_t i = 0; i < additiveNames.size(); i++) {
			org->dataMap.set(additiveNames[i], additiveCounts[i]);
		}
		for (size_t i = 0; i < assignedNames.size(); i++) {
			if (assignedCollections[i].size() == 0) {
				org->dataMap.set(assignedNames[i], -1);
			}
			else if (assignedDuplicatesRules[i] == "Unique") {
				if (assignedCollections[i].size() == 1) {
					org->dataMap.set(assignedNames[i], assignedCollections[i][0]);
					//std::cout << assignedNames[i] << " = " << org->dataMap.getAverage(assignedNames[i]) << std::endl;
				}
				else {
					org->dataMap.set(assignedNames[i], -1);
				}
			}
			else if (assignedDuplicatesRules[i] == "Most") {
				std::map<int, int> valueCounts; // value, number of times it appears
				int currentMost = assignedCollections[i][0]; // seed with first value
				for (size_t j = 0; j < assignedCollections[i].size(); j++) {
					valueCounts[assignedCollections[i][j]]++;
					if (valueCounts[assignedCollections[i][j]] >= valueCounts[currentMost]) {
						currentMost = assignedCollections[i][j];
					}
				}
				org->dataMap.set(assignedNames[i], currentMost);
			}
			else if (assignedDuplicatesRules[i] == "Average") {
				int total = 0;
				for (size_t j = 0; j < assignedCollections[i].size(); j++) {
					total += assignedCollections[i][j];
				}
				org->dataMap.set(assignedNames[i], (int)(total / (double)assignedCollections[i].size()));
			}
			else if (assignedDuplicatesRules[i] == "Greatest") {
				int Greatest = assignedCollections[i][0];
				for (size_t j = 0; j < assignedCollections[i].size(); j++) {
					Greatest = (assignedCollections[i][j] > Greatest) ? assignedCollections[i][j] : Greatest;
				}
				org->dataMap.set(assignedNames[i], Greatest);
			}
			else if (assignedDuplicatesRules[i] == "Least") {
				int Least = assignedCollections[i][0];
				for (size_t j = 0; j < assignedCollections[i].size(); j++) {
					Least = (assignedCollections[i][j] < Least) ? assignedCollections[i][j] : Least;
				}
				org->dataMap.set(assignedNames[i], Least);
			}
		}
	} // end write markers into DataMap
	else { // same thing, but print debug info
		std::cout << "these markers were found:" << std::endl;
		std::cout << "additive markers:" << std::endl;
		for (size_t i = 0; i < additiveCounts.size(); i++) {
			std::cout << "  " << additiveNames[i] << "  " << additiveCounts[i] << std::endl;
		}
		std::cout << "assigned markers:" << std::endl;
		for (size_t i = 0; i < assignedCollections.size(); i++) {
			std::cout << "  " << assignedNames[i] << "  values: ";
			for (size_t j = 0; j < assignedCollections[i].size(); j++) {
				std::cout << assignedCollections[i][j] << " ";
			}
			std::cout << std::endl;
		}

		// write markers information into orgs data map
		for (size_t i = 0; i < additiveNames.size(); i++) {
			std::cout << "ADD " << additiveNames[i] << " " << additiveCounts[i] << std::endl;
			org->dataMap.set(additiveNames[i], additiveCounts[i]);
		}
		for (size_t i = 0; i < assignedNames.size(); i++) {
			if (assignedCollections[i].size() == 0) {
				org->dataMap.set(assignedNames[i], -1);
				std::cout << "NONE " << assignedNames[i] << " " << assignedCollections[i].size() << std::endl;
			}
			else if (assignedDuplicatesRules[i] == "Unique") {
				if (assignedCollections[i].size() == 1) {
					org->dataMap.set(assignedNames[i], assignedCollections[i][0]);
					std::cout << "Unique " << assignedNames[i] << " found: " << assignedCollections[i].size() << " val: " << assignedCollections[i][0] << std::endl;
				}
				else {
					org->dataMap.set(assignedNames[i], -1);
					std::cout << "Unique " << assignedNames[i] << " found: " << assignedCollections[i].size() << " val: " << -1 << std::endl;
				}
			}
			else if (assignedDuplicatesRules[i] == "Most") {
				std::map<int, int> valueCounts; // value, number of times it appears
				int currentMost = assignedCollections[i][0]; // seed with first value
				for (size_t j = 0; j < assignedCollections[i].size(); j++) {
					valueCounts[assignedCollections[i][j]]++;
					if (valueCounts[assignedCollections[i][j]] >= valueCounts[currentMost]) {
						currentMost = assignedCollections[i][j];
					}
				}
				org->dataMap.set(assignedNames[i], currentMost);
				std::cout << "Most " << assignedNames[i] << " found: " << assignedCollections[i].size() << " val: " << currentMost << std::endl;
			}
			else if (assignedDuplicatesRules[i] == "Average") {
				int total = 0;
				for (size_t j = 0; j < assignedCollections[i].size(); j++) {
					total += assignedCollections[i][j];
				}
				org->dataMap.set(assignedNames[i], (int)(total / (double)assignedCollections[i].size()));
				std::cout << "Ave " << assignedNames[i] << " found: " << assignedCollections[i].size() << " val: " << (int)(total / (double)assignedCollections[i].size()) << std::endl;
			}
			else if (assignedDuplicatesRules[i] == "Greatest") {
				int Greatest = assignedCollections[i][0];
				for (size_t j = 0; j < assignedCollections[i].size(); j++) {
					Greatest = (assignedCollections[i][j] > Greatest) ? assignedCollections[i][j] : Greatest;
				}
				org->dataMap.set(assignedNames[i], Greatest);
				std::cout << "Greatest " << assignedNames[i] << " found: " << assignedCollections[i].size() << " val: " << Greatest << std::endl;
			}
			else if (assignedDuplicatesRules[i] == "Least") {
				int Least = assignedCollections[i][0];
				for (size_t j = 0; j < assignedCollections[i].size(); j++) {
					Least = (assignedCollections[i][j] < Least) ? assignedCollections[i][j] : Least;
				}
				org->dataMap.set(assignedNames[i], Least);
				std::cout << "Least " << assignedNames[i] << " found: " << assignedCollections[i].size() << " val: " << Least << std::endl;
			}
		}
	} // end write markers into DataMap - debug version
}

TwoSexesOptimizer::TwoSexesOptimizer(std::shared_ptr<ParametersTable> _PT) : AbstractOptimizer(_PT) {
	tournamentSizeMom = tournamentSizeMomPL->get(PT);
	tournamentSizeDad = tournamentSizeDadPL->get(PT);

	selectionMethod = selectionMethodPL->get(PT);

	sexFormula = stringToMTree(sexFormulaPL->get(PT));
	preferenceFormula = stringToMTree(preferenceFormulaPL->get(PT));

	surviveRateMT = stringToMTree(surviveRatePL->get(PT));

	optimizeValueMomMT = stringToMTree(optimizeValueMomPL->get(PT));
	optimizeValueDadMT = stringToMTree(optimizeValueDadPL->get(PT));

	nextPopSizeMT = stringToMTree(nextPopSizePL->get(PT));
	litterSize = litterSizePL->get(PT);

	// get sex trait formula
	convertCSVListToVector(sexTraitsPL->get(), traits, ',','"');
//	for (size_t i = 0; i < traits.size(); i++) {
//		std::cout << traits[i] << std::endl;
//	}

	// for each trait foruma string, determine if there are % values (error values) and create trait formula
	double tempDouble;
	for (size_t i = 0; i < traits.size(); i++) {
		//std::cout << "found trait: " << traits[i] << std::endl;
		if (traits[i] == "NONE") { // if one trait option is NONE, replace with 0 (i.e. random choice) formula
			traits[i] = "0";
		}
		// does this value contain %?
		if (traits[i].find("%") != std::string::npos) { // if entry contains % then a accuracy rate and (maybe) fail rate
														// accuracy determines if female experiances error when selecting from lek
														// fail rate is chance that female simply gets lek of size 1 (i.e. one male pick)
			traitFormulas.push_back(stringToMTree(traits[i].substr(0, traits[i].find("%")))); // add formula

			std::string afterPercent = traits[i].substr(traits[i].find("%") + 1);
			if (afterPercent.find("%") != std::string::npos) { // there are 2 "%" first is accuracy, second is fail rate
				convertString(afterPercent.substr(0, afterPercent.find("%")), tempDouble);
				traitsAccuracyRates.push_back(tempDouble);
				convertString(afterPercent.substr(afterPercent.find("%") + 1), tempDouble);
				traitsFailRates.push_back(tempDouble);
			}
			else { // there is only one %, accuracy; and fail rate = 0.0
				convertString(afterPercent, tempDouble);
				traitsAccuracyRates.push_back(tempDouble);
				traitsFailRates.push_back(0.0); // there is only one value, so there is no chance for a fail
			}
		}
		else { // no % provided, set accuracy and fail rate to 0
			traitFormulas.push_back(stringToMTree(traits[i])); // add formula
			traitsAccuracyRates.push_back(0.0); // no error
			traitsFailRates.push_back(0.0); // no chance for fail
		}

	}

	// get condition effects
	convertCSVListToVector(traitConditionEffectsPL->get(PT), traitConditionEffects, ',', '"');
	if (traitConditionEffects.size() == 1) {
		for (size_t i = 0; i < traitFormulas.size(); i++) { // if only one, then all will have that value
			traitConditionEffects.push_back(traitConditionEffects[0]);
		}
	}
	else if (traitConditionEffects.size() == 0) { // if none, then all will be 0, no effect
		traitConditionEffects.resize(traits.size());
	} 
	else if (traitConditionEffects.size() != traits.size()) {
		std::cout << "  in TwoSexesOptimizer construction:  traitConditionEffects must have 0 or 1 element or the same number as traits, but it has " << traitConditionEffects.size() << " elements.\n  exiting..." << std::endl;
		exit(1);
	} // else the list are the same length, all is good.

	// get selection system for sexual traits
	convertCSVListToVector(selectionSystemsPL->get(PT), selectionSystems, ',', '"');
	if (selectionSystems.size() == 1) {
		for (size_t i = 1; i < traitFormulas.size(); i++) {
			selectionSystems.push_back(selectionSystems[0]);
		}
	}
	else if (selectionSystems.size() != traitFormulas.size()) {
		std::cout << "  in TwoSexesOptimizer construction:  selectionSystems must have 1 element or the same number as traits, but it has " << selectionSystems.size() << " elements.\n  exiting..." << std::endl;
		exit(1);
	}

	int thresholdTraitCount = 0;
	for (size_t i = 0; i < selectionSystems.size(); i++) {
		if (selectionSystems[i] == "Threshold") {
			thresholdTraitCount++;
		}
	}


	// get MTrees for sexualTraitsThresholds
	std::vector<std::string> sexualTraitsThresholdsStrings;
	if (thresholdTraitCount == 0) { // there are no threshold traits
		sexualTraitsThresholdFormulas.resize(traitFormulas.size()); // placeholder formula
	}
	else { // there are threshold traits
		convertCSVListToVector(sexTraitsThresholdsPL->get(PT), sexualTraitsThresholdsStrings);
		for (auto sTTS : sexualTraitsThresholdsStrings) {
			sexualTraitsThresholdFormulas.push_back(stringToMTree(sTTS));
		}
		if (sexualTraitsThresholdFormulas.size() == 1) {
			for (size_t i = 1; i < traitFormulas.size(); i++) {
				sexualTraitsThresholdFormulas.push_back(sexualTraitsThresholdFormulas[0]);
				sexualTraitsThresholdsStrings.push_back(sexualTraitsThresholdsStrings[0]);
			}
		}
		else if (sexualTraitsThresholdFormulas.size() != traitFormulas.size()) {
			std::cout << "  in TwoSexesOptimizer construction:  sexualTraitsThresholds must have 1 element or the same number as traits, but it has " << sexualTraitsThresholdFormulas.size() << " elements.\n  exiting..." << std::endl;
			exit(1);
		}
	}

	std::cout << "\nWhile setting up TwoSexesOptimizer\n\n  found these sex traits:" << std::endl;
	for (size_t i = 0; i < traitFormulas.size(); i++) {
		std::string thresholdMessage = "";
		if (selectionSystems[i] == "Threshold") {
			thresholdMessage = " (" + sexualTraitsThresholdFormulas[i]->getFormula() + ")";
		}
		std::cout << "  * trait_" << i << "  " << traitFormulas[i]->getFormula() <<
			"\n      type: " << selectionSystems[i] << thresholdMessage <<
			"\n      traitError: " << traitsAccuracyRates[i] <<
			"    traitFailRate: " << traitsFailRates[i] << 
			"    traitConditionEffects: " << traitConditionEffects[i] << std::endl;
	}

	for (size_t i = 0; i < selectionSystems.size(); i++) {
		if (selectionSystems[i] == "Threshold") {
			thresholdIndices.push_back(i);
		}
	}

	selectionGroupSize = selectionGroupSizePL->get(PT);

	// optimizeFormulaDad and optimizeFormulaMom allow
	// for remaping of optimizeValueDad and optimizeValueDad to allow
	// for better control (dynamic to population ranges) on the
	// strength of selection. This usese the VECT type from MTree
	// When the formula is evaluated, the VECT object is loaded
	// with various values. Rather then force the user to know
	// that VECT[0,1] is the average mom score, we provide names
	// which the user can employ when createing a optimizeFormulaDad
	// or optimizeFormulaMom (i.e. a remapper). Here we substitute
	// the user readable names with the correct VECT entires.

	std::string temp;
	temp = reMapFormulaMomPL->get(PT);
	stringReplace(temp, "$score$", "VECT[0,6]");
	stringReplace(temp, "$scoreMomMin$", "VECT[0,0]");
	stringReplace(temp, "$scoreMomAve$", "VECT[0,1]");
	stringReplace(temp, "$scoreMomMax$", "VECT[0,2]");
	stringReplace(temp, "$scoreDadMin$", "VECT[0,3]");
	stringReplace(temp, "$scoreDadAve$", "VECT[0,4]");
	stringReplace(temp, "$scoreDadMax$", "VECT[0,5]");
	reMapFormulaMom = stringToMTree(temp);

	temp = reMapFormulaDadPL->get(PT);
	stringReplace(temp, "$score$", "VECT[0,6]");
	stringReplace(temp, "$scoreMomMin$", "VECT[0,0]");
	stringReplace(temp, "$scoreMomAve$", "VECT[0,1]");
	stringReplace(temp, "$scoreMomMax$", "VECT[0,2]");
	stringReplace(temp, "$scoreDadMin$", "VECT[0,3]");
	stringReplace(temp, "$scoreDadAve$", "VECT[0,4]");
	stringReplace(temp, "$scoreDadMax$", "VECT[0,5]");
	reMapFormulaDad = stringToMTree(temp);
	// optimizeFormulaDad and optimizeFormulaMom substitutions are done.



	// check for - in optimizeValueMateEffectDadPL and optimizeValueMateEffectMomPL
	updateDadScoresOnRepro = (optimizeValueMateEffectDadPL->get(PT) != "NONE");
	updateMomScoresOnRepro = (optimizeValueMateEffectMomPL->get(PT) != "NONE");
	// if optimizeValueMateEffectDadPL and/or optimizeValueMateEffectMomPL are being used...
	if (updateDadScoresOnRepro) {
		temp = optimizeValueMateEffectDadPL->get(PT);
		stringReplace(temp, "$score$", "VECT[0,0]");
		optimizeValueMateEffectDadMT = stringToMTree(temp);
	}
	if (updateMomScoresOnRepro) {
		temp = optimizeValueMateEffectMomPL->get(PT);
		stringReplace(temp, "$score$", "VECT[0,0]");
		optimizeValueMateEffectMomMT = stringToMTree(temp);
	}




	convertCSVListToVector(seedGenomesPL->get(PT), seedGenomes);
	convertCSVListToVector(seedRandomizeGenomesPL->get(PT), seedRandomizeGenomes);

	// seedCodons has format copiesxs1:s2:...,copies+s1:s2...
	std::vector<std::string> codonDefStrings;
	std::vector<std::string> oneCodonString;
	std::vector<int> defCodonValues;
	convertCSVListToVector(seedCodonsPL->get(PT), codonDefStrings);
	// temp data now has each codon, we must strip off the copies

	for (auto codonDef : codonDefStrings) {
		std::string newcondonString;

		if (codonDef.find('x') != std::string::npos){ // if codonDef has an x, this codon be seeded some number of times to random locaitons
			convertCSVListToVector(codonDef, oneCodonString,'x');
			newcondonString = oneCodonString[0] + ":" + oneCodonString[1];
		}
		else if (codonDef.find('<') != std::string::npos) { // if codonDef has an @, this codon be seeded once, this number of sites from start of genome
			convertCSVListToVector(codonDef, oneCodonString, '<');
			newcondonString = "-" + oneCodonString[0] + ":" + oneCodonString[1]; // negatiave count (first value) here indicates that this is a location in genome
		}
		else { // there is no x or <, not allowed
			std::cout << "found bad values in seedCodons parameter that had neither 'x' or '<'\n   paramater value: \"" << seedCodonsPL->get(PT) << "\".\n   exiting.";
			exit(1);
		}
		seedCodons.push_back({});
		convertCSVListToVector(newcondonString, seedCodons.back(),':');
	}

	for (auto sc : seedCodons) {
		for (auto elem : sc) {
			std::cout << elem << " ";
		}
		std::cout << " | ";
	}
	std::cout << std::endl;
	// localize genetic markers parameters
	convertCSVListToVector(TwoSexesOptimizer::additiveNamesPL->get(PT), additiveNames);

	convertCSVListToVector(TwoSexesOptimizer::additiveSourceGenomesPL->get(PT), additiveSourceGenomes);

	std::vector<std::string> additiveStartCodonsStrings;
	convertCSVListToVector(TwoSexesOptimizer::additiveStartCodonsPL->get(PT), additiveStartCodonsStrings);
	for (auto codon : additiveStartCodonsStrings) {
		additiveStartCodons.push_back({});
		convertCSVListToVector(codon, additiveStartCodons.back(), ':');
	}

	convertCSVListToVector(assignedNamesPL->get(PT), assignedNames);

	convertCSVListToVector(assignedSourceGenomesPL->get(PT), assignedSourceGenomes);

	std::cout << "assignedNamesPL->get(PT) = " << assignedNamesPL->get(PT) << std::endl;
	for (auto xx : assignedNames) {
		std::cout << xx << " : ";
	}
	std::cout << std::endl;
	std::cout << "assignedSourceGenomesPL->get(PT) = " << assignedSourceGenomesPL->get(PT) << std::endl;
	for (auto xx : assignedSourceGenomes) {
		std::cout << xx << " ";
	}
	std::cout << std::endl;

	std::vector<std::string> assignedStartCodonsStrings;
	convertCSVListToVector(TwoSexesOptimizer::assignedStartCodonsPL->get(PT), assignedStartCodonsStrings);
	for (auto codon : assignedStartCodonsStrings) {
		assignedStartCodons.push_back({});
		convertCSVListToVector(codon, assignedStartCodons.back(), ':');
	}

	convertCSVListToVector(TwoSexesOptimizer::assignedRangesPL->get(PT), assignedRanges);

	convertCSVListToVector(TwoSexesOptimizer::assignedDuplicatesRulesPL->get(PT), assignedDuplicatesRules);

	// make lists per type (additive/assigned) for each genome that is
	// associated with atleast one marker. generate a map with genome name
	// as keys and marker indices as values
	// while we are iterating over the markers, also use startcodons to generate hash values
	// durring genome reading, hash values will be generated and compaiered rather
	// then comparing all sites individually.
	for (size_t i = 0; i < additiveNames.size(); i++) {
		additiveMarkersSourceGenomes[additiveSourceGenomes[i]].push_back(i);
		additiveMarkersHashValues.push_back(0);
		for (auto v : additiveStartCodons[i]) {
			additiveMarkersHashValues.back() = (additiveMarkersHashValues.back() << 8) | v;
		}
	}

	for (size_t i = 0; i < assignedNames.size(); i++) {
		assignedMarkersSourceGenomes[assignedSourceGenomes[i]].push_back(i);
		assignedMarkersHashValues.push_back(0);
		for (auto v : assignedStartCodons[i]) {
			assignedMarkersHashValues.back() = (assignedMarkersHashValues.back() << 8) | v;
		}
	}

	// report what's been done
	std::cout << "\n  found the following additive markers" << std::endl;
	for (size_t m = 0; m < additiveNames.size(); m++) {
		std::cout << "  * " << additiveNames[m] << "  from genome: " << additiveSourceGenomes[m] << "  startCodon: ";
		for (auto v : additiveStartCodons[m]) {
			std::cout << v << " ";
		}
		std::cout << " hashValue(" << additiveMarkersHashValues[m] << ")" << std::endl;
	}
	std::cout << "\n  found the following assigned markers" << std::endl;
	for (size_t m = 0; m < assignedNames.size(); m++) {
		std::cout << "  * " << assignedNames[m] << "  from genome: " << assignedSourceGenomes[m] << "  startCodon: ";
		for (auto v : assignedStartCodons[m]) {
			std::cout << v << " ";
		}
		std::cout << "  hashValue(" << assignedMarkersHashValues[m] << ")" << std::endl;
	}

	// get genome names we need to search and put into searchGenomes vector
	for (auto g : additiveMarkersSourceGenomes) {
		if (std::find(searchGenomes.begin(), searchGenomes.end(), g.first) == searchGenomes.end()) {
			searchGenomes.push_back(g.first);
		}
	}
	for (auto g : assignedMarkersSourceGenomes) {
		if (std::find(searchGenomes.begin(), searchGenomes.end(), g.first) == searchGenomes.end()) {
			searchGenomes.push_back(g.first);
		}
	}



	popFileColumns.clear();
	popFileColumns.push_back("optimizeValueMom");
	popFileColumns.push_back("optimizeValueMom_VAR");
	popFileColumns.push_back("optimizeValueDad");
	popFileColumns.push_back("optimizeValueDad_VAR");
	popFileColumns.push_back("sex");
	popFileColumns.push_back("sex_VAR");
	popFileColumns.push_back("sexPrefrence");
	popFileColumns.push_back("sexPrefrence_VAR");
	for (int i = 0; i < (int)traits.size(); i++) {
		popFileColumns.push_back("sexPrefrence_" + std::to_string(i));
		popFileColumns.push_back("sexPrefrence_" + std::to_string(i) + "_VAR");
	}
	for (size_t j = 0; j < traitFormulas.size(); j++) { // get values for all traits for all orgs in current pop
		popFileColumns.push_back("trait_" + std::to_string(j));
		if (selectionSystems[j] == "Threshold") {
			popFileColumns.push_back("traitThreshold_" + std::to_string(j));
		}
	}

	
	for (auto marker : additiveNames) {
		popFileColumns.push_back(marker);
	}
	for (auto marker : assignedNames) {
		popFileColumns.push_back(marker);
	}
	

}

void TwoSexesOptimizer::optimize(std::vector<std::shared_ptr<Organism>> &population) {

	oldPopulationSize = (int)population.size();

	nextPopulationTargetSize = nextPopSizeMT->eval(PT)[0];// popSizeLPL->get(); // how big should the population be (after optimizer cleanup
	if (nextPopulationTargetSize == -1) {
		nextPopulationTargetSize = oldPopulationSize;
	}

	nextPopulationSize = 0;

	//vector<std::shared_ptr<SexyOrganism>> sexyPopulation;
	std::vector<std::shared_ptr<SexyOrganism>> sexyMalePopulation;
	std::vector<std::shared_ptr<SexyOrganism>> sexyFemalePopulation;

	// if this is the start of the run. prepare the genomes with display values, sex and preferances
	// NOTE: because this is not happening with an initializeGenome function (like brains) the 0 generation
	// organisms have been evaluated. The scores may not be correct, but this effect should not have any
	// siginifigent effects.

	if (Global::update == 0) {
	}
	for (auto org : population) {
		getGeneticMarkers(org, PT);
	}


	killList.clear(); // noone is dead (yet)
	int surviveCount = 0; // noone is surviving (yet)

	for (int i = 0; i < oldPopulationSize; i++) {

		std::shared_ptr<SexyOrganism> newSexyOrg = std::make_shared<SexyOrganism>();

		newSexyOrg->org = population[i]; // link to org
		newSexyOrg->sex = sexFormula->eval(newSexyOrg->org->dataMap, PT)[0];
		newSexyOrg->preference = preferenceFormula->eval(newSexyOrg->org->dataMap, PT)[0];


		//std::cout << " newSexyOrg " << i << " " << newSexyOrg->org->ID << " sex: " << newSexyOrg->sex << " pref: " << newSexyOrg->preference << std::endl;

		newSexyOrg->traitScores.clear();
		for (size_t traitIndex = 0; traitIndex < traitFormulas.size(); traitIndex++) { // get values for all traits for all orgs in current pop
			// get trait values
			newSexyOrg->traitScores.push_back(traitFormulas[traitIndex]->eval(population[i]->dataMap, PT)[0]);
			newSexyOrg->org->dataMap.set("trait_" + std::to_string(traitIndex), newSexyOrg->traitScores[traitIndex]);
			// get trait tresholds
			newSexyOrg->traitThresholds.push_back(0.0); // placeholder
			if (selectionSystems[traitIndex] == "Threshold") {
				newSexyOrg->traitThresholds[traitIndex] = sexualTraitsThresholdFormulas[traitIndex]->eval(newSexyOrg->org->dataMap, PT)[0];
				newSexyOrg->org->dataMap.append("traitThreshold_" + std::to_string(traitIndex), newSexyOrg->traitThresholds[traitIndex]);
			}
		}

		// org must have a sex and a pref and a threshold for each Threshold Trait

		bool hasThresholds = true;
		// if there are no threshold traits then this org has all needed thresholds

		for (auto thresholdIndex : thresholdIndices) { // if more then one thresholdIndices, then we must check to make sure this org has the needed thresholds
			if (newSexyOrg->traitThresholds[thresholdIndex] == -1) {
				hasThresholds = false;
			}
		}

		if (newSexyOrg->sex == 0 && newSexyOrg->preference != -1 && hasThresholds) {
			sexyMalePopulation.push_back(newSexyOrg);
			newSexyOrg->org->dataMap.set("sex", 1);
			//goodCount++;
			//(newSexyOrg->preference == 0) ? prefP++ : prefN++;
		}
		else if (newSexyOrg->sex == 1 && newSexyOrg->preference != -1 && hasThresholds) {
			sexyFemalePopulation.push_back(newSexyOrg);
			newSexyOrg->org->dataMap.set("sex", -1);
			//goodCount++;
			//(newSexyOrg->preference == 0) ? prefP++ : prefN++;
		}
		else {
			// not viable (missing sex and/or Preference)
			newSexyOrg->org->dataMap.set("sex", 0); // record no sex.
			//fatalCount++;
		}

		population[i]->dataMap.set("sexPrefrence", newSexyOrg->preference);
		for (int t = 0; t < (int)traits.size(); t++) {
			std::string name = "sexPrefrence_" + std::to_string(t);
			int thisPref = newSexyOrg->preference == t;
			newSexyOrg->org->dataMap.append(name, thisPref);
		}


		// get primary value
		newSexyOrg->scoreMom = optimizeValueMomMT->eval(population[i]->dataMap, PT)[0];
		newSexyOrg->scoreDad = optimizeValueDadMT->eval(population[i]->dataMap, PT)[0];

		newSexyOrg->org->dataMap.set("optimizeValueMom", newSexyOrg->scoreMom);
		newSexyOrg->org->dataMap.set("optimizeValueDad", newSexyOrg->scoreDad);


		if (Random::P(surviveRateMT->eval(newSexyOrg->org->dataMap, PT)[0])) {
			surviveCount++;
			nextPopulationSize++;
		}
		else {
			killList.insert(population[i]);// add all non-surviving orgs to kill list (to be cleaned up later)
		}


	}

	// parents have been sorted into sexes and all values have been generated
	// now get stats (min,max,ave,best...)

	int femalePopulationSize = sexyFemalePopulation.size();
	int malePopulationSize = sexyMalePopulation.size();


	/*
		if (sexyFemalePopulation.size() < 1 || sexyMalePopulation.size() < selectionGroupSize) {
			std::cout << "\n\n***********************\npopulation is not viable or less males then lek size" << std::endl;
			std::cout << "moms: " << sexyFemalePopulation.size() << "  dads: " << sexyMalePopulation.size() << std::endl;
			exit(0);
		}
	*/

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Possible issue when either no males or no males that can meet and females prefThreshold. 
	//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (sexyFemalePopulation.size() < 1) {
		std::cout << "\n\n***********************\npopulation is not viable there are no females." << std::endl;
		std::cout << "moms: " << sexyFemalePopulation.size() << "  dads: " << sexyMalePopulation.size() << std::endl;
		exit(0);
	}




	////////////
	// now cull
	////////////

	/*std::vector<double> allMomScores;
	std::vector<double> allDadScores;

	// get primary score ranges and trait ranges
	for (auto M : sexyFemalePopulation) {
		allMomScores.push_back(M->scoreMom);
	}
	for (auto D : sexyMalePopulation) {
		allDadScores.push_back(D->scoreDad);
	}

	double cull_index;


	std::cout << sexyFemalePopulation.size() << "  " << sexyMalePopulation.size() << std::endl;
	double momCullBelow = .05;
	cull_index = std::ceil(std::max(((momCullBelow * allMomScores.size()) - 1.0), 0.0));
	std::nth_element(std::begin(allMomScores),
		std::begin(allMomScores) + cull_index,
		std::end(allMomScores));
	double momCullBelowScore = allMomScores[cull_index];

	double dadCullBelow = .05;

	std::cout << dadCullBelow << "  " << allDadScores.size() << " -> " << std::ceil(std::max(((dadCullBelow * allDadScores.size()) - 1.0), 0.0)) << std::endl;
	cull_index = std::ceil(std::max(((dadCullBelow * allDadScores.size()) - 1.0), 0.0));
	std::nth_element(std::begin(allDadScores),
		std::begin(allDadScores) + cull_index,
		std::end(allDadScores));
	double dadCullBelowScore = allDadScores[cull_index];
	std::cout << dadCullBelowScore << " DCBS" << std::endl;

	int momScoresAllSame = -1;
	int dadScoresAllSame = -1;
	double aScore;
	std::vector<std::shared_ptr<SexyOrganism>> culled_sexyFemalePopulation;
	for (auto M : sexyFemalePopulation) {
		if (M->scoreMom >= momCullBelowScore) {
			// must check if all scores are the same
			if (momScoresAllSame == -1) {
				aScore = M->scoreMom;
				momScoresAllSame = 1;
			}
			else {
				if (aScore != M->scoreMom) {
					momScoresAllSame = 0;
				}
			}
			culled_sexyFemalePopulation.push_back(M);
		}
	}

	sexyFemalePopulation = culled_sexyFemalePopulation;

	std::vector<std::shared_ptr<SexyOrganism>> culled_sexyMalePopulation;
	for (auto D : sexyMalePopulation) {
		if (D->scoreDad >= dadCullBelowScore) {
			// must check if all scores are the same
			if (dadScoresAllSame == -1) {
				aScore = D->scoreDad;
				dadScoresAllSame = 1;
			}
			else {
				if (aScore != D->scoreDad) {
					dadScoresAllSame = 0;
				}
			}
			culled_sexyMalePopulation.push_back(D);
		}
	}

	sexyMalePopulation = culled_sexyMalePopulation;
	std::cout << sexyFemalePopulation.size() << "  " << sexyMalePopulation.size() << std::endl;
	std::cout << "sameDads = " << dadScoresAllSame << "sameMoms = " << momScoresAllSame << std::endl;

	femalePopulationSize = sexyFemalePopulation.size();
	malePopulationSize = sexyMalePopulation.size();

	if (sexyFemalePopulation.size() < 1) {
		std::cout << "\n\n*********************** AFTER CULL\npopulation is not viable there are no females." << std::endl;
		std::cout << "moms: " << sexyFemalePopulation.size() << "  dads: " << sexyMalePopulation.size() << std::endl;
		exit(0);
	}
	*/
	// done with cull
	// collect population data for orgs that did not get culled

	scoreMomAve = 0;
	scoreDadAve = 0;
	scoreMomMax = sexyFemalePopulation[0]->scoreMom;
	scoreMomMin = scoreMomMax;
	scoreDadMax = sexyMalePopulation[0]->scoreDad;
	scoreDadMin = scoreDadMax;

	// make container for traitRanges;
	std::vector<double> minMomTraitValues;
	std::vector<double> maxMomTraitValues;

	std::vector<double> minDadTraitValues;
	std::vector<double> maxDadTraitValues;

	for (auto t : sexyFemalePopulation[0]->traitScores) {
		minMomTraitValues.push_back(t);
		maxMomTraitValues.push_back(t);
	}
	for (auto t : sexyMalePopulation[0]->traitScores) {
		minDadTraitValues.push_back(t);
		maxDadTraitValues.push_back(t);
	}

	// setup vars for best mom and best dad
	auto bestMom = sexyFemalePopulation[0];
	auto bestDad = sexyMalePopulation[0];

	// get primary score ranges and trait ranges
	for (auto M : sexyFemalePopulation) {
		M->scoreMomCurrent = M->scoreMom;
		scoreMomMax = std::max(M->scoreMom, scoreMomMax);
		scoreMomMin = std::min(M->scoreMom, scoreMomMin);
		scoreMomAve += M->scoreMom;
		if (M->scoreMom == scoreMomMax) {
			bestMom = M;
		}
		for (size_t tIndex = 0; tIndex < M->traitScores.size(); tIndex++) {
			maxMomTraitValues[tIndex] = std::max(M->traitScores[tIndex], maxMomTraitValues[tIndex]);
			minMomTraitValues[tIndex] = std::min(M->traitScores[tIndex], minMomTraitValues[tIndex]);
		}
	}

	for (auto D : sexyMalePopulation) {
		D->scoreDadCurrent = D->scoreDad;
		scoreDadMax = std::max(D->scoreDad, scoreDadMax);
		scoreDadMin = std::min(D->scoreDad, scoreDadMin);
		scoreDadAve += D->scoreDad;
		if (D->scoreDad == scoreDadMax) {
			bestDad = D;
		}
		for (size_t tIndex = 0; tIndex < D->traitScores.size(); tIndex++) {
			maxDadTraitValues[tIndex] = std::max(D->traitScores[tIndex], maxDadTraitValues[tIndex]);
			minDadTraitValues[tIndex] = std::min(D->traitScores[tIndex], minDadTraitValues[tIndex]);
		}
	}

	scoreMomAve /= femalePopulationSize;
	scoreDadAve /= malePopulationSize;

	double scoreMomAve_SAVE = scoreMomAve;
	double scoreDadAve_SAVE = scoreDadAve;

	std::cout << std::endl;
	std::cout << "\tMoms(" << femalePopulationSize << "):: max = " << std::to_string(scoreMomMax) << "   ave = " << std::to_string(scoreMomAve) << std::endl;
	std::cout << "\tDads(" << malePopulationSize << "):: max = " << std::to_string(scoreDadMax) << "   ave = " << std::to_string(scoreDadAve) << std::endl;





	// if elitism, then make sure that the best mom and dad survive
	if (elitismPL->get()) {
		if (killList.find(bestMom->org) != killList.end()) { // if this org is not already surviving
			killList.erase(bestMom->org);
			surviveCount++;
			nextPopulationSize++;
		}
		if (killList.find(bestDad->org) != killList.end()) { // if this org is not already surviving
			killList.erase(bestDad->org);
			surviveCount++;
			nextPopulationSize++;
		}
	}

	// now select parents for remainder of population
	std::vector<int> offspringCountsMoms;
	std::vector<int> offspringCountsDads;
	offspringCountsMoms.resize(femalePopulationSize);
	offspringCountsDads.resize(malePopulationSize);

	int momIndex; // mom index in sexyFemalePopulation
	int dadIndex; // dad index in sexyMalePopulation
	int lekIndex; // used when picking from lek
	int challengerIndex; // used in tournaments

	int thisLekSize; // size of the current lek, will be determined by parameter, mom or set to 1 if random
	std::vector<int> lekIndexes; //container to hold indexes for dads in leks

	while (nextPopulationSize < nextPopulationTargetSize) {  // while we have not filled up the next generation
		// if there are two sexes then the male and female populations will have been set up correctly
		// if there is only one sex (i.e. no defined sexes) then sexyPopulation, sexyMalePopulation and
		// sexyFemalePopulation will all be the same (i.e. mom and dad will be being pulled from the same
		// populations)

		// select a mom (chooser)

		// Roulette selection uses score-min/max-min. This fixes strength of selection. In order
		// to adjust strenth of selection, Skew is used. When naturalSelectionStrength factor is
		// 1, there is no effect. As this value is incresed, the strength of selection is decresed
		// (by decreasing the effective min).

		//std::cout << "top of loop" << std::endl;
		std::vector<std::vector<double>> scoreRanges = { { scoreMomMin,scoreMomAve,scoreMomMax,scoreDadMin,scoreDadAve,scoreDadMax,0 } }; // used by MTrees values are {scoreMomMin,scoreMomMax,scoreDadMin,scoreDadMax,score}

		if (scoreMomMax == scoreMomMin) {
			momIndex = Random::getIndex(femalePopulationSize); // pick a random mom
		}
		else {
			if (selectionMethod == "Tournament") {
				momIndex = Random::getIndex(femalePopulationSize); // pick a random mom
				for (int c = 0; c < tournamentSizeMom - 1; c++) {
					challengerIndex = Random::getIndex(femalePopulationSize); // get a challanger tournamnet size times
					if (sexyFemalePopulation[challengerIndex]->scoreMomCurrent > sexyFemalePopulation[momIndex]->scoreMomCurrent) { // if challenger is better the mom...
						momIndex = challengerIndex; // replace mom of lek with challenger
					}
				}
			}
			else if (selectionMethod == "Roulette") {
				do {
					momIndex = Random::getIndex(femalePopulationSize);  //keep choosing a random genome from population until we get one that's good enough
					scoreRanges[0][6] = sexyFemalePopulation[momIndex]->scoreMomCurrent;
					//std::cout << momIndex << " s: " << sexyFemalePopulation[momIndex]->scoreMom << " / max: " << scoreRanges[0][1] << " - min: " << scoreRanges[0][0];
					//std::cout << "   " << (sexyFemalePopulation[momIndex]->scoreMom - scoreRanges[0][0]) / (scoreRanges[0][1] - scoreRanges[0][0]) << std::endl;
				} while (!Random::P(reMapFormulaMom->eval(sexyFemalePopulation[momIndex]->org->dataMap, PT, scoreRanges)[0]));
				//std::cout << "+" << flush;
			}
			else {
				std::cout << "in TwoSexesOptimizer:: = unknown selection method " << selectionMethod << ". exitting..." << std::endl;
				exit(1);
			}

		}
		//std::cout << "got mom" << std::endl;

		// based on moms prference set lek size (i.e. set to size 1 if preference fails)		
		if (Random::P(traitsFailRates[sexyFemalePopulation[momIndex]->preference])) {
			thisLekSize = 1; // mom failed to select on preference, set this lek to size 1 (random dad)
		}
		else {
			thisLekSize = selectionGroupSize; // did not fail so mom gets to select from a lek!
		}

		// now select a lek / selection group for Threshold

		lekIndexes.clear();
		for (int n = 0; n < thisLekSize; n++) {
			if (scoreDadMax == scoreDadMin) {
				lekIndex = Random::getIndex(malePopulationSize); // pick a random dad
			}
			else {
				if (selectionMethod == "Tournament") {
					lekIndex = Random::getIndex(malePopulationSize); // pick a random dad
					for (int c = 0; c < tournamentSizeDad - 1; c++) {
						challengerIndex = Random::getIndex(malePopulationSize); // get a challanger tournamnet size times
						if (sexyMalePopulation[challengerIndex]->scoreDad > sexyMalePopulation[lekIndex]->scoreDad) { // if challenger is better the mom...
							lekIndex = challengerIndex; // replace mom of lek with challenger
						}
					}
				}
				else if (selectionMethod == "Roulette") {
					do {
						lekIndex = Random::getIndex(malePopulationSize); // pick a random dad
						scoreRanges[0][6] = sexyMalePopulation[lekIndex]->scoreDadCurrent;
						//std::cout << scoreDadMin << "  "  << scoreDadMax << "      " << sexyMalePopulation[lekIndex]->scoreDadCurrent << " -> " << dadFitnessFormula->eval(sexyMalePopulation[dadIndex]->org->dataMap, PT, scoreRanges)[0] << std::endl;
					} while (!Random::P(reMapFormulaDad->eval(sexyFemalePopulation[momIndex]->org->dataMap, PT, scoreRanges)[0]));
					//std::cout << "-" << std::flush;
				}

			}
			lekIndexes.push_back(lekIndex);
		}

		//std::cout << "got dads" << std::endl;



		// get selection system for this trait
		if (selectionSystems[sexyFemalePopulation[momIndex]->preference] == "Lek") {
			// now select best from lek based on moms prefrence

			/*std::cout << "Lek  m: " << momIndex <<
				" Mpref: " << sexyFemalePopulation[momIndex]->prefrence <<
				" MThreshold: " << sexyFemalePopulation[momIndex]->traitThresholds[sexyFemalePopulation[momIndex]->prefrence] <<
				std::endl;
			*/
			double selectionError = traitsAccuracyRates[sexyFemalePopulation[momIndex]->preference];
			if (!Random::P(selectionError)) { // will female select perfectly?
				int lekPick = 0; // set current lek pick to first org in lek
				for (size_t lekChallenger = 1; lekChallenger < thisLekSize; lekChallenger++) { // for each other org
					//std::cout << sexyMalePopulation[lekIndexes[lekChallenger]]->traitScores[sexyFemalePopulation[momIndex]->prefrence] << "  ### " << sexyMalePopulation[lekIndexes[lekPick]]->traitScores[sexyFemalePopulation[momIndex]->prefrence] << std::endl;
					if (sexyMalePopulation[lekIndexes[lekChallenger]]->traitScores[sexyFemalePopulation[momIndex]->preference] >
						sexyMalePopulation[lekIndexes[lekPick]]->traitScores[sexyFemalePopulation[momIndex]->preference]) { // if they are better on moms pref
						lekPick = lekChallenger; // then they are now the lek pick
					}
				}
				dadIndex = lekIndexes[lekPick];
			}
			else { // female will select propotionaly
				double lowestTraitValue = sexyMalePopulation[lekIndexes[0]]->traitScores[sexyFemalePopulation[momIndex]->preference];
				for (int lekIndex : lekIndexes) {
					lowestTraitValue = std::min(lowestTraitValue, sexyMalePopulation[lekIndex]->traitScores[sexyFemalePopulation[momIndex]->preference]);
					//std::cout << "  lowestTraitValue = " << lowestTraitValue << "   this: " << sexyMalePopulation[lekIndex]->traitScores[sexyFemalePopulation[momIndex]->prefrence] << std::endl;
				}


				double totalLekTraitValue = 0;
				for (int lekIndex = 0; lekIndex < lekIndexes.size(); lekIndex++) {
					totalLekTraitValue += sexyMalePopulation[lekIndexes[lekIndex]]->traitScores[sexyFemalePopulation[momIndex]->preference] - lowestTraitValue;
					//std::cout << "  totalLekTraitValue = " << totalLekTraitValue << std::endl;
				}
				int lekPick = 0;
				double runningTotal = sexyMalePopulation[lekIndexes[0]]->traitScores[sexyFemalePopulation[momIndex]->preference] - lowestTraitValue;
				double randomPickMagic = Random::getDouble(0, totalLekTraitValue);
				//std::cout << lekIndexes[lekPick] << "  runningTotal = " << runningTotal << std::endl;
				while (runningTotal < randomPickMagic) {
					lekPick++;
					runningTotal += sexyMalePopulation[lekIndexes[lekPick]]->traitScores[sexyFemalePopulation[momIndex]->preference] - lowestTraitValue;
					//std::cout << lekIndexes[lekPick] << "  runningTotal = " << runningTotal << "   this: " << sexyMalePopulation[lekIndexes[lekPick]]->traitScores[sexyFemalePopulation[momIndex]->prefrence] << std::endl;
				}
				dadIndex = lekIndexes[lekPick];
			}
		}
		else if (selectionSystems[sexyFemalePopulation[momIndex]->preference] == "Threshold") {
			int lekPick = 0; // set current lek pick to first org in lek
			dadIndex = -1;

			/*std::cout << "Threshold  m: " << momIndex <<
				" Mpref: " << sexyFemalePopulation[momIndex]->preference <<
				" MThreshold: " << sexyFemalePopulation[momIndex]->traitThresholds[sexyFemalePopulation[momIndex]->preference] <<
				std::endl;
			*/

			while ((lekPick < thisLekSize) && (dadIndex == -1)) {

				/*std::cout << "   lekPick: " << lekPick << "  is dad: " <<
					lekIndexes[lekPick] << "  is org: " << sexyMalePopulation[lekIndexes[lekPick]]->org->ID <<
					"  with trait: " << sexyFemalePopulation[momIndex]->preference <<
					"  at value: " << sexyMalePopulation[lekIndexes[lekPick]]->traitScores[sexyFemalePopulation[momIndex]->preference] <<
					std::endl;
					*/

					//if (sexyMalePopulation[lekIndexes[lekPick]]->traitScores[sexyFemalePopulation[momIndex]->prefrence] > (sexyFemalePopulation[momIndex]->traitThresholds[sexyFemalePopulation[momIndex]->prefrence] * (Random::getDouble(10)))) {
				if (sexyMalePopulation[lekIndexes[lekPick]]->traitScores[sexyFemalePopulation[momIndex]->preference] > sexyFemalePopulation[momIndex]->traitThresholds[sexyFemalePopulation[momIndex]->preference]) {
					//auto delta = (sexyMalePopulation[lekIndexes[lekPick]]->traitScores[sexyFemalePopulation[momIndex]->prefrence] - sexyFemalePopulation[momIndex]->traitThresholds[sexyFemalePopulation[momIndex]->prefrence]);
					//if (delta < 10 && ((delta - 1) > Random::getDouble(10))) {
					dadIndex = lekIndexes[lekPick];
					std::cout << "+" << std::flush;
				}
				else {
					// this dad fails threshold
					//std::cout << "-";
				}
				lekPick++;
			}



			/*std::cout << "Threshold  m: " << momIndex << " d: " << dadIndex <<
				" Mpref: " << sexyFemalePopulation[momIndex]->prefrence <<
				" MThreshold: " << sexyFemalePopulation[momIndex]->traitThresholds[sexyFemalePopulation[momIndex]->prefrence];

			if (dadIndex != -1) {
				std::cout <<
					" DVal: " << sexyMalePopulation[dadIndex]->traitScores[sexyFemalePopulation[momIndex]->prefrence];
			}
			std::cout << std::endl;*/

			if (dadIndex == -1) {
				if (updateMomScoresOnRepro) {
					std::vector<std::vector<double>> scoreVec = { { sexyFemalePopulation[momIndex]->scoreMomCurrent } };
					sexyFemalePopulation[momIndex]->scoreMomCurrent = optimizeValueMateEffectMomMT->eval(sexyFemalePopulation[momIndex]->org->dataMap, PT, scoreVec)[0];

					scoreMomMax = sexyFemalePopulation[0]->scoreMomCurrent;
					scoreMomMin = scoreMomMax;
					scoreMomAve = 0;

					for (auto M : sexyFemalePopulation) {
						scoreMomMax = std::max(M->scoreMomCurrent, scoreMomMax);
						scoreMomMin = std::min(M->scoreMomCurrent, scoreMomMin);
						scoreMomAve += M->scoreMomCurrent;
					}
					scoreMomAve /= sexyFemalePopulation.size();
				}
				//	dadIndex = lekIndexes[Random::getIndex(lekIndexes.size())];
			}

		}
		else {
			std::cout << "  in TwoSexesOptimizer:: with mom: " << momIndex << "  (" << sexyFemalePopulation[momIndex]->org->ID << ") with preference " << sexyFemalePopulation[momIndex]->preference << " optimize found an unknown selectionSystems value: " << selectionSystems[sexyFemalePopulation[momIndex]->preference] << "\n  exitting." << std::endl;
			exit(1);
		}
		//std::cout << "got a dad (maybe): " << dadIndex << std::endl;




		// we have mom and dad - we need to see if dad is good enough to reproduce
		// now get dads reletive traitValues
		if (dadIndex != -1) {
			//std::cout << "XX" << std::endl;
			double maxTraitConditionEffect = 0;
			bool positiveTraits = false;
			for (auto tce : traitConditionEffects) {
				if (tce > 0) {
					maxTraitConditionEffect += tce;
					positiveTraits = true;
				}
			}
			//std::cout << "maxTraitConditionEffect = " << maxTraitConditionEffect << std::endl;
			double chanceToSucceed = 0.0; // the default is success! (or it will be after we add 1)
			for (size_t tIndex = 0; tIndex < sexyMalePopulation[dadIndex]->traitScores.size(); tIndex++) {
				double thisP = 0;
				if (traitConditionEffects[tIndex] != 0) {
					if (maxDadTraitValues[tIndex] == minDadTraitValues[tIndex]) {
						thisP = 0.01; // points for showing up
					}
					else {
						thisP = ((((sexyMalePopulation[dadIndex]->traitScores[tIndex] - minDadTraitValues[tIndex])) /
							((maxDadTraitValues[tIndex] - minDadTraitValues[tIndex]))) *
							traitConditionEffects[tIndex]); // relitive score * how important is this score
					}
				}
				chanceToSucceed += thisP + .01;
			}
			if (positiveTraits == false) { // if there are only negitive or 0 values then add 1.
										   // if all values were 0 then this will result in 1 (always succeed)
										   // if there are neg values, they will be down from 1
				chanceToSucceed++;
			}
			//std::cout << "chanceToSucceed = " << chanceToSucceed << std::endl;
			for (int i = 0; i < litterSize; i++) {
				if (Random::P(chanceToSucceed)) {
					//std::cout << nextPopulationSize << "  mom: " << sexyFemalePopulation[momIndex]->org->ID << " with pref: " << sexyFemalePopulation[momIndex]->prefrence << "  + dad: " << sexyMalePopulation[dadIndex]->org->ID << " with traitScore: " << sexyMalePopulation[dadIndex]->traitScores[sexyFemalePopulation[momIndex]->prefrence] << std::endl;
					population.push_back(sexyFemalePopulation[momIndex]->org->makeMutatedOffspringFromMany({ sexyFemalePopulation[momIndex]->org ,sexyMalePopulation[dadIndex]->org }));
					//if (population[population.size() - 1]->genomes["root::"]->countSites() > sexyFemalePopulation[momIndex]->org->genomes["root::"]->countSites() && population[population.size() - 1]->genomes["root::"]->countSites() > sexyMalePopulation[dadIndex]->org->genomes["root::"]->countSites()) {
					//	std::cout << sexyFemalePopulation[momIndex]->org->genomes["root::"]->countSites() << " + " << sexyMalePopulation[dadIndex]->org->genomes["root::"]->countSites() << " = " << population[population.size() - 1]->genomes["root::"]->countSites() << std::endl << std::endl;
					//}
					nextPopulationSize++;
					offspringCountsMoms[momIndex]++;
					offspringCountsDads[dadIndex]++;

					//std::cout << "org with " << chanceToSucceed << " succeded." << std::endl;
				} // end chanceToSucceed
				else {
					//std::cout << "org with " << chanceToSucceed << " failed." << std::endl;
				}
			} // litter size

			// if they had the opertunity to mate, even if they did not produce,
			// reduce scores if updateMomScoresOnRepro and/or updateDadScoresOnRepro
			if (updateMomScoresOnRepro) {
				std::vector<std::vector<double>> scoreVec = { { sexyFemalePopulation[momIndex]->scoreMomCurrent } };
				sexyFemalePopulation[momIndex]->scoreMomCurrent = optimizeValueMateEffectMomMT->eval(sexyFemalePopulation[momIndex]->org->dataMap, PT, scoreVec)[0];

				scoreMomMax = sexyFemalePopulation[0]->scoreMomCurrent;
				scoreMomMin = scoreMomMax;
				scoreMomAve = 0;

				for (auto M : sexyFemalePopulation) {
					scoreMomMax = std::max(M->scoreMomCurrent, scoreMomMax);
					scoreMomMin = std::min(M->scoreMomCurrent, scoreMomMin);
					scoreMomAve += M->scoreMomCurrent;
				}
				scoreMomAve /= sexyFemalePopulation.size();
			}

			if (updateDadScoresOnRepro) {
				std::vector<std::vector<double>> scoreVec = { { sexyMalePopulation[dadIndex]->scoreDadCurrent } };
				sexyMalePopulation[dadIndex]->scoreDadCurrent = optimizeValueMateEffectDadMT->eval(sexyMalePopulation[dadIndex]->org->dataMap, PT, scoreVec)[0];

				scoreDadMax = sexyMalePopulation[0]->scoreDadCurrent;
				scoreDadMin = scoreDadMax;
				scoreDadAve = 0;

				for (auto D : sexyMalePopulation) {
					scoreDadMax = std::max(D->scoreDadCurrent, scoreDadMax);
					scoreDadMin = std::min(D->scoreDadCurrent, scoreDadMin);
					scoreDadAve += D->scoreDadCurrent;
				}
				scoreDadAve /= sexyMalePopulation.size();
			}
		}// dadIndex was -1, no dad was good enough for mom.
		//std::cout << "done repro step" << std::endl;
	}

	std::vector<int> maleTotals;
	std::vector<int> femaleTotals;
	maleTotals.resize(nextPopulationTargetSize + 1);
	femaleTotals.resize(nextPopulationTargetSize + 1);


	std::vector<double> averageDadTraits(traitFormulas.size(), 0);
	std::vector<double> averageMomTraits(traitFormulas.size(), 0);
	std::vector<double> averageDadSelectionValues(traitFormulas.size(), 0);
	std::vector<double> averageMomSelectionValues(traitFormulas.size(), 0);

	for (size_t ii = 0; ii < malePopulationSize; ii++) {
		//std::cout << ii << "  " << offspringCountsDads[ii]  << "  " << maleTotals[offspringCountsDads[ii]] << std::endl;
		maleTotals[offspringCountsDads[ii]]++;
		for (size_t t_index = 0; t_index < traitFormulas.size(); t_index++) {
			averageDadTraits[t_index] = sexyMalePopulation[ii]->traitScores[t_index];
			if (sexyMalePopulation[ii]->preference == t_index) {
				averageDadSelectionValues[t_index]++;
			}
		}
		//sexyMalePopulation[0]->preference;
	}
	for (size_t ii = 0; ii < femalePopulationSize; ii++) {
		femaleTotals[offspringCountsMoms[ii]]++;
		for (size_t t_index = 0; t_index < traitFormulas.size(); t_index++) {
			averageMomTraits[t_index] = sexyFemalePopulation[ii]->traitScores[t_index];
			if (sexyFemalePopulation[ii]->preference == t_index) {
				averageMomSelectionValues[t_index]++;
			}
		}
	}
	std::string averageDadTraitsList = "";
	std::string averageMomTraitsList = "";
	std::string averageDadSelectionValuesList = "";
	std::string averageMomSelectionValuesList = "";

	for (size_t t_index = 0; t_index < traitFormulas.size(); t_index++) {
		//averageDadTraits[t_index] /= malePopulationSize;
		averageDadTraitsList += std::to_string(averageDadTraits[t_index] / malePopulationSize) + "|"; 
		//averageMomTraits[t_index] /= femalePopulationSize;
		averageMomTraitsList += std::to_string(averageMomTraits[t_index] / femalePopulationSize) + "|";
		//averageDadSelectionValues[t_index] /= malePopulationSize;
		averageDadSelectionValuesList += std::to_string(averageDadSelectionValues[t_index] / malePopulationSize) + "|";
		//averageMomSelectionValues[t_index] /= femalePopulationSize;
		averageMomSelectionValuesList += std::to_string(averageMomSelectionValues[t_index] / malePopulationSize) + "|";
	}

	//averageDadTraitsList += "\"";
    //averageMomTraitsList += "\"";
	//averageDadSelectionValuesList += "\"";
	//averageMomSelectionValuesList += "\"";


	std::string dadOffspringTotalsList = "";
	std::string momOffspringTotalsList = "";

	std::cout << "Dad offspring Counts:" << std::endl;
	for (size_t ii = 0; ii <= nextPopulationTargetSize; ii++) {
		dadOffspringTotalsList += std::to_string(maleTotals[ii]) + "|";
		if (maleTotals[ii] > 0) {
			std::cout << "  " << maleTotals[ii] << " had " << ii << " offspring" << std::endl;
		}
	}
	std::cout << "Mom offspring Counts:" << std::endl;
	for (size_t ii = 0; ii <= nextPopulationTargetSize; ii++) {
		momOffspringTotalsList += std::to_string(femaleTotals[ii]) + "|";
		if (femaleTotals[ii] > 0) {
			std::cout << "  " << femaleTotals[ii] << " had " << ii << " offspring" << std::endl;
		}
	}
	dadOffspringTotalsList += "";
	momOffspringTotalsList += "";




	std::cout << "survivors: " << surviveCount << std::endl;

	// save some data about optimization
	//  dadOffspringTotalsList; //string with list of were index is number of offspring and value is number of orgs with this number of offspring
	//  momOffspringTotalsList;
	//  scoreDadAve_SAVE; // ave dad score
	//  scoreMomAve_SAVE; // ave mom score
	//averageDadTraits[t_index]    average of traits on dads
	//averageMomTraits[t_index]    average of traits on moms
	//averageDadSelectionValues
	//averageMomSelectionValues

	if (0){	
	std::string headerLine = "dadOffCounts,momOffCounts,dadScoreAve,momScoreAve,dadTraitAves,momTraitAves,dadSelectionAves,momSelectionAves";
	std::string dataLine = dadOffspringTotalsList + "," + momOffspringTotalsList +
		"," + std::to_string(scoreDadAve_SAVE) + "," + std::to_string(scoreMomAve_SAVE) + "," +
		averageDadTraitsList + "," +
		averageMomTraitsList + "," +
		averageDadSelectionValuesList + "," +
		averageMomSelectionValuesList;

	FileManager::writeToFile("TwoSexesOptimizerData.csv", dataLine, headerLine);
	}




	//for (auto org : population) {
		//org->dataMap.Set("Simple_numOffspring", org->offspringCount);
	//}
	//std::cout << "BBB" << std::endl;

}

