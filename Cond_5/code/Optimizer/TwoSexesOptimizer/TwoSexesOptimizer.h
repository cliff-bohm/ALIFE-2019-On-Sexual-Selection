//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License


#include "../AbstractOptimizer.h"
#include "../../Utilities/MTree.h"
#include "../../Utilities/CSV.h"

#include <iostream>
#include <sstream>
#include <deque>

class TwoSexesOptimizer : public AbstractOptimizer {
 public:

	 static std::shared_ptr<ParameterLink<std::string>> additiveNamesPL;
	 static std::shared_ptr<ParameterLink<std::string>> additiveSourceGenomesPL;
	 static std::shared_ptr<ParameterLink<std::string>> additiveStartCodonsPL;
	 static std::shared_ptr<ParameterLink<std::string>> assignedNamesPL;
	 static std::shared_ptr<ParameterLink<std::string>> assignedSourceGenomesPL;
	 static std::shared_ptr<ParameterLink<std::string>> assignedStartCodonsPL;
	 static std::shared_ptr<ParameterLink<std::string>> assignedRangesPL;
	 static std::shared_ptr<ParameterLink<std::string>> assignedDuplicatesRulesPL;

	 static std::shared_ptr<ParameterLink<std::string>> seedGenomesPL;
	 static std::shared_ptr<ParameterLink<std::string>> seedCodonsPL;
	 static std::shared_ptr<ParameterLink<std::string>> seedRandomizeGenomesPL;

	 std::vector<std::string> seedGenomes;
	 std::vector<std::vector<int>> seedCodons;
	 std::vector<std::string> seedRandomizeGenomes;


	//static shared_ptr<ParameterLink<string>> selectionMethodPL; // Roullette([Exp=1.05 or Lin=1]),Tounament(size=5)
	static std::shared_ptr<ParameterLink<int>> tournamentSizeMomPL;
	static std::shared_ptr<ParameterLink<int>> tournamentSizeDadPL;
	//static std::shared_ptr<ParameterLink<int>> sexSitePL; // -1 = no (anyone can mate with anyone), 0+ = read site in genome to get sex
	//static std::shared_ptr<ParameterLink<int>> choiceSitePL; // 0+ = read site in genome to get choosy
	static std::shared_ptr<ParameterLink<int>> selectionGroupSizePL;
	static std::shared_ptr<ParameterLink<std::string>> selectionSystemsPL;
	//static std::shared_ptr<ParameterLink<bool>> choosersPL; // 0 no, anyone can be a chooser, 1 yes, "females" picked as first
	static std::shared_ptr<ParameterLink<std::string>> sexTraitsPL; // list of traits from DM - may include "NONE" (or const) for random mate
	static std::shared_ptr<ParameterLink<std::string>> sexTraitsThresholdsPL;
	static std::shared_ptr<ParameterLink<std::string>> sexFormulaPL; // what value is used to pick moms (in tournament)
	static std::shared_ptr<ParameterLink<std::string>> preferenceFormulaPL; // what value is used to pick moms (in tournament)
	std::shared_ptr<Abstract_MTree> sexFormula;
	std::shared_ptr<Abstract_MTree> preferenceFormula;


	static std::shared_ptr<ParameterLink<std::string>> optimizeValueMomPL; // what value is used to pick moms (in tournament)
	static std::shared_ptr<ParameterLink<std::string>> optimizeValueDadPL; // what value is used to pick dads (in tournament)

	static std::shared_ptr<ParameterLink<std::string>> traitConditionEffectsPL; // vector with condition effects for each trait

	//static std::shared_ptr<ParameterLink<int>> numGenDispsPL;
	//static std::shared_ptr<ParameterLink<std::string>> genDispsBoundsPL;
	//std::vector<double> genDispsBounds;

	static std::shared_ptr<ParameterLink<std::string>> selectionMethodPL;
	std::string selectionMethod;

	//static std::shared_ptr<ParameterLink<std::string>> sexGenomePL;
	//static std::shared_ptr<ParameterLink<std::string>> choiceGenomePL;
	//static std::shared_ptr<ParameterLink<std::string>> displayGenomePL;

	static std::shared_ptr<ParameterLink<std::string>> surviveRatePL; // value between 0 and 1 chance to self (if more then one parent)
	std::shared_ptr<Abstract_MTree> surviveRateMT;

	//static std::shared_ptr<ParameterLink<int>> displayOffsetPL;
	//static std::shared_ptr<ParameterLink<int>> displayInitialCountPL;
	//static std::shared_ptr<ParameterLink<bool>> displayLocatedRandomPL;

	static std::shared_ptr<ParameterLink<bool>> elitismPL;


	//static std::shared_ptr<ParameterLink<int>> sexCodon1PL;// = 33; // sexCodon is 33,34
	//static std::shared_ptr<ParameterLink<int>> sexCodon2PL;// = 34;

	//static std::shared_ptr<ParameterLink<int>> choiceCodon1PL;// = 44; // choiceCodon is 44,45
	//static std::shared_ptr<ParameterLink<int>> choiceCodon2PL;// = 45;

	//static std::shared_ptr<ParameterLink<int>> displayCodonLowPL;// = 0; // display0 codon is displayCodonLow,displayCodonHigh
	//static std::shared_ptr<ParameterLink<int>> displayCodonHighPL;// = 11; // display1 codon is displayCodonLow+1,displayCodonHigh-1
	//static std::shared_ptr<ParameterLink<bool>> displayOffsetOnePL;// = 1; // if this is 0 then the first value of all display codons will be the same

	static std::shared_ptr<ParameterLink<std::string>> reMapFormulaDadPL;
	static std::shared_ptr<ParameterLink<std::string>> reMapFormulaMomPL;

	static std::shared_ptr<ParameterLink<std::string>> nextPopSizePL;

	static std::shared_ptr<ParameterLink<std::string>> optimizeValueMateEffectDadPL;
	static std::shared_ptr<ParameterLink<std::string>> optimizeValueMateEffectMomPL;
	std::shared_ptr<Abstract_MTree> optimizeValueMateEffectDadMT;
	std::shared_ptr<Abstract_MTree> optimizeValueMateEffectMomMT;
	bool updateDadScoresOnRepro;
	bool updateMomScoresOnRepro;

	int tournamentSizeMom;
	int tournamentSizeDad;
	int oldPopulationSize;
	int nextPopulationTargetSize;
	int nextPopulationSize;

	double scoreMomAve;
	double scoreDadAve;
	double scoreMomMax;
	double scoreDadMax;
	double scoreMomMin;
	double scoreDadMin;

	static std::shared_ptr<ParameterLink<int>> litterSizePL;
	int litterSize;

	std::vector<std::vector<double>> scoreRanges; // used by MTrees, contains up to date min and max natural selection values for moms and dads

	std::shared_ptr<Abstract_MTree> reMapFormulaMom;
	std::shared_ptr<Abstract_MTree> reMapFormulaDad;




	std::vector<std::string> traits;
	std::vector<std::shared_ptr<Abstract_MTree>> sexualTraitsThresholdFormulas;
	std::vector<std::string> selectionSystems;
	int selectionGroupSize;

	std::vector<std::shared_ptr<Abstract_MTree>> traitFormulas;
	std::vector<double> traitsFailRates;
	std::vector<double> traitsAccuracyRates;

	std::vector<double> traitConditionEffects;

	std::vector<std::shared_ptr<Abstract_MTree>> traitSearchCostFormulas;

	std::vector<int> thresholdIndices; // will have list of traits that are theshold traits

	class SexyOrganism {
	public:
		std::shared_ptr<Organism> org;
		double scoreMom,scoreMomCurrent;
		double scoreDad,scoreDadCurrent;
		int sex;
		int preference;
		std::vector<double> traitScores;
		std::vector<double> traitThresholds; // used only with Thresholds
	};

	std::shared_ptr<Abstract_MTree> optimizeValueMomMT;
	std::shared_ptr<Abstract_MTree> optimizeValueDadMT;
	std::shared_ptr<Abstract_MTree> nextPopSizeMT;

	std::shared_ptr<ParameterLink<int>> popSizeLPL;


	std::vector<std::string> additiveNames;
	std::vector<std::string> additiveSourceGenomes;
	std::vector<std::vector<int>> additiveStartCodons;
	std::vector<std::string> assignedNames;
	std::vector<std::string> assignedSourceGenomes;
	std::vector<std::vector<int>> assignedStartCodons;
	std::vector<int> assignedRanges;
	std::vector<std::string> assignedDuplicatesRules;

	// additiveMarkersSourceGenomes and assignedMarkersSourceGenomes map genomes to indexes for markers found in those genomes
	std::map<std::string, std::vector<int>> additiveMarkersSourceGenomes;
	std::map<std::string, std::vector<int>> assignedMarkersSourceGenomes;
	std::vector<unsigned long long> additiveMarkersHashValues;
	std::vector<unsigned long long> assignedMarkersHashValues;
	std::vector<std::string> searchGenomes;

	TwoSexesOptimizer(std::shared_ptr<ParametersTable> _PT = nullptr);
	
	virtual void optimize(std::vector<std::shared_ptr<Organism>> &population) override;

	virtual std::unordered_set<std::string> requiredGenomes() override {
		//return { sexGenomePL->get(),choiceGenomePL->get(),displayGenomePL->get() };
		std::unordered_set<std::string> returnValue;

		for (auto g : additiveMarkersSourceGenomes) {
			returnValue.insert(g.first);
		}
		for (auto g : assignedMarkersSourceGenomes) {
			returnValue.insert(g.first);
		}

		return returnValue;
	}

	void getGeneticMarkers(std::shared_ptr<Organism> org, std::shared_ptr<ParametersTable> PT_);

	//virtual string maxValueName() override {
	//	return (PT == nullptr) ? optimizeValuePL->lookup() : PT->lookupString("OPTIMIZER_Simple-optimizeValue");
	//}
};

