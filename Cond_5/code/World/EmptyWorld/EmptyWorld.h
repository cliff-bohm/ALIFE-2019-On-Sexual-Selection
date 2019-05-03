//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

// Place holder world for when you don't want a world

#pragma once

#include "../AbstractWorld.h"

#include <cstdlib>
#include <thread>
#include <vector>

class EmptyWorld : public AbstractWorld {

public:
	static std::shared_ptr<ParameterLink<std::string>> groupNamePL;


	EmptyWorld(std::shared_ptr<ParametersTable> PT_ = nullptr) : AbstractWorld(PT_) {};

	virtual ~EmptyWorld() = default;

	void evaluateSolo(std::shared_ptr<Organism> org, int analyze,
	  int visualize, int debug) {
	  //do nothing
  }
  virtual void evaluate(std::map<std::string, std::shared_ptr<Group>> &groups,
                        int analyze, int visualize, int debug) {
	//do nothing
  }

  virtual std::unordered_map<std::string, std::unordered_set<std::string>>
  requiredGroups() override {
    return {{groupNamePL->get(PT),
             {}}};
    // requires nothing
  }
};

