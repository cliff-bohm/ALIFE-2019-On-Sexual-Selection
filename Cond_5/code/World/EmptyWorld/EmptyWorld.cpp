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


#include "EmptyWorld.h"

std::shared_ptr<ParameterLink<std::string>> EmptyWorld::groupNamePL =
	Parameters::register_parameter("WORLD_EMPTY-groupNameSpace",
	(std::string) "root::", "namespace a group you want created. Orgs in this group will be empty if no other module adds brains or genomes");

