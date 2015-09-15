/*
 * FeatureManager.cpp
 *
 *  Created on: Aug 5, 2015
 *      Author: joppich
 */

#include "FeatureManager.h"

FeatureManager::FeatureManager()
	: m_pFeatures(new std::map<std::string, std::vector<SplitEvent*>* >())
{

}

FeatureManager::~FeatureManager() {

	delete m_pFeatures;

}

