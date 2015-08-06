/*
 * FeatureManager.h
 *
 *  Created on: Aug 5, 2015
 *      Author: joppich
 */

#ifndef FEATUREMANAGER_H_
#define FEATUREMANAGER_H_

#include <map>
#include <vector>
#include <string.h>
#include <ostream>

#include "SplitEvent.h"

class FeatureManager {
public:
	FeatureManager();
	virtual ~FeatureManager();

	std::vector<SplitEvent*>* getFeatureVector(std::string sFeature)
	{
		if (m_pFeatures->find(sFeature) == m_pFeatures->end())
		{
			m_pFeatures->insert( m_pFeatures->end(), std::pair<std::string, std::vector<SplitEvent*>*>(sFeature, new std::vector<SplitEvent*>()));
		}

		return m_pFeatures->find(sFeature)->second;
	}

	SplitEvent* addFeatureEntry(std::string sFeature, bam1_t* pRead, uint32_t iStart, uint32_t iEnd, uint32_t iLeftContained, uint32_t iRightContained)
	{
		std::vector<SplitEvent*>* pFeatureVec = this->getFeatureVector(sFeature);

		SplitEvent* pSplitEvent = new SplitEvent( iStart, iEnd, 0 );

		pSplitEvent = this->getFeatureFromVec(pFeatureVec, pSplitEvent);

		pSplitEvent->addEntry(pRead, iLeftContained, iRightContained);

		return pSplitEvent;
	}

	SplitEvent* getFeatureFromVec(std::vector<SplitEvent*>* pFeatureVec, SplitEvent* pFeature)
	{

		// searches backwards
		for (uint32_t i = pFeatureVec->size(); i > 0; --i)
		{

			SplitEvent* pCurrentFeature = pFeatureVec->at(i-1);

			if ((pCurrentFeature->getStart() == pFeature->getStart()) && (pCurrentFeature->getEnd() == pFeature->getEnd()))
			{
				delete pFeature;
				return pCurrentFeature;
			}

		}

		pFeatureVec->push_back(pFeature);

		return pFeature;
	}
/*
	void print(std::ostream* pStream)
	{
		std::ostream oStream = *pStream;

		oStream << "Bla" << std::endl;

	}
*/
private:

	std::map<std::string, std::vector<SplitEvent*>* >* m_pFeatures;

};

#endif /* FEATUREMANAGER_H_ */
