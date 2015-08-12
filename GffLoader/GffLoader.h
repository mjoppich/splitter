/*
 * GffLoader.h
 *
 *  Created on: Jul 7, 2015
 *      Author: joppich
 */

#ifndef GFFLOADER_H_
#define GFFLOADER_H_

#include "GffEntry.h"

#include <map>
#include <vector>

class GffLoader {
public:
	GffLoader(
			std::string& sFileName,
			std::vector<std::string>* pIgnoreFeatures,
			std::string* pPrefix = new std::string("")
	);
	virtual ~GffLoader();

	std::vector<GffEntry*>* getEntriesForSeqName(std::string* pSeqName);

        GffEntry* getChromosome(std::string* pSeqName);
	std::vector<std::string>* getSeqNames();


protected:

	std::map<std::string, std::vector<GffEntry*>* >* pSortedGffEntries;
        std::vector<std::string>* m_pChromosomeNames;
        std::vector<GffEntry*>* m_pChromosomes;

private:

	std::vector<GffEntry*>* createEntriesForSeqName(std::string* pSeqName);
	std::vector<GffEntry*>* createIntrons( std::vector<GffEntry*>* pTranscriptElements );



	struct {
		bool operator()(GffEntry* pElem1, GffEntry* pElem2)
		{
			if (pElem1->getStart() < pElem2->getStart())
				return true;

			if (pElem1->getStart() > pElem2->getStart())
				return false;

			return pElem1->getLength() > pElem2->getLength();
		}
	} sSortEntriesAsc;

};

#endif /* GFFLOADER_H_ */
