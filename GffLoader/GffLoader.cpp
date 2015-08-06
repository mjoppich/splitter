/*
 * GffLoader.cpp
 *
 *  Created on: Jul 7, 2015
 *      Author: joppich
 */

#include "GffLoader.h"

#include <iostream>
#include <algorithm>
#include <fstream>

GffLoader::GffLoader(	std::string sFileName, std::vector<std::string>* pIgnoreFeatures, std::string sPrefix) {

	std::ifstream oInputStream;
	std::string sLine;

	pSortedGffEntries = new std::map<std::string, std::vector<GffEntry*>* >();
	std::map<std::string, std::vector<GffEntry*>* >::iterator oIt;

	std::vector<GffEntry*>* pCurrentTranscript = NULL;

	oInputStream.open(sFileName.c_str());

	GffEntry* pCurrentGene = NULL;


	while (! oInputStream.eof() )
	{

		std::getline(oInputStream, sLine);

		if ((sLine[0] == '#') || (sLine[0] == '\n') || (sLine[0] == '\0'))
			continue;

		GffEntry* pEntry = new GffEntry( &sLine );
		std::string* pSeqName = pEntry->getSeqName();

		if (pSeqName == NULL)
			continue;

		std::string* pPrefixedSeqName = new std::string(sPrefix.c_str());
		pPrefixedSeqName->append(pSeqName->c_str());
		pEntry->setSeqName(pPrefixedSeqName);
		pSeqName = pEntry->getSeqName();

		std::string* pFeature = pEntry->getFeature();

		if ( pIgnoreFeatures != NULL)
		{
			std::vector<std::string>::iterator oSIt = std::find(pIgnoreFeatures->begin(), pIgnoreFeatures->end(), *pFeature);

			// gff entry feature is in ignore features list
			if (oSIt != pIgnoreFeatures->end())
			{
				delete pEntry;
				continue;
			}

		}


		if (pFeature->compare("gene") == 0)
		{
			pCurrentGene = pEntry;

			std::vector<GffEntry*>* pGffEntries = createEntriesForSeqName(pSeqName);
			pGffEntries->push_back( pEntry );

		} else {

			pCurrentGene->addChild( pEntry );
			continue;

		}




	}

	// now sort the single vectors
	// TODO this is nicely parallelisable
	for (oIt = pSortedGffEntries->begin(); oIt != pSortedGffEntries->end(); ++oIt)
	{

		std::vector<GffEntry*>* pGffEntries = oIt->second; // contains genes

		for (size_t i = 0; i < pGffEntries->size(); ++i)
		{
			GffEntry* pGeneEntry = pGffEntries->at(i);
		}

		std::sort( pGffEntries->begin(), pGffEntries->end(), GffLoader::sSortEntriesAsc );

	}

}

std::vector<GffEntry*>* GffLoader::createEntriesForSeqName(std::string* pSeqName)
{
	std::map<std::string, std::vector<GffEntry*>* >::iterator oIt;
	oIt = pSortedGffEntries->find( pSeqName->c_str() );

	if (oIt == pSortedGffEntries->end() )
	{
		// needs to be inserted
		std::vector<GffEntry*>* pGffEntries = new std::vector<GffEntry*>();

		pSortedGffEntries->insert( std::pair<std::string, std::vector<GffEntry*>* >( *pSeqName, pGffEntries ) );

		return pGffEntries;
	}

	return oIt->second;

}

std::vector<GffEntry*>* GffLoader::getEntriesForSeqName(std::string* pSeqName)
{
	std::map<std::string, std::vector<GffEntry*>* >::iterator oIt;
	oIt = pSortedGffEntries->find( pSeqName->c_str() );

	if (oIt == pSortedGffEntries->end() )
		return NULL;

	return oIt->second;

}

std::vector<std::string>* GffLoader::getSeqNames()
{
	std::vector<std::string>* pSeqNames = new std::vector<std::string>();

	std::map<std::string, std::vector<GffEntry*>* >::iterator oIt;
	for (oIt = pSortedGffEntries->begin(); oIt != pSortedGffEntries->end(); ++oIt)
	{

		pSeqNames->push_back( oIt->first );

	}

	return pSeqNames;

}

std::vector<GffEntry*>* GffLoader::createIntrons( std::vector<GffEntry*>* pTranscriptElements )
{

	/*
	 * at [0] is always the transcript
	 *
	 *
	 */
	if (pTranscriptElements->size() <= 1)
		return NULL;

	std::vector<GffEntry*>* pIntrons = new std::vector<GffEntry*>();

	GffEntry* pTranscript = pTranscriptElements->at(0);

	GffEntry* pOldExon = pTranscriptElements->at(1);

	if (pTranscript->getStart() != pOldExon->getStart())
	{
		// add new intron
		GffEntry* pIntron = new GffEntry(pTranscript, true);

		pIntron->setFeature( new std::string("intron") );

		pIntron->setStart( pTranscript->getStart() );
		pIntron->setEnd( pOldExon->getStart() -1 );

		pIntrons->push_back( pIntron );
	}

	uint32_t iNext = 2;
	GffEntry* pExon = NULL;
	while (iNext < pTranscriptElements->size())
	{
		pExon = pTranscriptElements->at(2);

		uint32_t iDist = pExon->getStart() - pOldExon->getEnd();

		if (iDist > 1)
		{
			GffEntry* pIntron = new GffEntry(pTranscript, true);

			pIntron->setFeature( new std::string("intron") );

			pIntron->setStart( pOldExon->getEnd() +1 );
			pIntron->setEnd( pExon->getStart()-1 );

			pIntrons->push_back( pIntron );
		}

		pOldExon = pExon;

	}

	if (pOldExon->getEnd() != pTranscript->getEnd())
	{
		GffEntry* pIntron = new GffEntry(pTranscript, true);

		pIntron->setFeature( new std::string("intron") );

		pIntron->setStart( pOldExon->getEnd() +1 );
		pIntron->setEnd( pTranscript->getEnd() );

		pIntrons->push_back( pIntron );
	}

	return pIntrons;

}


GffLoader::~GffLoader() {

	// TODO for each key in pSortedGffEntries, for each element in vector, delete
	std::map<std::string, std::vector<GffEntry*>* >::iterator oIt;

	for (oIt = pSortedGffEntries->begin(); oIt != pSortedGffEntries->end(); ++oIt)
	{

		std::vector<GffEntry*>* pGffEntries = oIt->second;

		for (size_t i = 0; i < pGffEntries->size(); ++i)
		{
			delete pGffEntries->at(i);
		}


		delete pGffEntries;

	}

	delete pSortedGffEntries;
}

