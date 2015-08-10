/*
 * GffEntry.cpp
 *
 *  Created on: Jul 7, 2015
 *      Author: joppich
 */

#include "GffEntry.h"
#include <algorithm>
#include <string>
#include <cstdlib>
#include <iostream>
#include <limits>

#include "../Utils/Utils.h"

	/*
	 *
    Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'

    seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position of the feature, with sequence numbering starting at 1.
    end - End position of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

	 *
	 *
	 */

GffEntry::GffEntry(std::string* pLine )
	: GenomicRegion(0,0), m_pChildren(NULL), m_pAttributes(NULL), m_pRegions(NULL), m_pTranscripts(NULL)
{

	this->initialize();

	uint32_t iOldFound = 0;
	uint32_t iFound = pLine->find( GffEntry::cSep, 0 );
	m_pSeqName = new std::string( *pLine, iOldFound, iFound - iOldFound );
	this->checkForEmptyValue(m_pSeqName);

	iOldFound = iFound + 1;
	iFound = pLine->find( GffEntry::cSep, iOldFound );
	m_pSource = new std::string( *pLine, iOldFound, iFound - iOldFound );
	this->checkForEmptyValue(m_pSource);

	iOldFound = iFound + 1;
	iFound = pLine->find( GffEntry::cSep, iOldFound );
	m_pFeature = new std::string( *pLine, iOldFound, iFound - iOldFound );
	this->checkForEmptyValue(m_pFeature);


	iOldFound = iFound + 1;
	iFound = pLine->find( GffEntry::cSep, iOldFound );
	m_iStart = std::strtoul( pLine->substr( iOldFound, iFound - iOldFound ).c_str(), NULL, 0); // TODO this might not work

	iOldFound = iFound + 1;
	iFound = pLine->find( GffEntry::cSep, iOldFound );
	m_iEnd = std::strtoul( pLine->substr( iOldFound, iFound - iOldFound ).c_str(), NULL, 0); // TODO this might not work

	iOldFound = iFound + 1;
	iFound = pLine->find( GffEntry::cSep, iOldFound );

	std::string sScore = pLine->substr( iOldFound, iFound - iOldFound );
	if ((sScore.length() == 1) && (sScore[0] == '.'))
		m_fScore = std::numeric_limits<float>::quiet_NaN();
	else
		m_fScore = std::strtof( sScore.c_str(), NULL);

	iOldFound = iFound + 1;
	iFound = pLine->find( GffEntry::cSep, iOldFound );

	if (iFound - iOldFound != 1)
		std::cerr << "Error Parsing GTF v2.0 Line: Strand Information: " << *pLine << std::endl;

	char cStrandInfo = pLine->substr( iOldFound, iFound - iOldFound )[0];

	if (cStrandInfo == '.')
		m_iForwardStrand = -1;
	else
		m_iForwardStrand = (cStrandInfo == '+') ? 1 : 0;

	iOldFound = iFound + 1;
	iFound = pLine->find( GffEntry::cSep, iOldFound );

	if (iFound - iOldFound != 1)
		std::cerr << "Error Parsing GTF v2.0 Line: Codon: " << *pLine << std::endl;

	char cFrameInfo = pLine->substr(iOldFound, iFound - iOldFound)[0];
	if (cFrameInfo == '.')
		m_iFrame = -1;
	else
		m_iFrame = cFrameInfo - 48;

	iOldFound = iFound + 1;
	iFound = pLine->find( '\n', iOldFound );

	std::string sAttributes = pLine->substr(iOldFound, iFound - iOldFound);

	this->parseAttributes(sAttributes);

}

GffEntry::GffEntry(GffEntry* pEntry, bool bCopyAttributes)
	: GenomicRegion(pEntry->m_iStart, pEntry->m_iEnd), m_pChildren(NULL), m_pAttributes(NULL), m_pRegions(NULL), m_pTranscripts(NULL)
{

	this->initialize();

	this->m_pSeqName = new std::string( *(pEntry->m_pSeqName) );
	this->m_pFeature = new std::string( *(pEntry->m_pFeature) );
	this->m_pSource = new std::string( *(pEntry->m_pSource) );

	this->m_iStart = pEntry->m_iStart;
	this->m_iEnd = pEntry->m_iEnd;
	this->m_iForwardStrand = pEntry->m_iForwardStrand;
	this->m_iFrame = pEntry->m_iFrame;

	if (bCopyAttributes)
	{

		std::map<std::string, std::string>::iterator oIt;

		for (oIt = pEntry->m_pAttributes->begin(); oIt != pEntry->m_pAttributes->end(); ++oIt)
		{
			this->m_pAttributes->insert( std::pair<std::string, std::string>(oIt->first, oIt->second) );
		}

	}

}

GffEntry::GffEntry( std::string sSeqName, std::string sSource, std::string sFeature, uint32_t iStart, uint32_t iEnd )
: GenomicRegion(iStart, iEnd), m_pChildren(NULL), m_pAttributes(NULL), m_pRegions(NULL), m_pTranscripts(NULL)
{

	this->initialize();

	m_pSeqName = new std::string(sSeqName);
	m_pSource = new std::string(sSource);
	m_pFeature = new std::string(sFeature);

	m_iFrame = 0;
	m_fScore = 0;
	m_iForwardStrand = 1;

}

void GffEntry::initialize()
{

	m_pAttributes = new std::map<std::string, std::string>();

	m_pChildren = new std::vector<GffEntry*>();

	m_bAutoGenerated = false;
	m_bInTranscriptContained = false;

	m_pRegions = NULL;
	m_pTranscripts = NULL;

}

void GffEntry::setFeature(std::string* pNewFeature)
{
	m_pFeature = pNewFeature;
}

GffEntry::~GffEntry() {

	if (m_pAttributes != NULL)
		delete m_pAttributes;

	if (m_pFeature != NULL)
		delete m_pFeature;

	if (m_pSeqName != NULL)
		delete m_pSeqName;

	if (m_pSource != NULL)
		delete m_pSource;

	delete m_pChildren;

}

void GffEntry::printEntry(uint32_t iDepth)
{

	std::string* pGeneID = this->getAttribute(std::string("gene_id"));

	std::cout << *m_pSeqName << "(" << ((pGeneID != NULL) ? *pGeneID : "NULL") <<")" << " " << *m_pSource << " " << *m_pFeature << "["<< m_iStart << " " << m_iEnd << "]" << std::endl;

	if (iDepth == -1)
		return;

	for (uint32_t i = 0; i < m_pChildren->size(); ++i)
		m_pChildren->at(i)->printEntry(iDepth - 1);



}

bool GffEntry::sortChildren( std::vector< std::pair<std::string, std::string> >* pExpands )
{

	uint32_t iExpandIndex = Utils<std::string, std::string>::find(pExpands, this->m_pFeature, NULL);

	bool bPrint = false;

	if (m_pChildren->size() > 0)
	{
		for (size_t i = 0; i < m_pChildren->size(); ++i)
		{
			m_pChildren->at(i)->sortChildren( pExpands );
		}

		std::vector<GffEntry*>* pRegions = NULL;

		if (iExpandIndex != -1)
		{

			if (m_pChildren->size() > 1)
			{
				std::sort(m_pChildren->begin(), m_pChildren->end(), GffEntry::sSortEntriesAsc);
			}

			pRegions = this->fillEmptyRegion( pExpands->at(iExpandIndex).second );
		}

		if (pRegions != NULL)
		{

			if (pRegions->size() > 1)
			{
				this->printEntry();
				bPrint = true;
			}

			m_pChildren->insert(m_pChildren->end(), pRegions->begin(), pRegions->end());
		}

		if (m_pChildren->size() > 1)
		{
			std::sort(m_pChildren->begin(), m_pChildren->end(), GffEntry::sSortEntriesAsc);
		}

		if (bPrint)
			this->printEntry();
	}


	return true;
}

std::vector<GffEntry*>* GffEntry::fillEmptyRegion(std::string sFeature)
{
	if (m_pChildren->size() <= 1)
	{
		return NULL;
	}

	std::vector<GffEntry*>* pIntrons = new std::vector<GffEntry*>();

	GffEntry* pTranscript = this;

	GffEntry* pOldExon = m_pChildren->at(0);

	if (pTranscript->getStart() != pOldExon->getStart())
	{
		// add new intron
		GffEntry* pIntron = new GffEntry(pTranscript, true);

		pIntron->setFeature( new std::string( sFeature ) );

		pIntron->setSource( new std::string( "splitter(auto)" ));

		pIntron->setStart( pTranscript->getStart() );
		pIntron->setEnd( pOldExon->getStart() -1 );

		pIntron->setAutogenerated(true);

		pIntrons->push_back( pIntron );
	}

	uint32_t iNext = 1;
	GffEntry* pExon = NULL;
	while (iNext < m_pChildren->size())
	{
		pExon = m_pChildren->at( iNext );
		++iNext;

		uint32_t iDist = pExon->getStart() - pOldExon->getEnd();

		if (iDist > 1)
		{
			GffEntry* pIntron = new GffEntry(pTranscript, true);

			pIntron->setFeature( new std::string( sFeature ) );
			pIntron->setSource( new std::string( "splitter(auto)" ));

			pIntron->setStart( pOldExon->getEnd() +1 );
			pIntron->setEnd( pExon->getStart()-1 );

			pIntron->setAutogenerated(true);

			if (pIntron->getStart() <= pIntron->getEnd())
				pIntrons->push_back( pIntron );
		}

		pOldExon = pExon;

	}

	if (pOldExon->getEnd() != pTranscript->getEnd())
	{
		GffEntry* pIntron = new GffEntry(pTranscript, true);

		pIntron->setFeature( new std::string( sFeature ) );
		pIntron->setSource( new std::string( "splitter(auto)" ));

		pIntron->setStart( pOldExon->getEnd() +1 );
		pIntron->setEnd( pTranscript->getEnd() );

		pIntron->setAutogenerated(true);

		if (pIntron->getStart() <= pIntron->getEnd())
			pIntrons->push_back( pIntron );
	}

	return pIntrons;
}

void GffEntry::parseAttributes(std::string& sAttributes)
{

	size_t iFound = 0;
	size_t iOldFound = 0;

	std::replace(sAttributes.begin(), sAttributes.end(), '\t', ' ');

	bool bVisitedLoop = false;

	while ( ( iFound = sAttributes.find(';', iOldFound) ) != std::string::npos )
	{

		bVisitedLoop = true;

		std::string sAttribute = sAttributes.substr(iOldFound, iFound - iOldFound);
		sAttribute = sAttribute.substr( sAttribute.find_first_not_of(' '), -1 );

		this->parseAttribute(sAttribute);


		iOldFound = iFound + 1;

	}

	if (!bVisitedLoop)
	{
		if (sAttributes.length() > 0)
			this->parseAttribute(sAttributes);
	}

}

void GffEntry::parseAttribute(std::string& sAttribute)
{

	std::string sKey, sValue;
	size_t iSplit = 0;

	iSplit = sAttribute.find("=");

	if (iSplit != std::string::npos)
	{

		sKey = sAttribute.substr(0, iSplit);
		sValue = sAttribute.substr(iSplit+1, -1);

	} else {
		iSplit = sAttribute.find(" ");

		if (iSplit != std::string::npos)
		{
			sKey = sAttribute.substr(0, iSplit);
			sValue = sAttribute.substr(iSplit+1, -1);

			if (sValue[0] == '\"')
			{
				sValue = sValue.substr(1, sValue.length()-2);
			}

		}
	}

	m_pAttributes->insert( std::pair<std::string, std::string>(sKey, sValue) );

}

void GffEntry::setInTranscriptContained(bool bValue)
{
	m_bInTranscriptContained = bValue;
}

bool GffEntry::getInTranscriptContained()
{
	return m_bInTranscriptContained;
}

void GffEntry::flatten(std::string* pFlattenLevel)
{

	std::cout << this << std::endl;

	std::cout << this->m_pFeature->compare( pFlattenLevel->c_str() ) << std::endl;

	if (m_pFeature->compare( pFlattenLevel->c_str() ) != 0)
	{
		for (uint32_t i = 0; i < m_pChildren->size(); ++i)
		{
			GffEntry* pElem = m_pChildren->at(i);
			pElem->flatten(pFlattenLevel);
		}

		return;
	}

	GffEntry* pAllRegions = new GffEntry( *this->getSeqName(), *this->getSource(), "region", this->getStart(), this->getEnd());
	std::vector<GffEntry*>* pRegions = new std::vector<GffEntry*>();
	pRegions->push_back(pAllRegions);


	for (uint32_t i = 0; i < m_pChildren->size(); ++i)
	{

		GffEntry* pElement = m_pChildren->at(i);

		std::vector<GffEntry*>::iterator oIt = pRegions->begin();
		for ( ; oIt != pRegions->end(); ++oIt)
		{

			GffEntry* pRegion = *oIt;;

			std::vector<GffEntry*>* pNewRegions = pRegion->split(pElement->getStart());

			// problem: if split needed, end must also be split

			if (pNewRegions != NULL)
			{
				std::vector<GffEntry*>::iterator oJt = pNewRegions->begin();
				for (; oJt != pNewRegions->end(); ++oJt)
				{

					GffEntry* pNewElement = *oJt;

					std::vector<GffEntry*>* pDoubleSplit = pNewElement->split(pElement->getEnd());

					if (pDoubleSplit != NULL)
					{
						pNewRegions->erase(oJt);
						pNewRegions->insert(pNewRegions->end(), pDoubleSplit->begin(), pDoubleSplit->end());
						break;
					}
				}

			} else {

				pNewRegions = pRegion->split(pElement->getEnd());

			}

			// pnewregions is null!!

			if (pNewRegions != NULL)
			{
				pRegions->erase(oIt);
				pRegions->insert(pRegions->end(), pNewRegions->begin(), pNewRegions->end());
			}

		}

	}

	std::sort(pRegions->begin(), pRegions->end(), GffEntry::sSortEntriesAsc);

	std::vector<GffTranscript*>* pTranscripts = new std::vector<GffTranscript*>();

	std::vector<GffEntry*>::iterator oIt;
	for (oIt = m_pChildren->begin(); oIt != m_pChildren->end(); ++oIt)
	{

		GffEntry* pElem = *oIt;
		std::string* pTID = pElem->getAttribute("transcript_id");

		if (pTID == NULL )
			continue;

                GffTranscript* pTranscript = new GffTranscript( *pTID, pElem->getStart(), pElem->getEnd() );

		std::vector<GffEntry*>::iterator oJt;
		for ( oJt = pElem->getChildren()->begin(); oJt != pElem->getChildren()->end(); ++oJt)
		{

			GffEntry* pElem = *oJt;

                        
			std::vector<GffEntry*>* pExons = this->find(pRegions, pElem->getStart(), pElem->getEnd());
                        
                        if (pExons->size() != 1)
                        {
                            std::cerr << "Problem with Transcript ID: " << *pTID << std::endl;
                        }
                        

			for (uint32_t i = 0; i < pExons->size(); ++i)
			{
				GffEntry* pElem = pExons->at(i);
				pElem->setInTranscriptContained( true );
			}

			pTranscript->addExons( pExons );
                        
		}
                
                pTranscripts->push_back(pTranscript);


	}
        
        if ((pRegions->size() == 0) || (pTranscripts->size() == 0))
        {
            std::cout << "problem " << pRegions->size() << " " << pTranscripts->size() << std::endl;
        }

	m_pRegions = pRegions;
	m_pTranscripts = pTranscripts;

}

