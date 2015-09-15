/*
 * GffLoader.cpp
 *
 *  Created on: Jul 7, 2015
 *      Author: joppich
 */

#include "GffLoader.h"
#include "GffEntry.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <malloc.h>
#include <omp.h>

GffLoader::GffLoader(std::string& sFileName, std::vector<std::string>* pIgnoreFeatures, std::string* pPrefix) {

    std::ifstream oInputStream;
    std::string sLine;

    pSortedGffEntries = new std::map<std::string, std::vector<GffEntry*>* >();
    m_pChromosomeNames = new std::vector<std::string>();
    m_pChromosomes = new std::vector<GffEntry*>();


    std::map<std::string, std::vector<GffEntry*>* >::iterator oIt;

    std::vector<GffEntry*>* pCurrentTranscript = NULL;

    oInputStream.open(sFileName.c_str());

    GffEntry* pCurrentGene = NULL;

    std::cout << "Start read in for: " << sFileName << " on " << omp_get_max_threads() << " threads" << std::endl;
    while (!oInputStream.eof()) {

        std::getline(oInputStream, sLine);

        if ((sLine[0] == '#') || (sLine[0] == '\n') || (sLine[0] == '\0'))
            continue;

        GffEntry* pEntry = new GffEntry(&sLine);
        std::string* pSeqName = pEntry->getSeqName();

        if (pSeqName == NULL)
            continue;

        if (std::find(m_pChromosomeNames->begin(), m_pChromosomeNames->end(), *pSeqName) == m_pChromosomeNames->end())
            m_pChromosomeNames->push_back(*pSeqName);

        std::string* pPrefixedSeqName = new std::string(pPrefix->c_str());
        pPrefixedSeqName->append(pSeqName->c_str());
        pEntry->setSeqName(pPrefixedSeqName);
        pSeqName = pEntry->getSeqName();

        std::string* pFeature = pEntry->getFeature();

        if (pIgnoreFeatures != NULL) {
            std::vector<std::string>::iterator oSIt = std::find(pIgnoreFeatures->begin(), pIgnoreFeatures->end(), *pFeature);

            // gff entry feature is in ignore features list
            if (oSIt != pIgnoreFeatures->end()) {
                delete pEntry;
                continue;
            }

        }


        if (pFeature->compare("gene") == 0) {
            pCurrentGene = pEntry;

            std::vector<GffEntry*>* pGffEntries = createEntriesForSeqName(pSeqName);
            pGffEntries->push_back(pEntry);

        } else {

            pCurrentGene->addChild(pEntry);
            continue;

        }




    }

#pragma omp parallel
    {

#pragma omp single
        {

            // now sort the single vectors
            // TODO this is nicely parallelisable
            for (oIt = pSortedGffEntries->begin(); oIt != pSortedGffEntries->end(); ++oIt) {

#pragma omp task firstprivate(oIt)
                {

                    GffEntry::sSortEntriesAsc oSorter;

                    std::vector<GffEntry*>* pGffEntries = oIt->second; // contains genes

                    for (size_t i = 0; i < pGffEntries->size(); ++i) {
                        GffEntry* pGeneEntry = pGffEntries->at(i);
                    }
                    std::sort(pGffEntries->begin(), pGffEntries->end(), oSorter);
                }

            }

        }

    }

    std::string* pFlattenLevel = new std::string("gene");

    std::vector<std::string>::iterator oJt;


    uint32_t i = 0;
#pragma omp parallel for schedule(dynamic,1)
    for (i = 0; i < m_pChromosomeNames->size(); ++i) {
        std::string sChromName = m_pChromosomeNames->at(i);

        std::vector<GffEntry*>* pGffEntries = this->getEntriesForSeqName(&sChromName);

        GffEntry* pChromosome = new GffEntry(sChromName, "splitter", "chromosome", 1, 1);
        pChromosome->addChildren(pGffEntries);
        pChromosome->sortChildren(NULL);

        uint32_t iEnd = pChromosome->getMaxLocation();
        pChromosome->setEnd(iEnd);

        pChromosome->flatten(pFlattenLevel);

#pragma omp critical
        {
            m_pChromosomes->push_back(pChromosome);
        }

    }

    delete pFlattenLevel;

}

std::vector<GffEntry*>* GffLoader::createEntriesForSeqName(std::string* pSeqName) {
    std::map<std::string, std::vector<GffEntry*>* >::iterator oIt;
    oIt = pSortedGffEntries->find(pSeqName->c_str());

    if (oIt == pSortedGffEntries->end()) {
        // needs to be inserted
        std::vector<GffEntry*>* pGffEntries = new std::vector<GffEntry*>();

        pSortedGffEntries->insert(std::pair<std::string, std::vector<GffEntry*>* >(*pSeqName, pGffEntries));

        return pGffEntries;
    }

    return oIt->second;

}

std::vector<GffEntry*>* GffLoader::getEntriesForSeqName(std::string* pSeqName) {
    std::map<std::string, std::vector<GffEntry*>* >::iterator oIt;
    oIt = pSortedGffEntries->find(pSeqName->c_str());

    if (oIt == pSortedGffEntries->end())
        return NULL;

    return oIt->second;

}

GffEntry* GffLoader::getChromosome(std::string* pSeqName) {
    std::vector< GffEntry* > ::iterator oIt;
    for (oIt = m_pChromosomes->begin(); oIt != m_pChromosomes->end(); ++oIt) {
        GffEntry* pChrom = *oIt;

        if (pChrom->getSeqName()->compare(*pSeqName) == 0)
            return pChrom;

    }

    return NULL;

}

std::vector<std::string>* GffLoader::getSeqNames() {
    std::vector<std::string>* pSeqNames = new std::vector<std::string>();

    std::map<std::string, std::vector<GffEntry*>* >::iterator oIt;
    for (oIt = pSortedGffEntries->begin(); oIt != pSortedGffEntries->end(); ++oIt) {

        pSeqNames->push_back(oIt->first);

    }

    return pSeqNames;

}

void GffLoader::printStatistics(std::vector<GffLoader::sStatisticElement*>* pStatPairs) {

    if (pStatPairs == NULL) {

        pStatPairs = new std::vector< GffLoader::sStatisticElement* >();
        pStatPairs->push_back(new sStatisticElement("chromosome", "gene", 0));
        pStatPairs->push_back(new sStatisticElement("gene", "transcript", 0));
        pStatPairs->push_back(new sStatisticElement("transcript", "exon", 0));
        pStatPairs->push_back(new sStatisticElement("transcript", "exon", 1));

    }

    std::vector< uint32_t* > vStatValues;
    std::vector< uint32_t* > vStatInstances;

    uint32_t iResultSets = m_pChromosomes->size() * pStatPairs->size();
    sStatisticResult** pResults = (sStatisticResult**) calloc(sizeof (sStatisticResult*), iResultSets);
    sStatisticResult** pTotalResult = (sStatisticResult**) calloc(sizeof (sStatisticResult*), pStatPairs->size());

    size_t c = 0;
    uint32_t iMax = m_pChromosomes->size();

#pragma omp parallel for schedule(dynamic,1) shared(iMax, pResults, pStatPairs)
    for (c = 0; c < iMax; ++c) {

        m_pChromosomes->at(c)->printHierarchy(0);

        for (uint32_t i = 0; i < pStatPairs->size(); ++i) {
            GffEntry* pChrom = m_pChromosomes->at(c);
            
            
            
            uint32_t iIndex = c * pStatPairs->size() + i;

            sStatisticElement* pElement = pStatPairs->at(i);

            pResults[iIndex] = new sStatisticResult();

            sStatisticResult* pResult = pResults[iIndex];
            pChrom->getStatistics(pElement, pResult);
            pResult->prepareResults();

        }
    }

    for (c = 0; c < iMax; ++c) {
        GffEntry* pChrom = m_pChromosomes->at(c);

        for (uint32_t i = 0; i < pStatPairs->size(); ++i) {

            uint32_t iIndex = c * pStatPairs->size() + i;
            sStatisticElement* pElement = pStatPairs->at(i);
            sStatisticResult* pResult = pResults[iIndex];

            this->printStatisticResult(*(pChrom->getSeqName()), pElement, pResult);
        }

    }

    uint32_t iIndex = 0;

    for (uint32_t i = 0; i < pStatPairs->size(); ++i) {
        iIndex = i;
        sStatisticElement* pElement = pStatPairs->at(i);
        pTotalResult[i] = new sStatisticResult();

        sStatisticResult* pResult = pTotalResult[i];

        while (iIndex < iResultSets) {
            pResult->addResult(pResults[iIndex]);

            iIndex += pStatPairs->size();
        }

        this->printStatisticResult("total", pElement, pResult);

    }


    free(pResults);
    free(pTotalResult);

}

std::vector<GffEntry*>* GffLoader::createIntrons(std::vector<GffEntry*>* pTranscriptElements) {

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

    if (pTranscript->getStart() != pOldExon->getStart()) {
        // add new intron
        GffEntry* pIntron = new GffEntry(pTranscript, true);

        pIntron->setFeature(new std::string("intron"));

        pIntron->setStart(pTranscript->getStart());
        pIntron->setEnd(pOldExon->getStart() - 1);

        pIntrons->push_back(pIntron);
    }

    uint32_t iNext = 2;
    GffEntry* pExon = NULL;
    while (iNext < pTranscriptElements->size()) {
        pExon = pTranscriptElements->at(2);

        uint32_t iDist = pExon->getStart() - pOldExon->getEnd();

        if (iDist > 1) {
            GffEntry* pIntron = new GffEntry(pTranscript, true);

            pIntron->setFeature(new std::string("intron"));

            pIntron->setStart(pOldExon->getEnd() + 1);
            pIntron->setEnd(pExon->getStart() - 1);

            pIntrons->push_back(pIntron);
        }

        pOldExon = pExon;

    }

    if (pOldExon->getEnd() != pTranscript->getEnd()) {
        GffEntry* pIntron = new GffEntry(pTranscript, true);

        pIntron->setFeature(new std::string("intron"));

        pIntron->setStart(pOldExon->getEnd() + 1);
        pIntron->setEnd(pTranscript->getEnd());

        pIntrons->push_back(pIntron);
    }

    return pIntrons;

}

GffLoader::~GffLoader() {

    // TODO for each key in pSortedGffEntries, for each element in vector, delete
    std::map<std::string, std::vector<GffEntry*>* >::iterator oIt;

    for (oIt = pSortedGffEntries->begin(); oIt != pSortedGffEntries->end(); ++oIt) {

        std::vector<GffEntry*>* pGffEntries = oIt->second;

        for (size_t i = 0; i < pGffEntries->size(); ++i) {
            delete pGffEntries->at(i);
        }


        delete pGffEntries;

    }

    delete pSortedGffEntries;
}

