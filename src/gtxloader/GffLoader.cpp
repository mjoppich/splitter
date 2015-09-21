/*
 * GffLoader.cpp
 *
 *  Created on: Jul 7, 2015
 *      Author: joppich
 */

#include "GffLoader.h"
#include "GffEntry.h"

#include <utils/Utils.h>

#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <malloc.h>
#include <omp.h>

std::vector<GffLoader::sStatisticElement *> *GffLoader::parseStatFile() {

    std::vector<GffLoader::sStatisticElement *> *pStatPairs = NULL;
    if (m_pStats != NULL) {

        std::vector<std::string> *pStats = Utils::readByLine(m_pStats);
        pStatPairs = new std::vector<GffLoader::sStatisticElement *>();

        for (uint32_t i = 0; i < pStats->size(); ++i) {
            GffLoader::sStatisticElement *pElement = new GffLoader::sStatisticElement();

            std::string sLine = pStats->at(i);

            std::vector<std::string> vElems = Utils::split(sLine, '\t');

            pElement->sParent = vElems.at(0);
            pElement->sBase = vElems.at(1);
            pElement->iModifier = atoi(vElems.at(2).c_str());

            pStatPairs->push_back(pElement);

        }

        pStats->clear();
        delete pStats;

    } else {

        pStatPairs = new std::vector<GffLoader::sStatisticElement *>();
        pStatPairs->push_back(new sStatisticElement("chromosome", "gene", 0));
        pStatPairs->push_back(new sStatisticElement("gene", "transcript", 0));
        pStatPairs->push_back(new sStatisticElement("transcript", "exon", 0));
        pStatPairs->push_back(new sStatisticElement("transcript", "exon", 1));

    }

    return pStatPairs;


}

bool GffLoader::checkConfig()
{

    bool m_bError = false;

    m_pGTxFile = this->m_pParser->getArgument("gtx");

    if (m_pGTxFile == NULL) {
        std::cerr << "You must set the -gtx option to the gt(x) file you want to load." << std::endl;
        m_bError = true;
    }

    if (this->m_pParser->isSet("stats") == true)
    {
        // stats file must exist
        m_pStats = this->m_pParser->getArgument("stats");

        std::cerr << "Stats File: " << *m_pStats << std::endl;

    } else {
        m_pStats = NULL;
    }

    if (this->m_pParser->isSet("prefix") == true)
    {
        m_pPrefix = this->m_pParser->getArgument("prefix");
    } else {
        m_pPrefix = new std::string("");
    }

    m_bValidate = this->m_pParser->isSet("validate");


    if (this->m_pParser->isSet("ignorefeatures") == true)
    {
        // TODO this must still be coded
    } else {
        m_pIgnoreFeatures = NULL;
    }


    return !m_bError;

}

uint32_t GffLoader::prepareRun(CLParser *pParser)
{

    bool m_bPassed = this->checkConfig();

    if (!m_bPassed)
        return 1;

    return 0;

}

void GffLoader::loadGTxFile()
{
    std::ifstream oInputStream;
    std::string sLine;

    pSortedGffEntries = new std::map<std::string, std::vector<GffEntry*>* >();
    m_pChromosomeNames = new std::vector<std::string>();
    m_pChromosomes = new std::vector<GffEntry*>();


    std::map<std::string, std::vector<GffEntry*>* >::iterator oIt;
    oInputStream.open(m_pGTxFile->c_str());

#pragma omp parallel
    {

#pragma omp single
        {
            std::vector<GffEntry *> *pGenes = new std::vector<GffEntry *>();

            std::cout << "Start read in for: " << *m_pGTxFile << " on " << omp_get_max_threads() << " threads" <<
            std::endl;
            while (!oInputStream.eof()) {

                std::getline(oInputStream, sLine);

                if ((sLine[0] == '#') || (sLine[0] == '\n') || (sLine[0] == '\0'))
                    continue;

                GffEntry *pEntry = new GffEntry(&sLine);
                std::string *pSeqName = pEntry->getSeqName();

                if (pSeqName == NULL)
                    continue;

                if (std::find(m_pChromosomeNames->begin(), m_pChromosomeNames->end(), *pSeqName) ==
                    m_pChromosomeNames->end())
                    m_pChromosomeNames->push_back(*pSeqName);

                std::string *pPrefixedSeqName = new std::string(m_pPrefix->c_str());
                pPrefixedSeqName->append(pSeqName->c_str());
                pEntry->setSeqName(pPrefixedSeqName);

                std::string *pFeature = pEntry->getFeature();

                if (m_pIgnoreFeatures != NULL) {
                    std::vector<std::string>::iterator oSIt = std::find(m_pIgnoreFeatures->begin(),
                                                                        m_pIgnoreFeatures->end(), *pFeature);

                    // gff entry feature is in ignore features list
                    if (oSIt != m_pIgnoreFeatures->end()) {
                        delete pEntry;
                        continue;
                    }

                }

                if (pFeature->compare("gene") == 0) {

                    // if size of intermediate vector too large -> add stuff in parallel

                    if (pGenes->size() > 50) {
                        pGenes = this->addGenesInParallel(pGenes);
                    }
                }

                // add entry to current gene
                pGenes->push_back(pEntry);

            }

            pGenes = this->addGenesInParallel(pGenes);

            delete pGenes;

            // now sort the single vectors
            // TODO this is nicely parallelisable
            for (oIt = pSortedGffEntries->begin(); oIt != pSortedGffEntries->end(); ++oIt) {

#pragma omp task firstprivate(oIt)
                {

                    GffEntry::sSortEntriesAsc oSorter;

                    std::vector<GffEntry*>* pGffEntries = oIt->second; // contains genes
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

void GffLoader::run()
{

    uint32_t iReturn = this->prepareRun( this->m_pParser );

    if (iReturn > 0)
        return;

    this->loadGTxFile();

    if (m_pParser->isSet("stats") == true)
    {

        std::vector<GffLoader::sStatisticElement *> *pStatPairs = this->parseStatFile();

        this->printStatistics(pStatPairs);

        if (pStatPairs != NULL) {

            for (uint32_t i = 0; i < pStatPairs->size(); ++i)
                delete pStatPairs->at(i);

            delete pStatPairs;

        }

    }



}

GffLoader::GffLoader(CLParser* pParser)
        : CLRunnable( pParser )
{

}

GffLoader::GffLoader(std::string sArguments)
        : CLRunnable( new CLParser( sArguments ) )
{


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

    delete m_pGTxFile;
    delete m_pStats;
    delete m_pPrefix;

    // TODO for each key in pSortedGffEntries, for each element in vector, delete
    std::map<std::string, std::vector<GffEntry*>* >::iterator oIt;

    if (pSortedGffEntries != NULL)
    {

        for (oIt = pSortedGffEntries->begin(); oIt != pSortedGffEntries->end(); ++oIt) {

            std::vector<GffEntry*>* pGffEntries = oIt->second;

            for (size_t i = 0; i < pGffEntries->size(); ++i) {
                delete pGffEntries->at(i);
            }


            delete pGffEntries;

        }

        delete pSortedGffEntries;

    }

    if (m_pIgnoreFeatures != NULL)
        delete m_pIgnoreFeatures;
}

std::vector<GffEntry *> *GffLoader::addGenesInParallel(std::vector<GffEntry *> *pGenes) {

#pragma omp task firstprivate(pGenes)
    {

        GffEntry *pCurrentGene = NULL;

        for (uint32_t i = 0; i < pGenes->size(); ++i) {

            GffEntry *pProcEntry = pGenes->at(i);

            if (pProcEntry->getFeature()->compare("gene") == 0) {
                pCurrentGene = pProcEntry;
#pragma omp critical
                {
                    std::vector<GffEntry *> *pGffEntries = createEntriesForSeqName(pProcEntry->getSeqName());
                    pGffEntries->push_back(pProcEntry);
                }

            } else {

                pCurrentGene->addChild(pProcEntry);
                continue;

            }

        }

        delete pGenes;

    }

    return new std::vector<GffEntry *>();

}
