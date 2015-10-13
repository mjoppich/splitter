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

GffLoader::GffLoader(CLParser* pParser)
        : CLRunnable( pParser )
{

}

GffLoader::GffLoader(std::string sArguments)
        : CLRunnable( new CLParser( sArguments ) )
{


}

std::vector<GffLoader::sStatisticElement *> *GffLoader::parseStatFile(std::string* pStatFileLocation) {

    std::vector<GffLoader::sStatisticElement *> *pStatPairs = NULL;
    if (pStatFileLocation != NULL) {

        std::vector<std::string> *pStats = Utils::readByLine(pStatFileLocation);
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

    std::ios_base::sync_with_stdio (false);

    std::ifstream oInputStream;
    std::string sLine;


    pSortedGffEntries = new std::map<std::string, std::vector<GffEntry*>* >();
    m_pChromosomeNames = new std::vector<std::string>();
    m_pChromosomes = new std::vector<GffEntry*>();

    std::map<std::string, std::vector<GffEntry*>* >::iterator oIt;
    oInputStream.open(m_pGTxFile->c_str());

    if (oInputStream.eof())
        return;

#pragma omp parallel
    {

#pragma omp single
        {
            std::vector<GffEntry *> *pGenes = new std::vector<GffEntry *>();

            std::cout << "Start read in for: " << *m_pGTxFile << " on " << omp_get_max_threads() << " threads" <<
            std::endl;

            std::vector< std::string >* pLinePackage = new std::vector<std::string>();

            while (std::getline(oInputStream, sLine) ) {

                if ((sLine[0] == '#') || (sLine[0] == '\n') || (sLine[0] == '\0'))
                    continue;

                /*
                 *
                 * does it suffice to search for '\tgene\t' ?
                 *
                 *
                 */

                if (sLine.find( "\tgene\t" ) != std::string::npos )
                {

                    if (pLinePackage->size() > 200)
                    {

#pragma omp task untied firstprivate(pLinePackage)
                        {
                            this->loadLines( pLinePackage );

                            delete pLinePackage;
                        }

                        pLinePackage = new std::vector<std::string>();
                        pLinePackage->reserve(250);
                    }
                }

                pLinePackage->push_back( sLine );

                continue;

            }
#pragma omp task untied firstprivate(pLinePackage)
            {
                this->loadLines( pLinePackage );

                delete pLinePackage;

            }

#pragma omp taskwait


            // now sort the single vectors
            for (oIt = pSortedGffEntries->begin(); oIt != pSortedGffEntries->end(); ++oIt) {

#pragma omp task untied firstprivate(oIt)
                {

                    std::string sChromName = oIt->first;
                    std::vector<GffEntry*>* pGffEntries = oIt->second; // contains genes
                    GffEntry* pChromosome = new GffEntry(sChromName, "splitter", "chromosome", 1, 1);
                    pChromosome->addChildren(pGffEntries);
                    pChromosome->sortChildren(NULL);

                    uint32_t iEnd = pChromosome->getMaxLocation();
                    pChromosome->setEnd(iEnd);

#pragma omp critical
                    {
                        m_pChromosomes->push_back(pChromosome);
                    }
                }

            }

            //delete pFlattenLevel;


        }

    }

}

void GffLoader::run()
{

    uint32_t iReturn = this->prepareRun( this->m_pParser );

    if (iReturn > 0)
        return;

    this->loadGTxFile();

    if (m_pParser->isSet("stats") == true)
    {

        std::vector<GffLoader::sStatisticElement *> *pStatPairs = this->parseStatFile( m_pStats );

        this->printStatistics(pStatPairs);

        if (pStatPairs != NULL) {

            for (uint32_t i = 0; i < pStatPairs->size(); ++i)
                delete pStatPairs->at(i);

            delete pStatPairs;

        }

    }



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

    uint32_t iResultSets = m_pChromosomes->size() * pStatPairs->size();
    sStatisticResult** pResults = (sStatisticResult**) calloc(sizeof (sStatisticResult*), iResultSets);
    sStatisticResult** pTotalResult = (sStatisticResult**) calloc(sizeof (sStatisticResult*), pStatPairs->size());

    size_t c = 0;
    uint32_t iMax = m_pChromosomes->size();
    uint32_t iWorkPackageSize = 0;
    uint32_t iWorkPackageStart = 0;
    uint32_t iWorkPackageEnd = 0;

#pragma omp parallel
    {
#pragma omp single
        {
            for (c = 0; c < iMax; ++c) {

                //m_pChromosomes->at(c)->printHierarchy(0);
                GffEntry* pChrom = m_pChromosomes->at(c);
                uint32_t iChildren = pChrom->getChildren()->size();
                iWorkPackageSize += iChildren;

                if (iWorkPackageSize > 200)
                {

                    iWorkPackageEnd = c;

#pragma omp task firstprivate(iWorkPackageStart, iWorkPackageEnd, pStatPairs, pResults)
                    {
                        for (uint32_t c = iWorkPackageStart; c < iWorkPackageEnd+1; ++c)
                        {
                            GffEntry* pChrom = m_pChromosomes->at(c);

                            for (uint32_t i = 0; i < pStatPairs->size(); ++i) {

                                uint32_t iIndex = c * pStatPairs->size() + i;

                                sStatisticElement* pElement = pStatPairs->at(i);

                                pResults[iIndex] = new sStatisticResult();

                                sStatisticResult* pResult = pResults[iIndex];
                                pChrom->getStatistics(pElement, pResult);
                                pResult->prepareResults();

                            }
                        }

                    }

                    iWorkPackageSize = 0;
                    iWorkPackageStart = c+1;
                }


            }

#pragma omp task firstprivate(iWorkPackageStart, iMax, pStatPairs, pResults)
            {
                for (uint32_t c = iWorkPackageStart; c < iMax; ++c)
                {
                    GffEntry* pChrom = m_pChromosomes->at(c);

                    for (uint32_t i = 0; i < pStatPairs->size(); ++i) {

                        uint32_t iIndex = c * pStatPairs->size() + i;

                        sStatisticElement* pElement = pStatPairs->at(i);

                        pResults[iIndex] = new sStatisticResult();

                        sStatisticResult* pResult = pResults[iIndex];
                        pChrom->getStatistics(pElement, pResult);
                        pResult->prepareResults();

                    }
                }

            }
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

void GffLoader::loadLines(std::vector<std::string>* pLines) {

    std::string *pSeqName;
    std::string *pFeature;
    std::string *pPrefixedSeqName;

    GffEntry* pCurrentGene = NULL;
    std::vector< GffEntry* > vGenes;

    for (uint32_t i = 0; i < pLines->size(); ++i)
    {

        std::string sLine = pLines->at(i);

        GffEntry *pEntry = new GffEntry(&sLine);
        pSeqName = pEntry->getSeqName();

        if (pSeqName == NULL)
        {
            delete pEntry;
            continue;
        }

        pFeature = pEntry->getFeature();

        if (m_pIgnoreFeatures != NULL) {
            bool bIsIgnoreFeature = (m_pIgnoreFeatures->end() != std::find(m_pIgnoreFeatures->begin(), m_pIgnoreFeatures->end(), *pFeature));

            // gff entry feature is in ignore features list
            if ( bIsIgnoreFeature ) {
                delete pEntry;
                continue;
            }

        }

        pPrefixedSeqName = new std::string(m_pPrefix->c_str());
        pPrefixedSeqName->append(pSeqName->c_str());
        pEntry->setSeqName(pPrefixedSeqName);

#pragma omp critical
        {
            if (std::find(m_pChromosomeNames->begin(), m_pChromosomeNames->end(), *(pEntry->getSeqName())) == m_pChromosomeNames->end())
            {
                m_pChromosomeNames->push_back( *(pEntry->getSeqName()) );
            }
        }



        if (pEntry->getFeature()->compare("gene") == 0) {
            pCurrentGene = pEntry;
            vGenes.push_back(pCurrentGene);
        } else {

            pCurrentGene->addChild(pEntry);
            continue;

        }

    }


    std::string* pFlattenLevel = new std::string("gene");
    std::string* pChildLevel = new std::string("exon");

    for (uint32_t i = 0; i < vGenes.size(); ++i)
    {

        GffEntry* pProcEntry = vGenes.at(i);

        pProcEntry->flatten(pFlattenLevel, pChildLevel);

#pragma omp taskyield
#pragma omp critical
        {
            std::vector<GffEntry *> *pGffEntries = createEntriesForSeqName(pProcEntry->getSeqName());
            pGffEntries->push_back(pProcEntry);
        }

    }

    delete pFlattenLevel;
    delete pChildLevel;


}

void GffLoader::printStatisticResult(std::string sPrefix, sStatisticElement* pElement, sStatisticResult* pResult) {
    char cDel = '\t';

    pResult->prepareResults();

    if (sPrefix.size() > 0 )
        std::cout << sPrefix << cDel;

    std::cout << pElement->sParent << cDel << pElement->sBase << cDel << (char)(48+pElement->iModifier) << cDel << pResult->iParentCount << cDel << pResult->iBaseCount << cDel << pResult->fAvgCount << cDel << pResult->iLengthBase << cDel << pResult->fAvgLength << std::endl;
}