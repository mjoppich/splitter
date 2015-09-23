/*
 * GffLoader.h
 *
 *  Created on: Jul 7, 2015
 *      Author: joppich
 */

#ifndef GFFLOADER_H_
#define GFFLOADER_H_

#include <map>
#include <vector>
#include <inttypes.h>
#include <string>
#include <iostream>
#include <utils/CLRunnable.h>

class GffEntry;

class GffLoader : public CLRunnable {
public:

    GffLoader(CLParser* pParser);
    GffLoader(std::string sArguments);
    virtual ~GffLoader();

void run();

    struct sStatisticElement {
        std::string sParent;
        std::string sBase;
        uint8_t iModifier;
        
        sStatisticElement()
        {
            
        }

        sStatisticElement(std::string sParent, std::string sBase, uint8_t iModifier) {
            this->sParent = sParent;
            this->sBase = sBase;
            this->iModifier = iModifier;
        }

    };

    struct sStatisticResult {

        uint32_t iLengthBase;
        uint32_t iBaseCount;
        uint32_t iParentCount;
        
        uint32_t iModifier;
        
        float fAvgLength;
        float fAvgCount;
        
        void prepareResults()
        {
            fAvgCount = (float) iBaseCount / (float) iParentCount;
            fAvgLength = (float) iLengthBase / (float) iBaseCount;
            
            if (iParentCount == 0)
                fAvgCount = 0;
            
            if (iBaseCount == 0)
                fAvgLength = 0;
        }
        
        
        void addResult( sStatisticResult* pOther)
        {
            
            iLengthBase = iLengthBase + pOther->iLengthBase;
            iBaseCount = iBaseCount + pOther->iBaseCount;
            iParentCount = iParentCount + pOther->iParentCount;
            
            this->prepareResults();
            
        }
    };

    std::vector<GffEntry *> *addGenesInParallel(std::vector<GffEntry *> *pGenes);
    
    void printStatisticResult(std::string sPrefix, sStatisticElement* pElement, sStatisticResult* pResult)
    {
        char cDel = '\t';

        pResult->prepareResults();
        
        if (sPrefix.size() > 0 )
            std::cout << sPrefix << cDel;
        
        std::cout << pElement->sParent << cDel << pElement->sBase << cDel << (char)(48+pElement->iModifier) << cDel << pResult->iParentCount << cDel << pResult->iBaseCount << cDel << pResult->fAvgCount << cDel << pResult->iLengthBase << cDel << pResult->fAvgLength << std::endl;
    }

    std::vector<GffEntry*>* getEntriesForSeqName(std::string* pSeqName);

    GffEntry* getChromosome(std::string* pSeqName);
    std::vector<std::string>* getSeqNames();

    void printStatistics(std::vector< sStatisticElement* >* pStatPairs);

    void printValidation() {
        // http://mblab.wustl.edu/GTF22.html
    }

protected:

    void loadLines(std::vector<std::string>* pLines);

    std::vector<GffLoader::sStatisticElement *> *parseStatFile();

    std::map<std::string, std::vector<GffEntry*>* >* pSortedGffEntries = NULL;
    std::vector<GffEntry*>* m_pChromosomes;
    std::vector<std::string>* m_pIgnoreFeatures;
    std::vector<std::string>* m_pChromosomeNames;

private:

    bool checkConfig();
    uint32_t prepareRun(CLParser* pParser);
    void loadGTxFile();

    std::string* m_pGTxFile;
    std::string* m_pStats;
    bool m_bValidate;
    std::string* m_pPrefix;

    std::vector<GffEntry*>* createEntriesForSeqName(std::string* pSeqName);
    std::vector<GffEntry*>* createIntrons(std::vector<GffEntry*>* pTranscriptElements);

    /*
    struct {

        bool operator()(GffEntry* pElem1, GffEntry* pElem2) {
            if (pElem1->getStart() < pElem2->getStart())
                return true;

            if (pElem1->getStart() > pElem2->getStart())
                return false;

            return pElem1->getLength() > pElem2->getLength();
        }
    } sSortEntriesAscASD;
    
    */

};

#endif /* GFFLOADER_H_ */

