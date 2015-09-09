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

    void printStatistics() {

        std::vector< std::pair<std::string, std::string> > vStatPairs;
        vStatPairs.push_back(std::pair<std::string, std::string>("chromosome", "gene"));
        vStatPairs.push_back(std::pair<std::string, std::string>("gene", "transcript"));
        vStatPairs.push_back(std::pair<std::string, std::string>("transcript", "exon"));
        vStatPairs.push_back(std::pair<std::string, std::string>("transcript", "intron"));

        std::vector< uint32_t* > vStatValues;
        std::vector< uint32_t* > vStatInstances;

        const uint32_t iFieldSize = 2;

        for (uint32_t i = 0; i < vStatPairs.size(); ++i) {

            vStatValues.push_back((uint32_t*) calloc(iFieldSize * vStatPairs.size(), sizeof (uint32_t)));
            vStatInstances.push_back((uint32_t*) calloc(iFieldSize * vStatPairs.size(), sizeof (uint32_t)));

        }

        std::vector<GffEntry*>::iterator oIt = m_pChromosomes->begin();
        uint32_t iChromLength = 0;
        for (; oIt != m_pChromosomes->end(); ++oIt) {

            GffEntry* pChrom = *oIt;
            
            iChromLength += pChrom->getLength();

            for (uint32_t i = 0; i < vStatPairs.size(); ++i) {

                std::string* pInstance = &(vStatPairs.at(i).first);
                std::string* pChild = &(vStatPairs.at(i).second);

                uint32_t* pInstances = new uint32_t();
                uint32_t* pLengthInstances = new uint32_t();
                *pLengthInstances = 0;
                *pInstances = 0;

                uint32_t iFound = pChrom->getCount(pInstance, pChild, pInstances);
                uint32_t iLength = pChrom->getTotalLength(pChild, pLengthInstances);

                vStatValues.at(i)[0] += iFound;
                vStatValues.at(i)[1] += iLength;

                vStatInstances.at(i)[0] += *pInstances;
                vStatInstances.at(i)[1] += *pLengthInstances;

                float fInstanceFactor = 0.0f;
                float fLengthFactor = 0.0f;

                if (*pInstances > 0)
                    fInstanceFactor = (float) iFound / (float) *pInstances;

                if (*pLengthInstances > 0)
                    fLengthFactor = (float) iLength / (float) *pLengthInstances;

                //std::cout << "I found " << iFound << " " << *pChild << " in " << *pInstances << " " << *pInstance << "( " << fInstanceFactor << " ) with an average length of " << fLengthFactor << "bp in " << *pChrom->getSeqName() << std::endl;
                std::cout << *pChrom->getSeqName() << "\t" << *pInstance << "\t" << *pChild << "\t" << *pInstances << "\t" << iFound << "\t" << fInstanceFactor;
                std::cout << "\t" << *pLengthInstances << "\t" << iLength << "\t" << fLengthFactor << std::endl;

                delete pLengthInstances;
                delete pInstances;

            }


        }
        
        
        
        std::cout << "Total" << "\t" << "all" << "\t" << "chromosome" << "\t" << 1 << "\t" << m_pChromosomes->size() << "\t" << m_pChromosomes->size() ;
        std::cout << "\t" << m_pChromosomes->size() << "\t" << iChromLength << "\t" << (float)(((float) iChromLength) / (float) m_pChromosomes->size()) << std::endl;

        for (uint32_t i = 0; i < vStatPairs.size(); ++i) {

            std::string* pInstance = &(vStatPairs.at(i).first);
            std::string* pChild = &(vStatPairs.at(i).second);

            uint32_t* pValues = vStatValues.at(i);
            uint32_t* pInstances = vStatInstances.at(i);

            float fInstanceFactor = 0.0f;
            float fLengthFactor = 0.0f;

            if (pInstances[0] > 0)
                fInstanceFactor = (float) pValues[0] / (float) pInstances[0];

            if (pInstances[1] > 0)
                fLengthFactor = (float) pValues[1] / (float) pInstances[1];

            //std::cout << "Total  " << pValues[0] << " " << *pChild << " in " << pInstances[0] << " " << *pInstance << "( " << fInstanceFactor << " ) with an average length of " << fLengthFactor << "bp in total" << std::endl;

            std::cout << "Total" << "\t" << *pInstance << "\t" << *pChild << "\t" << pInstances[0] << "\t" << pValues[0] << "\t" << fInstanceFactor;
            std::cout << "\t" << pInstances[1] << "\t" << pValues[1] << "\t" << fLengthFactor << std::endl;

            free(pValues);
            free(pInstances);

        }


    }

    void printValidation() {
        // http://mblab.wustl.edu/GTF22.html
    }


protected:

    std::map<std::string, std::vector<GffEntry*>* >* pSortedGffEntries;
    std::vector<std::string>* m_pChromosomeNames;
    std::vector<GffEntry*>* m_pChromosomes;

private:

    std::vector<GffEntry*>* createEntriesForSeqName(std::string* pSeqName);
    std::vector<GffEntry*>* createIntrons(std::vector<GffEntry*>* pTranscriptElements);

    struct {

        bool operator()(GffEntry* pElem1, GffEntry* pElem2) {
            if (pElem1->getStart() < pElem2->getStart())
                return true;

            if (pElem1->getStart() > pElem2->getStart())
                return false;

            return pElem1->getLength() > pElem2->getLength();
        }
    } sSortEntriesAsc;

};

#endif /* GFFLOADER_H_ */
