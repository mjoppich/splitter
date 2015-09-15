/*
 * GffEntry.h
 *
 *  Created on: Jul 7, 2015
 *      Author: joppich
 */

#ifndef GFFENTRY_H_
#define GFFENTRY_H_

#include <vector>
#include <inttypes.h>
#include <map>
#include <algorithm>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "GenomicRegion.h"
#include "GffLoader.h"

class GffTranscript;

class GffEntry : public GenomicRegion {
public:

    GffEntry(std::string* pLine);
    GffEntry(std::string* pLine, uint8_t cSep);
    GffEntry(GffEntry* pEntry, bool bCopyAttributes);
    GffEntry(std::string sSeqName, std::string sSource, std::string sFeature, uint32_t iStart, uint32_t iEnd);

    struct sSortEntriesAsc {

        bool operator()(GffEntry* pElem1, GffEntry* pElem2) {
            if (pElem1->getStart() < pElem2->getStart())
                return true;

            if (pElem1->getStart() > pElem2->getStart())
                return false;

            return pElem1->getLength() > pElem2->getLength();
        }
    } ;

    virtual ~GffEntry();

    void parseLine(std::string* pLine, uint8_t cSep);

    //static const uint8_t cSep = '\t';

    void setAutogenerated(bool bAuto) {
        m_bAutoGenerated = bAuto;
    }

    bool getAutogenerated() {
        return m_bAutoGenerated;
    }

    std::string* getSeqName() {
        return m_pSeqName;
    }

    void setSeqName(std::string* pNewSeqName) {
        if (m_pSeqName != NULL)
            delete m_pSeqName;

        m_pSeqName = pNewSeqName;
    }

    std::string* getSource() {
        return m_pSource;
    }

    std::string* getFeature() {
        return this->m_pFeature;
    }

    void setFeature(std::string* pNewFeature);

    void setSource(std::string* pNewSource) {
        m_pSource = pNewSource;
    }

    void setAttribute(std::string& sKey, std::string sValue) {
        m_pAttributes->insert(std::pair<std::string, std::string>(sKey, sValue));
    }

    std::string* getAttribute(std::string sKey) {
        std::map<std::string, std::string>::iterator oIt = m_pAttributes->find(sKey);

        if (oIt != m_pAttributes->end())
            return new std::string(oIt->second);

        return NULL;
    }

    GffEntry* getParent() {
        return m_pParent;
    }

    void setParent(GffEntry* pParent) {
        m_pParent = pParent;
    }

    std::vector<GffTranscript*>* getTranscripts() {
        return m_pTranscripts;
    }

    uint32_t charToStrand(char cStrandInfo) {
        uint8_t iStrandInfo;
        if (cStrandInfo == '.')
            iStrandInfo = -1;
        else
            iStrandInfo = (cStrandInfo == '+') ? 1 : 0;

        return iStrandInfo;
    }

    // flattens entry at when feature == flattenlevel
    void flatten(std::string* pFlattenLevel);

    GffEntry* getRegion(uint32_t iPosition) {

        std::vector<GffEntry*>::iterator oIt;
        for (oIt = m_pRegions->begin(); oIt != m_pRegions->end(); ++oIt) {

            GffEntry* pRegion = *oIt;

            if (pRegion->contains(iPosition) == true)
                return pRegion;

        }

        return NULL;

    }

    std::vector<GffTranscript*>* hasTranscript(std::vector<uint32_t>* pPositions, bool bHasPartialContainment);

    std::vector<GffEntry*>* findLevels(std::vector<uint32_t>* pPositions, std::string* pLevel, bool bPartialContainment = false) {

        std::vector<GffEntry*>* pLevelContained = new std::vector<GffEntry*>();

        if (this->m_pFeature->compare(*pLevel) != 0) {

            for (uint32_t i = 0; i < m_pChildren->size(); ++i) {

                GffEntry* pChild = m_pChildren->at(i);

                std::vector<GffEntry*>* pChildResults = pChild->findLevels(pPositions, pLevel, bPartialContainment);

                if (pChildResults != NULL) {
                    pLevelContained->insert(pLevelContained->end(), pChildResults->begin(), pChildResults->end());
                    delete pChildResults;
                }

            }

            return pLevelContained;

        } else {

            bool bPartiallyContained = false;
            bool bFullyContained = true;
            // correct level
            for (uint32_t i = 0; i < pPositions->size(); ++i) {

                bool bContained = this->contains(pPositions->at(i));

                bPartiallyContained |= bContained;
                bFullyContained &= bContained;

            }

            bool bAdd = (bPartialContainment == true) ? bPartiallyContained : bFullyContained;

            if (bAdd)
                pLevelContained->push_back(this);


            return pLevelContained;

        }

        return pLevelContained;

    }

    std::vector<GffEntry*>* find(std::vector<GffEntry*>* pElements, uint32_t iStart, uint32_t iEnd) {

        std::vector<GffEntry*>* pReturn = new std::vector<GffEntry*>();
        std::vector<GffEntry*>::iterator oIt;

        for (oIt = pElements->begin(); oIt != pElements->end(); ++oIt) {

            GffEntry* pElem = *oIt;

            if ((pElem->getStart() >= iStart) && (pElem->getEnd() <= iEnd)) {
                pReturn->push_back(pElem);
            }

        }

        if (pReturn->size() == 0) {
            delete pReturn;
            return NULL;
        }

        return pReturn;

    }

    GffEntry* getRegion(uint32_t iStart, uint32_t iEnd) {
        return new GffEntry(*this->getSeqName(), *this->getSource(), "region", iStart, iEnd);
    }

    int32_t compare(GffEntry* pOtherEntry) {

        int32_t iDiff = 0;

        if (pOtherEntry == NULL)
            return -1;

        if (this->m_iForwardStrand != pOtherEntry->m_iForwardStrand) {
            ++iDiff;
        }

        iDiff += std::min((uint32_t) abs(pOtherEntry->m_iStart - this->m_iStart), this->getLength());
        iDiff += std::min((uint32_t) abs(pOtherEntry->m_iEnd - this->m_iEnd), pOtherEntry->getLength());

        return iDiff;

    }

    GffEntry* addChildren(std::vector<GffEntry*>* pChildren) {

        std::cout << "Adding Children to " << *m_pFeature << " named " << *m_pSeqName << " : " << pChildren->size() << std::endl;

        std::vector<GffEntry*>::iterator oIt = pChildren->begin();

        for (; oIt < pChildren->end(); ++oIt) {
            GffEntry* pChild = (*oIt);

            pChild->setParent(this);

        }

        m_pChildren->reserve(m_pChildren->size() + pChildren->size());
        m_pChildren->insert(m_pChildren->end(), pChildren->begin(), pChildren->end());

        return this;

    }

    std::vector<GffEntry*>* getAllChildren(std::string* pFeature = NULL);

    GffEntry* addChild(GffEntry* pCandidate) {

        std::vector<GffEntry*>::iterator oIt = m_pChildren->begin();

        for (; oIt < m_pChildren->end(); ++oIt) {
            GffEntry* pChild = (*oIt);


            if (pChild->getFeature()->compare(*(pCandidate->getFeature())) == 0) {
                break;
            }

            if (pChild->contains(pCandidate)) {

                std::string sFeatureID = *(pChild->getFeature());
                sFeatureID.append("_id");

                std::string* pChildAttribute = pChild->getAttribute(sFeatureID);

                if (pChildAttribute != NULL) {
                    std::string* pCandidateAttribute = pCandidate->getAttribute(sFeatureID);

                    if (!(pCandidateAttribute == NULL)) {

                        if (pCandidateAttribute->compare(*pChildAttribute) == 0) {
                            return pChild->addChild(pCandidate);
                        }

                    } else {

                        return pChild->addChild(pCandidate);

                    }
                } else {

                    return pChild->addChild(pCandidate);

                }
            }
        }

        pCandidate->setParent(this);
        m_pChildren->push_back(pCandidate);

        return this;

    }

    static uint32_t getLengthVec(std::vector<GffEntry*>* pElements) {
        uint32_t iLength = 0;
        for (uint32_t i = 0; i < pElements->size(); ++i) {

            GffEntry* pElem = pElements->at(i);

            iLength += pElem->getLength();

        }

        return iLength;

    }

    static void deleteVec(std::vector<GffEntry*>* pElements) {

        for (uint32_t i = 0; i < pElements->size(); ++i) {

            GffEntry* pElem = pElements->at(i);

            delete pElem;

        }

        delete pElements;

    }

    static GenomicRegion* getBoundaries(std::vector<GffEntry*>* pElements) {

        uint32_t iStart = -1;
        uint32_t iEnd = 0;


        for (uint32_t i = 0; i < pElements->size(); ++i) {

            GffEntry* pElem = pElements->at(i);

            if (pElem->getStart() < iStart)
                iStart = pElem->getStart();

            if (pElem->getEnd() > iEnd)
                iEnd = pElem->getEnd();

        }

        GenomicRegion* pBounds = new GenomicRegion(iStart, iEnd);
        return pBounds;
    }

    uint32_t getStatistics(GffLoader::sStatisticElement* pElement, GffLoader::sStatisticResult* pResults);

    bool hasChild(std::string* pFeature) {
        for (std::vector<GffEntry*>::iterator oIt = m_pChildren->begin(); oIt != m_pChildren->end(); ++oIt) {
            GffEntry* pChild = *oIt;

            if (pChild->getFeature()->compare(*pFeature) == 0)
                return true;
        }

        return false;
    }



    // getCount("gene", "transcript", 0)

    // parent must be pThisFeature, how many children with feature pChildFeature exist?

    uint32_t getCount(std::string* pThisFeature, std::string* pChildFeature, uint32_t* pFoundInstances) {

        // this is not parent feature, go down further
        if (!this->m_pFeature->compare(*pThisFeature) == 0) {
            uint32_t iTotalCount = 0;

            std::vector<GffEntry*>::iterator oIt = m_pChildren->begin();

            for (; oIt < m_pChildren->end(); ++oIt) {
                GffEntry* pChild = (*oIt);

                uint32_t iReturnCount = pChild->getCount(pThisFeature, pChildFeature, pFoundInstances);

                iTotalCount += iReturnCount;
            }


            return iTotalCount;
        }

        // it's this feature

        if (pFoundInstances != NULL)
            *pFoundInstances = *pFoundInstances + 1;

        std::vector<GffEntry*>* pElements = this->getChildren(pChildFeature);
        if (pElements == NULL)
            return 0;

        return pElements->size();

    }

    std::vector<GffEntry*>* getChildren(std::string sFeature) {
        return this->getChildren(&sFeature);
    }

    std::vector<GffEntry*>* getChildren(std::string* pFeature = NULL) {

        if (pFeature == NULL)
            return m_pChildren;


        std::vector<GffEntry*>* pFoundChildren = new std::vector<GffEntry*>();
        std::vector<GffEntry*>::iterator oIt = m_pChildren->begin();


        for (; oIt < m_pChildren->end(); ++oIt) {
            GffEntry* pChild = (*oIt);

            if (pFeature->compare(*pChild->getFeature()) == 0) {
                pFoundChildren->push_back(pChild);
            }
        }

        if (pFoundChildren->size() > 0) {
            return pFoundChildren;
        }

        if (pFeature->compare("intron") == 0) {

            std::string* psExon = new std::string("exon");
            std::vector<GffEntry*>* pExons = this->getChildren(psExon);

            delete pFoundChildren;
            pFoundChildren = this->fillEmptyRegion(this, pExons, std::string("intron"));

            delete pExons;
            delete psExon;

            return pFoundChildren;
        }

        return pFoundChildren;

    }

    std::vector<GffEntry*>* fillEmptyRegion(GffEntry* pRoot, std::vector<GffEntry*>* pNegElements, std::string sFeature) {
        if (pNegElements->size() <= 1) {
            return NULL;
        }

        std::vector<GffEntry*>* pReturn = new std::vector<GffEntry*>();

        GffEntry* pOldElement = pNegElements->at(0);

        if (pRoot->getStart() != pOldElement->getStart()) {
            // add new intron
            GffEntry* pGap = new GffEntry(pRoot, true);

            pGap->setFeature(new std::string(sFeature));
            pGap->setSource(new std::string("splitter(auto)"));

            pGap->setStart(pRoot->getStart());
            pGap->setEnd(pOldElement->getStart() - 1);

            pGap->setAutogenerated(true);

            pReturn->push_back(pGap);
        }

        uint32_t iNext = 1;
        GffEntry* pElement = NULL;
        while (iNext < pNegElements->size()) {
            pElement = pNegElements->at(iNext);
            ++iNext;

            uint32_t iDist = pElement->getStart() - pOldElement->getEnd();

            if (iDist > 1) {
                GffEntry* pGap = new GffEntry(pRoot, true);

                pGap->setFeature(new std::string(sFeature));
                pGap->setSource(new std::string("splitter(auto)"));

                pGap->setStart(pOldElement->getEnd() + 1);
                pGap->setEnd(pElement->getStart() - 1);

                pGap->setAutogenerated(true);

                if (pGap->getStart() <= pGap->getEnd())
                    pReturn->push_back(pGap);
            }

            pOldElement = pElement;

        }

        if (pOldElement->getEnd() != pRoot->getEnd()) {
            GffEntry* pGap = new GffEntry(pRoot, true);

            pGap->setFeature(new std::string(sFeature));
            pGap->setSource(new std::string("splitter(auto)"));

            pGap->setStart(pOldElement->getEnd() + 1);
            pGap->setEnd(pRoot->getEnd());

            pGap->setAutogenerated(true);

            if (pGap->getStart() <= pGap->getEnd())
                pReturn->push_back(pGap);
        }

        return pReturn;
    }

    uint32_t getMaxLocation() {

        uint32_t iMax = this->getEnd();

        std::vector<GffEntry*>::iterator oIt = m_pChildren->begin();

        for (; oIt < m_pChildren->end(); ++oIt) {
            GffEntry* pChild = (*oIt);

            if (pChild->getEnd() > iMax)
                iMax = pChild->getEnd();

        }

        return iMax;

    }

    void printHierarchy(uint32_t iLevel = 0) {

        this->sortChildren();

        std::string sPrefix = "";
        sPrefix.insert(0, iLevel, '-');

        std::cout << sPrefix << *m_pFeature << "\t" << m_iStart << "\t" << m_iEnd << "\t" << this->getLength() << std::endl;

        std::vector<GffEntry*>::iterator oIt = m_pChildren->begin();

        for (; oIt < m_pChildren->end(); ++oIt) {
            GffEntry* pChild = (*oIt);

            pChild->printHierarchy(iLevel + 1);

        }


    }

    void printEntry(uint32_t iDepth = -1);
    bool sortChildren(std::vector< std::pair<std::string, std::string> >* pExpands = NULL);
    void setInTranscriptContained(bool bValue);
    bool getInTranscriptContained();

protected:


    // NULL if empty!
    std::string* m_pSeqName;
    std::string* m_pSource;
    std::string* m_pFeature;

    // -1 if empty
    uint8_t m_iFrame;

    float m_fScore;
    int8_t m_iForwardStrand;

    bool m_bAutoGenerated;
    bool m_bInTranscriptContained;

    std::map<std::string, std::string>* m_pAttributes;

    std::vector<GffEntry*>* m_pChildren;
    std::vector<GffEntry*>* m_pRegions;
    std::vector<GffTranscript*>* m_pTranscripts;

    GffEntry* m_pParent;

    void parseAttributes(std::string& sAttributes);
    void parseAttribute(std::string& sAttribute);



private:

    void initialize();
    std::vector<GffEntry*>* fillEmptyRegion(std::string sFeature);

    bool checkForEmptyValue(std::string* pString) {
        if ((pString->length() == 1) && (pString->compare(".") == 0)) {
            delete m_pSeqName;
            m_pSeqName = NULL;

            return true;
        }

        return false;
    }

};

#endif /* GFFENTRY_H_ */
