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
#include "GffTranscript.h"

class GffEntry : public GenomicRegion {
public:

    GffEntry(std::string* pLine);
    GffEntry(GffEntry* pEntry, bool bCopyAttributes);
    GffEntry(std::string sSeqName, std::string sSource, std::string sFeature, uint32_t iStart, uint32_t iEnd);

    virtual ~GffEntry();

    static const uint8_t cSep = '\t';

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

    bool hasTranscript(std::vector<uint32_t>* pPositions) {
        if (pPositions == NULL)
            return true;

        if (m_pTranscripts == NULL)
            return false;

        bool bTranscriptFound = false;

        for (uint32_t i = 0; i < m_pTranscripts->size(); ++i) {

            GffTranscript* pTranscript = m_pTranscripts->at(i);

            bool bStartContained = pTranscript->contains(pPositions->at(0));
            bool bEndContained = pTranscript->contains(pPositions->at(pPositions->size() - 1));

            if (!bStartContained || !bEndContained)
                continue;

            bTranscriptFound = true;

        }

        return bTranscriptFound;
    }

    std::vector<GffTranscript*>* hasTranscript(std::vector<uint32_t>* pPositions, bool bHasPartialContainment) {

        std::vector<GffTranscript*>* pReturn = new std::vector<GffTranscript*>();

        bool bPartiallyContained = false;
        bool bFullyContained = false;

        for (uint32_t i = 0; i < m_pTranscripts->size(); ++i) {
            GffTranscript* pTranscript = m_pTranscripts->at(i);

            for (uint32_t i = 0; i < pPositions->size(); ++i) {

                bool bContained = pTranscript->contains(pPositions->at(i));

                bPartiallyContained |= bContained;
                bFullyContained &= bContained;

            }
            bool bAdd = (bHasPartialContainment == true) ? bPartiallyContained : bFullyContained;

            if (bAdd)
                pReturn->push_back(pTranscript);

        }

        return pReturn;

    }

    std::vector<GffEntry*>* findLevels(std::vector<uint32_t>* pPositions, std::string* pLevel, bool bPartialContainment = false) {

        std::vector<GffEntry*>* pLevelContained = new std::vector<GffEntry*>();

        if (this->m_pFeature->compare(*pLevel) != 0) {

            for (uint32_t i = 0; i < m_pChildren->size(); ++i) {

                GffEntry* pChild = m_pChildren->at(i);

                std::vector<GffEntry*>* pChildResults = pChild->findLevels(pPositions, pLevel);

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

        delete pLevelContained;

        return NULL;

    }

    std::vector<GffEntry*>* find(std::vector<GffEntry*>* pElements, uint32_t iStart, uint32_t iEnd) {

        std::vector<GffEntry*>* pReturn = new std::vector<GffEntry*>();
        std::vector<GffEntry*>::iterator oIt;

        for (oIt = pElements->begin(); oIt != pElements->end(); ++oIt) {

            GffEntry* pElem = *oIt;

            if ((iStart >= pElem->getStart()) && (iEnd <= pElem->getEnd())) {
                pReturn->push_back(pElem);
            }

        }

        if (pReturn->size() == 0) {
            delete pReturn;
            return NULL;
        }

        return pReturn;

    }

    std::vector<GffEntry*>* split(uint32_t iPosition) {
        if ((iPosition <= m_iStart) || (iPosition >= m_iEnd))
            return 0;

        std::vector<GffEntry*>* pNewRegions = new std::vector<GffEntry*>();

        GffEntry* pRegion1 = new GffEntry(*this->getSeqName(), *this->getSource(), "region", m_iStart, iPosition);
        GffEntry* pRegion2 = new GffEntry(*this->getSeqName(), *this->getSource(), "region", iPosition + 1, m_iEnd);

        pNewRegions->push_back(pRegion1);
        pNewRegions->push_back(pRegion2);

        return pNewRegions;
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

    GffEntry* addChild(GffEntry* pCandidate) {

        std::vector<GffEntry*>::iterator oIt = m_pChildren->begin();

        for (; oIt < m_pChildren->end(); ++oIt) {
            GffEntry* pChild = (*oIt);

            if (pChild->contains(pCandidate)) {
                return pChild->addChild(pCandidate);
            }

        }

        pCandidate->setParent(this);
        m_pChildren->push_back(pCandidate);

        return this;

    }

    std::vector<GffEntry*>* getChildren() {
        return m_pChildren;
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

    struct {

        bool operator()(GffEntry* pElem1, GffEntry* pElem2) {
            if (pElem1->getStart() < pElem2->getStart())
                return true;

            if (pElem1->getStart() > pElem2->getStart())
                return false;

            return pElem1->getLength() > pElem2->getLength();
        }
    } sSortEntriesAsc;

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
