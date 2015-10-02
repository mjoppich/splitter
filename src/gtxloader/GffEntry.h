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

    static std::vector< std::string* >* g_pFeatures;

    struct sSortEntriesAsc {

        bool operator()(GffEntry* pElem1, GffEntry* pElem2) {
            if (pElem1->getStart() < pElem2->getStart())
                return true;

            if (pElem1->getStart() > pElem2->getStart())
                return false;

            return pElem1->getLength() > pElem2->getLength();
        }
    };

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

    std::string* getFeatureString( std::string* pNewFeature)
    {
        for (uint32_t i = 0; i < g_pFeatures->size(); ++i)
            if (g_pFeatures->at(i)->compare(*pNewFeature) == 0)
                return g_pFeatures->at(i);
        std::string* pFeature = NULL;

#pragma omp critical
        {

            // has it been inserted in the mean time?
            for (uint32_t i = 0; i < g_pFeatures->size(); ++i)
                if (g_pFeatures->at(i)->compare(*pNewFeature) == 0)
                {
                    pFeature = g_pFeatures->at(i);
                    break;
                }


            if (pFeature == NULL)
            {
                pFeature = new std::string(*pNewFeature);
                g_pFeatures->push_back( pFeature );
            }

        }

        return pFeature;
    }

    void setSource(std::string* pNewSource) {
        m_pSource = pNewSource;
    }

    void setAttribute(std::string &sKey, std::string sValue);

    std::string *getAttribute(std::string sKey);
    void printAttribute(std::string* pAttrib = NULL);
    std::string *getID();

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

    std::vector<GffEntry *> *findLevels(std::vector<uint32_t> *pPositions, std::string *pLevel,
                                        bool bPartialContainment = false);

    std::vector<GffEntry *> *find(std::vector<GffEntry *> *pElements, uint32_t iStart, uint32_t iEnd);

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

    GffEntry *addChildren(std::vector<GffEntry *> *pChildren);

    std::vector<GffEntry*>* getAllChildren(std::string* pFeature = NULL);

    GffEntry *findChildWithAttribute(std::string *pAttrib, std::string *pCandidateAttribute);

    GffEntry *addChildSimple(GffEntry *pCandidate);

    GffEntry *addChild(GffEntry *pCandidate, bool bValidate = false);

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

    static uint32_t getLengthVec(std::vector<GffEntry *> *pElements) {
        uint32_t iLength = 0;
        for (uint32_t i = 0; i < pElements->size(); ++i) {

            GffEntry *pElem = pElements->at(i);

            iLength += pElem->getLength();

        }

        return iLength;

    }

    static void deleteVec(std::vector<GffEntry *> *pElements) {

        for (uint32_t i = 0; i < pElements->size(); ++i) {

            GffEntry *pElem = pElements->at(i);

            delete pElem;

        }

        delete pElements;

    }

    uint32_t getStatistics(GffLoader::sStatisticElement* pElement, GffLoader::sStatisticResult* pResults);

    bool hasChild(std::string *pFeature);



    // getCount("gene", "transcript", 0)

    // parent must be pThisFeature, how many children with feature pChildFeature exist?

    uint32_t getCount(std::string *pThisFeature, std::string *pChildFeature, uint32_t *pFoundInstances);

    std::vector<GffEntry*>* getChildren(std::string sFeature) {
        return this->getChildren(&sFeature);
    }

    std::vector<GffEntry *> *getChildren(std::string *pFeature = NULL);

    std::vector<GffEntry *> *fillEmptyRegion(GffEntry *pRoot, std::vector<GffEntry *> *pNegElements,
                                             std::string sFeature);

    uint32_t getMaxLocation();

    void printHierarchy(uint32_t iLevel = 0, uint32_t iMaxLevel = -1);

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
