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

#include "GffTranscript.h"

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

GffEntry::GffEntry(std::string* pLine)
: GenomicRegion(0, 0), m_pChildren(NULL), m_pAttributes(NULL), m_pRegions(NULL), m_pTranscripts(NULL) {

    this->initialize();

    if (pLine->find('\t') != std::string::npos) {
        this->parseLine(pLine, '\t');
        return;
    }

    if (pLine->find(' ') != std::string::npos) {
        this->parseLine(pLine, ' ');
        return;
    }

}

GffEntry::GffEntry(std::string* pLine, uint8_t cSep)
: GenomicRegion(0, 0), m_pChildren(NULL), m_pAttributes(NULL), m_pRegions(NULL), m_pTranscripts(NULL) {

    this->initialize();

    this->parseLine(pLine, cSep);

}

void GffEntry::parseLine(std::string* pLine, uint8_t cSep) {

    uint32_t iOldFound = 0;
    uint32_t iFound = pLine->find(cSep, 0);
    m_pSeqName = new std::string(*pLine, iOldFound, iFound - iOldFound);
    this->checkForEmptyValue(m_pSeqName);

    iOldFound = iFound + 1;
    iFound = pLine->find(cSep, iOldFound);
    m_pSource = new std::string(*pLine, iOldFound, iFound - iOldFound);
    this->checkForEmptyValue(m_pSource);

    iOldFound = iFound + 1;
    iFound = pLine->find(cSep, iOldFound);
    m_pFeature = new std::string(*pLine, iOldFound, iFound - iOldFound);
    this->checkForEmptyValue(m_pFeature);


    iOldFound = iFound + 1;
    iFound = pLine->find(cSep, iOldFound);
    m_iStart = std::strtoul(pLine->substr(iOldFound, iFound - iOldFound).c_str(), NULL, 0); // TODO this might not work

    iOldFound = iFound + 1;
    iFound = pLine->find(cSep, iOldFound);
    m_iEnd = std::strtoul(pLine->substr(iOldFound, iFound - iOldFound).c_str(), NULL, 0); // TODO this might not work

    iOldFound = iFound + 1;
    iFound = pLine->find(cSep, iOldFound);

    std::string sScore = pLine->substr(iOldFound, iFound - iOldFound);
    if ((sScore.length() == 1) && (sScore[0] == '.'))
        m_fScore = std::numeric_limits<float>::quiet_NaN();
    else
        m_fScore = std::strtof(sScore.c_str(), NULL);

    iOldFound = iFound + 1;
    iFound = pLine->find(cSep, iOldFound);

    if (iFound - iOldFound != 1)
        std::cerr << "Error Parsing GTF v2.0 Line: Strand Information: " << *pLine << std::endl;

    char cStrandInfo = pLine->substr(iOldFound, iFound - iOldFound)[0];

    if (cStrandInfo == '.')
        m_iForwardStrand = -1;
    else
        m_iForwardStrand = (cStrandInfo == '+') ? 1 : 0;

    iOldFound = iFound + 1;
    iFound = pLine->find(cSep, iOldFound);

    if (iFound - iOldFound != 1)
        std::cerr << "Error Parsing GTF v2.0 Line: Codon: " << *pLine << std::endl;

    char cFrameInfo = pLine->substr(iOldFound, iFound - iOldFound)[0];
    if (cFrameInfo == '.')
        m_iFrame = -1;
    else
        m_iFrame = cFrameInfo - 48;

    iOldFound = iFound + 1;
    iFound = pLine->find('\n', iOldFound);

    std::string sAttributes = pLine->substr(iOldFound, iFound - iOldFound);

    this->parseAttributes(sAttributes);

}

GffEntry::GffEntry(GffEntry* pEntry, bool bCopyAttributes)
: GenomicRegion(pEntry->m_iStart, pEntry->m_iEnd), m_pChildren(NULL), m_pAttributes(NULL), m_pRegions(NULL), m_pTranscripts(NULL) {

    this->initialize();

    this->m_pSeqName = new std::string(*(pEntry->m_pSeqName));
    this->m_pFeature = new std::string(*(pEntry->m_pFeature));
    this->m_pSource = new std::string(*(pEntry->m_pSource));

    this->m_iStart = pEntry->m_iStart;
    this->m_iEnd = pEntry->m_iEnd;
    this->m_iForwardStrand = pEntry->m_iForwardStrand;
    this->m_iFrame = pEntry->m_iFrame;

    if (bCopyAttributes) {

        std::map<std::string, std::string>::iterator oIt;

        for (oIt = pEntry->m_pAttributes->begin(); oIt != pEntry->m_pAttributes->end(); ++oIt) {
            this->m_pAttributes->insert(std::pair<std::string, std::string>(oIt->first, oIt->second));
        }

    }

}

GffEntry::GffEntry(std::string sSeqName, std::string sSource, std::string sFeature, uint32_t iStart, uint32_t iEnd)
: GenomicRegion(iStart, iEnd), m_pChildren(NULL), m_pAttributes(NULL), m_pRegions(NULL), m_pTranscripts(NULL) {

    this->initialize();

    m_pSeqName = new std::string(sSeqName);
    m_pSource = new std::string(sSource);
    m_pFeature = new std::string(sFeature);

    m_iFrame = 0;
    m_fScore = 0;
    m_iForwardStrand = 1;

}

void GffEntry::initialize() {

    m_pAttributes = new std::map<std::string, std::string>();

    m_pChildren = new std::vector<GffEntry*>();

    m_bAutoGenerated = false;
    m_bInTranscriptContained = false;

    m_pRegions = NULL;
    m_pTranscripts = NULL;

}

void GffEntry::setFeature(std::string* pNewFeature) {
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

void GffEntry::printEntry(uint32_t iDepth) {

    std::string* pGeneID = this->getAttribute(std::string("gene_id"));

    std::cout << *m_pSeqName << "(" << ((pGeneID != NULL) ? *pGeneID : "NULL") << ")" << " " << *m_pSource << " " << *m_pFeature << "[" << m_iStart << " " << m_iEnd << "]" << std::endl;

    if (iDepth == -1)
        return;

    for (uint32_t i = 0; i < m_pChildren->size(); ++i)
        m_pChildren->at(i)->printEntry(iDepth - 1);



}

bool GffEntry::sortChildren(std::vector< std::pair<std::string, std::string> >* pExpands) {

    uint32_t iExpandIndex = Utils<std::string, std::string>::find(pExpands, this->m_pFeature, NULL);

    bool bPrint = false;

    if (m_pChildren->size() > 0) {
        for (size_t i = 0; i < m_pChildren->size(); ++i) {
            m_pChildren->at(i)->sortChildren(pExpands);
        }

        std::vector<GffEntry*>* pRegions = NULL;
        GffEntry::sSortEntriesAsc oSorter;

        if (iExpandIndex != -1) {


            if (m_pChildren->size() > 1) {
                std::sort(m_pChildren->begin(), m_pChildren->end(), oSorter);
            }

            pRegions = this->fillEmptyRegion(pExpands->at(iExpandIndex).second);
        }

        if (pRegions != NULL) {

            if (pRegions->size() > 1) {
                this->printEntry();
                bPrint = true;
            }

            m_pChildren->insert(m_pChildren->end(), pRegions->begin(), pRegions->end());
        }

        if (m_pChildren->size() > 1) {
            std::sort(m_pChildren->begin(), m_pChildren->end(), oSorter);
        }

        if (bPrint)
            this->printEntry();
    }


    return true;
}

std::vector<GffEntry*>* GffEntry::fillEmptyRegion(std::string sFeature) {
    if (m_pChildren->size() <= 1) {
        return NULL;
    }

    std::vector<GffEntry*>* pIntrons = new std::vector<GffEntry*>();

    GffEntry* pTranscript = this;

    GffEntry* pOldExon = m_pChildren->at(0);

    if (pTranscript->getStart() != pOldExon->getStart()) {
        // add new intron
        GffEntry* pIntron = new GffEntry(pTranscript, true);

        pIntron->setFeature(new std::string(sFeature));

        pIntron->setSource(new std::string("splitter(auto)"));

        pIntron->setStart(pTranscript->getStart());
        pIntron->setEnd(pOldExon->getStart() - 1);

        pIntron->setAutogenerated(true);

        pIntrons->push_back(pIntron);
    }

    uint32_t iNext = 1;
    GffEntry* pExon = NULL;
    while (iNext < m_pChildren->size()) {
        pExon = m_pChildren->at(iNext);
        ++iNext;

        uint32_t iDist = pExon->getStart() - pOldExon->getEnd();

        if (iDist > 1) {
            GffEntry* pIntron = new GffEntry(pTranscript, true);

            pIntron->setFeature(new std::string(sFeature));
            pIntron->setSource(new std::string("splitter(auto)"));

            pIntron->setStart(pOldExon->getEnd() + 1);
            pIntron->setEnd(pExon->getStart() - 1);

            pIntron->setAutogenerated(true);

            if (pIntron->getStart() <= pIntron->getEnd())
                pIntrons->push_back(pIntron);
        }

        pOldExon = pExon;

    }

    if (pOldExon->getEnd() != pTranscript->getEnd()) {
        GffEntry* pIntron = new GffEntry(pTranscript, true);

        pIntron->setFeature(new std::string(sFeature));
        pIntron->setSource(new std::string("splitter(auto)"));

        pIntron->setStart(pOldExon->getEnd() + 1);
        pIntron->setEnd(pTranscript->getEnd());

        pIntron->setAutogenerated(true);

        if (pIntron->getStart() <= pIntron->getEnd())
            pIntrons->push_back(pIntron);
    }

    return pIntrons;
}

void GffEntry::parseAttributes(std::string& sAttributes) {

    size_t iFound = 0;
    size_t iOldFound = 0;

    std::replace(sAttributes.begin(), sAttributes.end(), '\t', ' ');

    bool bVisitedLoop = false;

    while ((iFound = sAttributes.find(';', iOldFound)) != std::string::npos) {

        bVisitedLoop = true;

        std::string sAttribute = sAttributes.substr(iOldFound, iFound - iOldFound);
        sAttribute = sAttribute.substr(sAttribute.find_first_not_of(' '), -1);

        this->parseAttribute(sAttribute);


        iOldFound = iFound + 1;

    }

    if (!bVisitedLoop) {
        if (sAttributes.length() > 0)
            this->parseAttribute(sAttributes);
    }

}

void GffEntry::parseAttribute(std::string& sAttribute) {

    std::string sKey, sValue;
    size_t iSplit = 0;

    iSplit = sAttribute.find("=");

    if (iSplit != std::string::npos) {

        sKey = sAttribute.substr(0, iSplit);
        sValue = sAttribute.substr(iSplit + 1, -1);

    } else {
        iSplit = sAttribute.find(" ");

        if (iSplit != std::string::npos) {
            sKey = sAttribute.substr(0, iSplit);
            sValue = sAttribute.substr(iSplit + 1, -1);

            if (sValue[0] == '\"') {
                sValue = sValue.substr(1, sValue.length() - 2);
            }

        }
    }

    m_pAttributes->insert(std::pair<std::string, std::string>(sKey, sValue));

}

void GffEntry::setInTranscriptContained(bool bValue) {
    m_bInTranscriptContained = bValue;
}

bool GffEntry::getInTranscriptContained() {
    return m_bInTranscriptContained;
}

std::vector<GffEntry*>* GffEntry::getAllChildren(std::string* pFeature) {

    if (pFeature == NULL)
        return m_pChildren;


    std::vector<GffEntry*>* pFoundChildren = new std::vector<GffEntry*>();
    std::vector<GffEntry*>::iterator oIt = m_pChildren->begin();


    for (; oIt < m_pChildren->end(); ++oIt) {
        GffEntry* pChild = (*oIt);

        if (pFeature->compare(*pChild->getFeature()) == 0) {
            pFoundChildren->push_back(pChild);
        } else {

            std::vector<GffEntry*>* pFound = pChild->getAllChildren(pFeature);
            pFoundChildren->insert(pFoundChildren->end(), pFound->begin(), pFound->end());
            delete pFound;

        }
    }

    //pFoundChildren = this->fillEmptyRegion(this, pExons, std::string("intron"));

    return pFoundChildren;

}

uint32_t GffEntry::getStatistics(GffLoader::sStatisticElement* pElement, GffLoader::sStatisticResult* pResults) {

    if (!this->m_pFeature->compare(pElement->sParent) == 0) {

        std::vector<GffEntry*>::iterator oIt = m_pChildren->begin();

        for (; oIt < m_pChildren->end(); ++oIt) {
            GffEntry* pChild = (*oIt);

            pChild->getStatistics(pElement, pResults);
        }

        return 0;
    }

    std::vector<GffEntry*>* pChildren = this->getAllChildren(&(pElement->sBase));
    std::vector<GffEntry*>* pRevChildren;
    GenomicRegion* pBounds;
    GffEntry* pBoundary;
    
    switch (pElement->iModifier) {

        case 0:
            // normal case
            ++(pResults->iParentCount);
            pResults->iBaseCount += pChildren->size();
            pResults->iLengthBase += GffEntry::getLengthVec(pChildren);

            break;
        case 1:

            // reverse case

            ++(pResults->iParentCount);

            pRevChildren = this->fillEmptyRegion(this, pChildren, "intron");

            pResults->iBaseCount += pRevChildren->size();
            pResults->iLengthBase += GffEntry::getLengthVec(pRevChildren);

            GffEntry::deleteVec(pRevChildren);

            break;

        case 2:

            // reverse case, tight

            ++(pResults->iParentCount);

            pBounds = GffEntry::getBoundaries(pChildren);
            pBoundary = new GffEntry(this, true);
            pBoundary->setStart(pBounds->getStart());
            pBoundary->setEnd(pBounds->getEnd());
            pRevChildren = this->fillEmptyRegion(pBoundary, pChildren, "intron");

            pResults->iBaseCount += pRevChildren->size();
            pResults->iLengthBase += GffEntry::getLengthVec(pRevChildren);

            GffEntry::deleteVec(pRevChildren);
            delete pBounds;
            delete pBoundary;



            break;

        default: break;

    }


    return 0;

}

void GffEntry::flatten(std::string* pFlattenLevel) {

    if (m_pFeature->compare(pFlattenLevel->c_str()) != 0) {
        for (uint32_t i = 0; i < m_pChildren->size(); ++i) {
            GffEntry* pElem = m_pChildren->at(i);
            pElem->flatten(pFlattenLevel);
        }

        return;
    }

    std::vector<GffEntry*>* pRegions = new std::vector<GffEntry*>();


    std::vector<uint32_t> vExonBounds;

    vExonBounds.push_back(this->getStart());
    vExonBounds.push_back(this->getEnd() + 1);
    for (uint32_t i = 0; i < m_pChildren->size(); ++i) {

        GffEntry* pTranscript = m_pChildren->at(i);
        std::vector<GffEntry*>* pExons = pTranscript->getChildren();

        GffEntry* pLastExon;

        for (uint32_t i = 0; i < pExons->size(); ++i) {

            GffEntry* pExon = pExons->at(i);
            pLastExon = pExon;

            vExonBounds.push_back(pExon->getStart());
            vExonBounds.push_back(pExon->getEnd() + 1);

        }

    }

    std::sort(vExonBounds.begin(), vExonBounds.end());
    // deletes duplicates
    vExonBounds.erase(std::unique(vExonBounds.begin(), vExonBounds.end()), vExonBounds.end());


    //    std::cout << "THIS: " << this->getStart() << " " << this->getEnd() << " " << this->getLength() << std::endl;
    uint32_t i;
    for (i = 0; i < vExonBounds.size() - 1; ++i) {
        if (this->getStart() > vExonBounds.at(i))
            continue;

        if (this->getEnd() <= vExonBounds.at(i + 1) - 1)
            break;

        GffEntry* pRegion = new GffEntry(*(this->m_pSeqName), *(this->m_pSource), "region", vExonBounds.at(i), vExonBounds.at(i + 1) - 1);
        pRegions->push_back(pRegion);
        //std::cout << pRegion->toString() << std::endl;

    }

    GffEntry* pRegion = new GffEntry(*(this->m_pSeqName), *(this->m_pSource), "region", vExonBounds.at(i), vExonBounds.at(i + 1) - 1);
    pRegions->push_back(pRegion);
    //std::cout << pRegion->toString() << std::endl;

    std::vector<GffTranscript*>* pTranscripts = new std::vector<GffTranscript*>();

    std::vector<GffEntry*>::iterator oIt;
    for (oIt = m_pChildren->begin(); oIt != m_pChildren->end(); ++oIt) {

        GffEntry* pEntryTranscript = *oIt;
        std::string* pTID = pEntryTranscript->getAttribute("transcript_id");

        if (pTID == NULL)
        {
            
            pTID = pEntryTranscript->getAttribute("ID");
            
            
            if (pTID == NULL)
                continue;
        }
            

        GffTranscript* pTranscript = new GffTranscript(*pTID, pEntryTranscript->getStart(), pEntryTranscript->getEnd());

        std::vector<GffEntry*>::iterator oJt;
        for (oJt = pEntryTranscript->getChildren()->begin(); oJt != pEntryTranscript->getChildren()->end(); ++oJt) {


            GffEntry* pExon = *oJt;

            std::vector<GffEntry*>* pExons = this->find(pRegions, pExon->getStart(), pExon->getEnd());

            for (uint32_t i = 0; i < pExons->size(); ++i) {
                GffEntry* pElem = pExons->at(i);
                pElem->setInTranscriptContained(true);
            }

            pTranscript->addExons(pExons);

        }

        pTranscripts->push_back(pTranscript);


    }

    if ((pRegions->size() == 0) || (pTranscripts->size() == 0)) {
        std::cout << "problem " << pRegions->size() << " " << pTranscripts->size() << std::endl;
    }

    m_pRegions = pRegions;
    m_pTranscripts = pTranscripts;

}

std::vector<GffTranscript*>* GffEntry::hasTranscript(std::vector<uint32_t>* pPositions, bool bHasPartialContainment) {

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
