/*
 * GffTranscript.h
 *
 *  Created on: Jul 31, 2015
 *      Author: joppich
 */

#ifndef GFFTRANSCRIPT_H_
#define GFFTRANSCRIPT_H_

#include <string>
#include <vector>
#include <inttypes.h>

#include "GenomicRegion.h"

class GffEntry;

class GffTranscript : public GenomicRegion {
public:

    GffTranscript(std::string sTransID, uint32_t iStart, uint32_t iEnd)
    : GenomicRegion(iStart, iEnd) {
        m_pTranscriptID = new std::string(sTransID);

        m_pExons = new std::vector<GffEntry*>();
    }

    void addExons(std::vector<GffEntry*>* pExons) {
        if (pExons == NULL)
            return;

        m_pExons->insert(m_pExons->end(), pExons->begin(), pExons->end());
    }

    void addExon(GffEntry* pExon) {
        if (pExon == NULL)
            return;

        m_pExons->push_back(pExon);
    }

    std::vector<GffEntry*>* getExons() {
        return m_pExons;
    }
    
    std::string* getTranscriptID()
    {
        return m_pTranscriptID;
    }
    
    bool hasSplit(uint32_t iStart, uint32_t iEnd, uint32_t iConfInterval = 0);

private:

    std::vector<GffEntry*>* m_pExons;

    std::string* m_pTranscriptID;

};

#endif /* GFFTRANSCRIPT_H_ */
