/*
 * GffTranscript.cpp
 *
 *  Created on: Jul 31, 2015
 *      Author: joppich
 */

#include "GffTranscript.h"
#include "GffEntry.h"

bool GffTranscript::hasSplit(uint32_t iStart, uint32_t iEnd, uint32_t iConfInterval)
    {
        
        for (uint32_t i = 0; i < m_pExons->size()-1; ++i)
        {
            GffEntry* pExon = m_pExons->at(i);
            GffEntry* pNextExon = m_pExons->at(i+1);
            
            uint32_t iStartDiff = std::abs( (int64_t)pExon->getEnd() - (int64_t)iStart);
            uint32_t iEndDiff = std::abs( (int64_t)pNextExon->getStart() - (int64_t)iEnd);
            
            bool bStart = !( iStartDiff > iConfInterval );
            bool bEnd = !( iEndDiff > iConfInterval );
            
            if ( bStart && bEnd)
                return true;
        }
        
        return false;
        
    }