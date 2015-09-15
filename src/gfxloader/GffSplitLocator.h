/*
 * Copyright 2015 <copyright holder> <email>
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#ifndef GFFSPLITLOCATOR_H
#define GFFSPLITLOCATOR_H

#include "../Utils/LineProcessor.h"
#include "../Utils/FileWriter.h"
#include "GffEntry.h"
#include "GffTranscript.h"
#include "GffLoader.h"
#include <string>
#include <vector>
#include <sstream>

class GffSplitLocator : public LineProcessor {
public:

    GffSplitLocator(std::string sFileName, std::string sOutName, GffLoader* pGffLoader)
    : LineProcessor(sFileName), m_pGffLoader(pGffLoader)
    {
        m_pLineElements = new std::vector<std::string>();
        m_pSearchLevel = new std::string("gene");
        m_pFileWriter = new FileWriter(sOutName);
        
        m_iLines = 0;
        
    }
    
    ~GffSplitLocator()
    {
        
        std::cout << "Finished " << m_iLines << " lines" << std::endl;
        
        delete m_pFileWriter;
    }

private:

    void process(std::string& sLine, void* pData) {

        if (m_iLines++ == 0)
            return;
        
        this->split(sLine, '\t', m_pLineElements);

        // chrom = 0, strand = 1, start = 2, end = 3
        
        uint32_t iStart = std::stoi( m_pLineElements->at(1) );
        uint32_t iEnd = std::stoi( m_pLineElements->at(2) );
                
        GffEntry* pChrom = m_pGffLoader->getChromosome( &(m_pLineElements->at(0)) );
        
        std::vector<uint32_t> vPos;
        vPos.push_back(iStart);
        vPos.push_back(iEnd);
        
        std::vector<GffEntry*>* pGenes = pChrom->findLevels(&vPos, m_pSearchLevel, true);
        
        uint32_t iCoveredTranscripts = 0;
        uint32_t iTranscripts = 0;
        
        std::vector<GffEntry*>::iterator oIt;
        std::vector<GffTranscript*>* pAllTranscripts = new std::vector<GffTranscript*>();
        std::vector<GffTranscript*>* pSplicedTranscripts = new std::vector<GffTranscript*>();
        
        for (oIt = pGenes->begin(); oIt != pGenes->end(); ++oIt)
        {
            
            GffEntry* pGene = *oIt;
            
            bool bCoveredTranscript = pGene->hasTranscript(&vPos, false);
            
            if (bCoveredTranscript)
                ++iCoveredTranscripts;
            
            std::vector<GffTranscript*>* pTranscripts = pGene->hasTranscript(&vPos, true);
            
            if ( pTranscripts->size() > 0)
            {
                ++iTranscripts;
                
                pAllTranscripts->insert(pAllTranscripts->end(), pTranscripts->begin(), pTranscripts->end());
                
                std::vector<GffTranscript*>::iterator oTr = pTranscripts->begin();
                bool bSplitFound = false;
                for (oTr; oTr != pTranscripts->end(); ++oTr)
                {
                    
                    GffTranscript* pTrans = *oTr;
                    
                    bSplitFound |= pTrans->hasSplit(iStart, iEnd, 1);
                    
                    if (bSplitFound)
                    {                      
                        pSplicedTranscripts->push_back(pTrans);
                    }
                    
                }    
            }       
        }
        
        std::stringstream oOutLine; // get rid of endline
        oOutLine << sLine.substr(0, sLine.length()-1) << '\t' << pGenes->size() << '\t' << iCoveredTranscripts << '\t' << iTranscripts << '\t' << pSplicedTranscripts->size() << '\t';
        
        for (uint32_t i = 0; i < pGenes->size(); ++i)
        {
            
            if (i > 0)
                oOutLine << ",";
            oOutLine << *(pGenes->at(i)->getAttribute("gene_id"));
            
        }
        
        if (pGenes->size() == 0)
            oOutLine << "";
        
        oOutLine << '\t';
        
        for (uint32_t i = 0; i < pAllTranscripts->size(); ++i)
        {
            
            if (i > 0)
                oOutLine << ",";
            
            GffTranscript* pTrans = pAllTranscripts->at(i);
            
            oOutLine << "["<< pTrans->getStart() << "," << pTrans->getEnd() << "]";
            
        }
        
        if (pAllTranscripts->size() == 0)
            oOutLine << "";
        
        oOutLine << '\t';
        
        for (uint32_t i = 0; i < pSplicedTranscripts->size(); ++i)
        {
            
            if (i > 0)
                oOutLine << ",";
            
            GffTranscript* pTrans = pSplicedTranscripts->at(i);
            
            oOutLine << *pTrans->getTranscriptID();
            
        }
        
        if (pSplicedTranscripts->size() == 0)
            oOutLine << "";
        
        oOutLine << std::endl;
        
        std::string sOutLine = oOutLine.str();
        
        m_pFileWriter->writeDirect( sOutLine );
        
        delete pGenes;
        delete pAllTranscripts;
        delete pSplicedTranscripts;
        
        m_pLineElements->clear();
    }

    std::vector<std::string>* m_pLineElements;
    
    uint32_t m_iLines;
    std::string* m_pSearchLevel;
    GffLoader* m_pGffLoader;
    FileWriter* m_pFileWriter;
    
};

#endif // GFFSPLITLOCATOR_H
