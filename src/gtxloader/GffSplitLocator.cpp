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

#include "GffSplitLocator.h"

void GffSplitLocator::process(std::string& sLine, void* pData) {

    if (m_iLines++ == 0)
        return;

    m_pLineElements->clear();

    this->split(sLine, '\t', m_pLineElements);
    // chrom = 0, strand = 1, start = 2, end = 3

    uint32_t iStart = std::stoi( m_pLineElements->at(1) );
    uint32_t iEnd = std::stoi( m_pLineElements->at(2) );

    GffEntry* pChrom = m_pGffLoader->getChromosome( &(m_pLineElements->at(0)) );


    GenomicRegion* pSplit = new GenomicRegion(iStart, iEnd);

    std::vector<GffEntry*>* pGenes = pChrom->findChildrenAt(pSplit, m_pSearchLevel, false);
    sSplitResults* pResult = NULL;
    std::stringstream oOutLine;

    if (pGenes->size() == 1)
    {

        /**
         * A single gene explains this split
         */

        pResult = this->queryGeneSame(pGenes->at(0), pSplit);

    }

    if (pGenes->size() > 1)
    {

        for (uint32_t i = 0; i < pGenes->size(); ++i)
        {

            GffEntry* pGene = pGenes->at(i);
            sSplitResults* pGeneResult = this->queryGeneSame(pGene, pSplit);

            if (pResult == NULL) {
                pResult = pGeneResult;
            } else {

                if ((pResult->bSameGene == false) && (pGeneResult->bSameGene == true))
                {
                    pResult = pGeneResult;
                    continue;
                }

                if ((pResult->bExonsOfGene == false) && (pGeneResult->bExonsOfGene == true))
                {
                    pResult = pGeneResult;
                    continue;
                }

                if ((pResult->bSameTranscriptInOneGene == false) && (pGeneResult->bSameTranscriptInOneGene == true))
                {
                    pResult = pGeneResult;
                    continue;
                }

            }


        }


    }

    if (pResult != NULL)
    {

        oOutLine << sLine.substr(0, sLine.length()-1) << pResult->toString('\t') << std::endl;
        std::string sFullLine = oOutLine.str();
        m_pFileWriter->writeDirect( sFullLine );

        return;

    }

    pGenes = pChrom->findChildrenAt(pSplit, m_pSearchLevel, true);
    std::vector<GffEntry*>* pRegions = NULL;
    GffEntry* pRegion = NULL;


    if (pGenes->size() == 1)
    {

        /**
         * A single gene explains this split
         */
        GffEntry* pGene = pGenes->at(0);

        if (pGene->contains( pSplit->getStart() ))
            pRegions = pChrom->findChildrenAt(NULL, pSplit->getEnd());
        else
            pRegions = pChrom->findChildrenAt(NULL, pSplit->getStart());

        if (pRegions->size() > 0)
            pRegion = pRegions->at(0);

        pResult = this->queryGeneDifferent(pGene, pRegion, pSplit);

    }

    if (pGenes->size() > 1)
    {

        for (uint32_t i = 0; i < pGenes->size(); ++i)
        {

            GffEntry* pGene = pGenes->at(i);

            if (pGene->contains( pSplit->getStart() ))
                pRegions = pChrom->findChildrenAt(NULL, pSplit->getEnd());
            else
                pRegions = pChrom->findChildrenAt(NULL, pSplit->getStart());

            if (pRegions->size() > 0)
                pRegion = pRegions->at(0);

            sSplitResults* pGeneResult = this->queryGeneDifferent(pGene, pRegion, pSplit);

            if (pResult == NULL) {
                pResult = pGeneResult;
            } else {

                if ((pResult->bSameGene == false) && (pGeneResult->bSameGene == true))
                {
                    pResult = pGeneResult;
                    continue;
                }

                if ((pResult->bExonsOfGene == false) && (pGeneResult->bExonsOfGene == true))
                {
                    pResult = pGeneResult;
                    continue;
                }

                if ((pResult->bSameTranscriptInOneGene == false) && (pGeneResult->bSameTranscriptInOneGene == true))
                {
                    pResult = pGeneResult;
                    continue;
                }
            }


        }


    }

    if (pResult != NULL)
    {

        oOutLine << sLine.substr(0, sLine.length()-1) << pResult->toString('\t') << std::endl;
        std::string sFullLine = oOutLine.str();
        m_pFileWriter->writeDirect( sFullLine );

        return;

    }

    m_pLineElements->clear();
}

GffSplitLocator::sSplitResults *GffSplitLocator::queryGeneSame(GffEntry* pGene, GenomicRegion* pSplit) {
    GffSplitLocator::sSplitResults* pResult = new GffSplitLocator::sSplitResults();

    pResult->pMatchingGene = pGene;
    pResult->bSameGene = true;

    if (pGene->contains( pSplit->getStart() ))
        pResult->b5pOfGene = pGene->is5pLocated( pSplit->getEnd() );
    else
        pResult->b5pOfGene = pGene->is5pLocated( pSplit->getStart() );

    std::vector<GffEntry*>* pExonsStart = pGene->findChildrenAt( "exon", pSplit->getStart());
    std::vector<GffEntry*>* pExonsEnd = pGene->findChildrenAt( "exon", pSplit->getEnd());

    if ((pExonsStart->size() > 0) && (pExonsEnd->size() > 0))
        pResult->bExonsOfGene = true;
    else
        pResult->bExonsOfGene = false;

    delete pExonsStart;
    delete pExonsEnd;

    std::vector<GffTranscript*>* pTranscripts = pGene->findTranscript( pSplit );
    pResult->bSameTranscriptInOneGene = (pTranscripts->size() > 0);
    delete pTranscripts;

    return pResult;
}

GffSplitLocator::sSplitResults *GffSplitLocator::queryGeneDifferent(GffEntry* pGene, GffEntry* pRegion, GenomicRegion* pSplit) {
    GffSplitLocator::sSplitResults* pResult = new GffSplitLocator::sSplitResults();

    pResult->pMatchingGene = pGene;
    pResult->bSameGene = false;
    pResult->pNonmatchingPart = pRegion;

    if (pGene->contains( pSplit->getStart() ))
        pResult->b5pOfGene = pGene->is5pLocated( pSplit->getEnd() );
    else
        pResult->b5pOfGene = pGene->is5pLocated( pSplit->getStart() );

    std::vector<GffEntry*>* pExonsStart = pGene->findChildrenAt( "exon", pSplit->getStart());
    std::vector<GffEntry*>* pExonsEnd = pGene->findChildrenAt( "exon", pSplit->getEnd());

    if ((pExonsStart->size() > 0) && (pExonsEnd->size() > 0))
        pResult->bExonsOfGene = true;
    else
        pResult->bExonsOfGene = false;

    delete pExonsStart;
    delete pExonsEnd;

    std::vector<GffTranscript*>* pTranscripts = pGene->findTranscript( pSplit );
    pResult->bSameTranscriptInOneGene = (pTranscripts->size() > 0);
    delete pTranscripts;

    return pResult;
}