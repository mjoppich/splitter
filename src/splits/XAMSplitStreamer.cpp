/*
 * XAMSplitStreamer.cpp
 *
 *  Created on: Jul 7, 2015
 *      Author: joppich
 */

#include "XAMSplitStreamer.h"


#include <string>
#include <iostream>
#include <set>
#include <vector>

#include "../GffLoader/GffEntry.h"
#include "SplitHashMap.h"

#include <htslib/sam.h>
#include "SplitEvent.h"

#include "FeatureManager.h"



XAMSplitStreamer::XAMSplitStreamer(std::string sBAMFile, std::string sBAMidxFile)
: m_pBAMFile( new std::string(sBAMFile)), m_pBAMidxFile(new std::string(sBAMidxFile)), m_pFeatures(NULL)
{

	m_pFeatureHierarchy = new std::vector<std::string>();

	m_pFeatureHierarchy->push_back("intron");
	m_pFeatureHierarchy->push_back("exon");
	m_pFeatureHierarchy->push_back("transcript");
	m_pFeatureHierarchy->push_back("gene");

}

void XAMSplitStreamer::processInputFiles(std::string* pBAMFile, std::string* pBAMidxFile)
{
	samFile *pInFile = sam_open(pBAMFile->c_str(), "r");
	bam_hdr_t* pHeader = sam_hdr_read(pInFile);
	kstring_t ks = { 0, 0, NULL };

	for (size_t i = 0; i < pHeader->n_targets; ++i)
	{
		std::cout << i << " " << pHeader->target_name[i] << " " << pHeader->target_len[i] << std::endl;
	}


	hts_idx_t* pIdx = bam_index_load(pBAMidxFile->c_str());

	uint32_t iTID = 0; // extracts all alignments with name header->target_name[iTID]

	hts_itr_t* potherit = bam_itr_queryi( pIdx ,iTID ,0 ,10000);

	bam1_t* pRes = bam_init1();

	while (bam_itr_next(pInFile, potherit, pRes) > 0)
	{

		uint8_t iBase = bam_seqi(bam_get_seq(pRes), 1);
		std::cout << "bla: " << pRes->core.tid << " " << (int) iBase << " " << pRes->core.pos << std::endl;

		if (sam_format1(pHeader, pRes, &ks) >= 0)
		{
			std::cout << ks.s << std::endl;
		}

	}

	hts_itr_t* pIt = sam_itr_queryi(NULL, HTS_IDX_REST, 0, 0);

	hts_itr_destroy(pIt);
	hts_itr_destroy(sam_itr_queryi(NULL, HTS_IDX_NONE, 0, 0));
}

SplitEvent* XAMSplitStreamer::evaluateAlignment(bam_hdr_t* pHeader, uint32_t iSeqID, bam1_t* pAlignment, GffEntry* pChromosome)
{
	uint32_t iLeftStart = pAlignment->core.pos;
	uint32_t iCIGARitems = pAlignment->core.n_cigar;

	uint32_t* pCIGARs = bam_get_cigar(pAlignment);

	int32_t iSeqLength = 0;

	std::vector<GffEntry*> pMatchingGffEntriesPStrand;
	std::vector<GffEntry*> pMatchingGffEntriesNStrand;


	kstring_t ks = { 0, 0, NULL };

	if (sam_format1(pHeader, pAlignment, &ks) >= 0)
	{
		std::cout << ks.s << std::endl;
	}

	std::vector<uint32_t> vPositions, vOperations;
	vPositions.push_back( iLeftStart );
	vOperations.push_back( 0 );


	// TODO this might be too much sampling! not all introns might be recognized
	for (size_t i = 0; i < iCIGARitems; ++i)
	{
		uint32_t iCIGARop = bam_cigar_op(pCIGARs[i]);
		char cCIGARop = BAM_CIGAR_STR[ iCIGARop ];
		uint32_t iCIGARlen = bam_cigar_oplen(pCIGARs[i]);

		if ((iCIGARop == 1) || (iCIGARop == 4))
			continue;

		if (iCIGARop == 0) // MATCH M
			iSeqLength = iSeqLength + iCIGARlen;

		if (iCIGARop == 1) // INSERTION TO REFERENCE I
			iSeqLength = iSeqLength + 0;

		if (iCIGARop == 2) // DELETION FROM REFERENCE D
			iSeqLength = iSeqLength + iCIGARlen;

		if (iCIGARop == 3) // SKIPPED REGION IN REFERENCE N
			iSeqLength = iSeqLength + iCIGARlen;

		if (iCIGARop == 4) // SOFT CLIPPING S
			iSeqLength = iSeqLength + 0;

		if (iCIGARop == 5) // HARD CLIPPING H
			iSeqLength = iSeqLength + iCIGARlen;

		if (iCIGARop == 6) // PADDING P
			iSeqLength = iSeqLength + iCIGARlen;

		if (iCIGARop == 7) // SEQUENCE MATCH =
			iSeqLength = iSeqLength + iCIGARlen;

		if (iCIGARop == 8) // SEQUENCE MISMATCH X
			iSeqLength = iSeqLength + iCIGARlen;

		std::cout << "ITEM = " << i << " " << iCIGARop << "=" << cCIGARop << " " << iCIGARlen << " " << iLeftStart+iSeqLength << std::endl;

		vPositions.push_back( iLeftStart + iSeqLength );
		vOperations.push_back( iCIGARop );

	}


	std::string sStopDepth = "gene";

	// is this read fully contained in some gene?
	std::vector<GffEntry*>* pTouchedGenes = pChromosome->findLevels(&vPositions, &sStopDepth);

	// if no -> intragenic
	if ((pTouchedGenes == NULL) || (pTouchedGenes->size() == 0))
	{
		// TODO something if intergenic

		pTouchedGenes = pChromosome->findLevels(&vPositions, &sStopDepth, true);

		SplitEvent* pReturn = NULL;

		if ((pTouchedGenes != NULL) && (pTouchedGenes->size() > 0))
		{

			uint32_t iStart = 0;
			uint32_t iCIGAROpPos = 0;
			uint32_t iMaxContained = 0;
			uint32_t iMinOutside = 0;
			bool bFirstInside = false;

			for (uint32_t i = 0; i < pTouchedGenes->size(); ++i)
			{

				GffEntry* pGene = pTouchedGenes->at(i);
				bool bFirstIn = pGene->contains( vPositions.at(0) );

				for (uint32_t j = 1; j < vPositions.size(); ++j)
				{

					bool bContained = pGene->contains( vPositions.at(j) );

					if (bContained == bFirstIn)
						continue;

					uint32_t iCurContained = (bFirstIn) ? vPositions.at(j) : iSeqLength - vPositions.at(j);
					uint32_t iCurOutside = iSeqLength - iCurContained;

					if (iCurContained > iMaxContained)
					{
						iMaxContained = iCurContained;
						iMinOutside = iCurOutside;
						bFirstInside = bFirstIn;
						iStart = iLeftStart + vPositions.at(j);
						iCIGAROpPos = j;
					}
				}

			}

			if (iMaxContained > 0)
			{
				uint32_t iEnd = iStart;
				if (bFirstInside)
				{

					if (vOperations.at(iCIGAROpPos+1) == 3 ) // jump
						iEnd = vPositions.at(iCIGAROpPos+1);

				} else {

					if (vOperations.at(iCIGAROpPos-1) == 3 ) // jump
						iEnd = vPositions.at(iCIGAROpPos-1);

				}

				if (bFirstInside)
				{
					pReturn = m_pFeatures->addFeatureEntry("intergenic_split", pAlignment, iStart, iEnd, iMaxContained, iMinOutside);
				} else {
					pReturn = m_pFeatures->addFeatureEntry("intergenic_split", pAlignment, iStart, iEnd, iMinOutside, iMaxContained);
				}
			}

		} else {

			pReturn = m_pFeatures->addFeatureEntry("intergenic_read", pAlignment, iLeftStart, iLeftStart+iSeqLength, 0,0);

		}


		return pReturn;
	}

	// if yes
	// does there exist any transcript such that all positions on transcript?
	bool bCoveredTranscript = false;
        bool bReadFullyCovered = false;
        
        std::vector< GffTranscript* > * pFullTranscripts = new std::vector<GffTranscript*>();
        std::vector< GffTranscript* > * pTranscripts = new std::vector<GffTranscript*>();
                
	for (uint32_t i = 0; i < pTouchedGenes->size(); ++i)
	{
		GffEntry* pGene = pTouchedGenes->at(i);
                
                std::vector<GffTranscript*>* pReturnTranscripts = pGene->hasTranscript( &vPositions, false );
                
                if (pReturnTranscripts->size() > 0)
                {
                    bReadFullyCovered |= true;
                    pFullTranscripts->insert(pFullTranscripts->end(), pReturnTranscripts->begin(), pReturnTranscripts->end());
                    break;
                }
                
                delete pReturnTranscripts;
                
                // not fully covered!
                pReturnTranscripts = pGene->hasTranscript( &vPositions, true );
                
                if (pReturnTranscripts->size() > 0)
                    pTranscripts->insert(pTranscripts->end(), pReturnTranscripts->begin(), pReturnTranscripts->end());
                
                delete pReturnTranscripts;
	}
        
	// if yes -> all ok
	// one might have a split, but everything is explained by a single transcript
	if (bReadFullyCovered)
        {
            delete pTranscripts;
            delete pFullTranscripts;
            
            std::cout << std::endl << std::endl << "Valid Transcript" << std::endl << std::endl;
            
		return NULL;
        } 

	//if no either intron retention or split

	// 1) exists position such that position is exon and next position not an exon in any case (or vice versa)? -> IR
	// 1a) if N on CIGAR position before intron: SP
	// 2) exists position such that anything before in gene, but then intergenic -> IRI
	// 2a) if N on CIGAR before leaving gene: SPI

            std::cout << std::endl << std::endl << "Non Valid Transcript" << std::endl << std::endl;

}

FeatureManager* XAMSplitStreamer::process(uint32_t iSeqID, GffEntry* pChromosomeEntry, std::vector<std::string>* pFeatureHierarchy)
{
	if (pFeatureHierarchy == NULL)
		pFeatureHierarchy = m_pFeatureHierarchy;

	samFile *pInFile = sam_open(m_pBAMFile->c_str(), "r");
	bam_hdr_t* pHeader = sam_hdr_read(pInFile);
	kstring_t ks = { 0, 0, NULL };

	std::cout << iSeqID << " " << pHeader->target_name[iSeqID] << " " << pHeader->target_len[iSeqID] << std::endl;

	hts_idx_t* pIdx = bam_index_load(m_pBAMidxFile->c_str());
	hts_itr_t* potherit = bam_itr_queryi( pIdx ,iSeqID ,0 ,10000);
	bam1_t* pRes = bam_init1();

	if (m_pFeatures != NULL)
		delete m_pFeatures;

	m_pFeatures = new FeatureManager();
        uint32_t iReads = 0;

	while (bam_itr_next(pInFile, potherit, pRes) > 0)
	{

		this->evaluateAlignment(pHeader, iSeqID, pRes, pChromosomeEntry);
                ++iReads;
	}
        
        
        std::cout << "Reads processed: " << iReads << std::endl;
        m_pFeatures->printFeatures( std::cout );

	hts_itr_t* pIt = sam_itr_queryi(NULL, HTS_IDX_REST, 0, 0);

	hts_itr_destroy(pIt);
	hts_itr_destroy(sam_itr_queryi(NULL, HTS_IDX_NONE, 0, 0));

	sam_close(pInFile);
	bam_hdr_destroy(pHeader);

	return m_pFeatures;
}

void XAMSplitStreamer::printSeqNames()
{

	samFile *pInFile = sam_open(m_pBAMFile->c_str(), "r");
	bam_hdr_t* pHeader = sam_hdr_read(pInFile);

	uint32_t iResult = 0;

	for (iResult = 0; iResult < pHeader->n_targets; ++iResult)
	{
		std::cout << pHeader->target_name[iResult] << std::endl;
	}


	sam_close(pInFile);
	bam_hdr_destroy(pHeader);

}

uint32_t XAMSplitStreamer::getSeqNameID(std::string* pString)
{
	samFile *pInFile = sam_open(m_pBAMFile->c_str(), "r");
	bam_hdr_t* pHeader = sam_hdr_read(pInFile);

	uint32_t iResult = 0;

	for (iResult = 0; iResult < pHeader->n_targets; ++iResult)
	{
		if (pString->compare(pHeader->target_name[iResult]) == 0)
		{
			break;
		}
	}

	if (iResult == pHeader->n_targets)
		iResult = -1;

	sam_close(pInFile);
	bam_hdr_destroy(pHeader);

	return iResult;
}

uint32_t XAMSplitStreamer::getSeqLength(std::string* pString)
{
	samFile *pInFile = sam_open(m_pBAMFile->c_str(), "r");
	bam_hdr_t* pHeader = sam_hdr_read(pInFile);

	uint32_t iResult = 0;

	for (iResult = 0; iResult < pHeader->n_targets; ++iResult)
	{
		if (pString->compare(pHeader->target_name[iResult]) == 0)
		{
			return pHeader->target_len[iResult];
		}
	}

	sam_close(pInFile);
	bam_hdr_destroy(pHeader);

	return -1;
}

XAMSplitStreamer::~XAMSplitStreamer() {

	delete m_pFeatureHierarchy;
	delete m_pBAMFile;
	delete m_pBAMidxFile;

}

