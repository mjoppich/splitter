/*
 * XAMSplitStreamer.h
 *
 *  Created on: Jul 7, 2015
 *      Author: joppich
 */

#ifndef XAMSPLITSTREAMER_H_
#define XAMSPLITSTREAMER_H_

#include <string>
#include <inttypes.h>
#include <vector>
#include <map>
#include <htslib/sam.h>
#include "SplitEvent.h"

class GffEntry;
class SplitEvent;
class SplitHashMap;
class FeatureManager;

class XAMSplitStreamer {
public:
	XAMSplitStreamer(std::string sBAMFile, std::string sBAMidxFile);

	uint32_t getSeqNameID(std::string* pString);
	uint32_t getSeqLength(std::string* pString);
	FeatureManager* process(uint32_t iSeqID, GffEntry* pChromosomeEntry, std::vector<std::string>* pFeatureHierarchy = NULL);

	virtual ~XAMSplitStreamer();

	void printSeqNames();

private:

	void processInputFiles(std::string* pBAMFile, std::string* pBAMidxFile);

	SplitEvent* evaluateAlignment(bam_hdr_t* pHeader, uint32_t iSeqID, bam1_t* pAlignment, GffEntry* pChromosome);

	std::string* m_pBAMFile;
	std::string* m_pBAMidxFile;

	std::vector<std::string>* m_pFeatureHierarchy;
	SplitHashMap* m_pSplits;

	FeatureManager* m_pFeatures;


};

#endif /* XAMSPLITSTREAMER_H_ */
