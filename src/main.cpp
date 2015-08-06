//============================================================================
// Name        : splitter.cpp
// Author      : Markus Joppich
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <iostream>

#include <omp.h>

#include "../GffLoader/GffEntry.h"
#include "../GffLoader/GffLoader.h"
#include "../SplitFinder/XAMSplitStreamer.h"
#include "../SplitFinder/FeatureManager.h"

void printArray(double* pArray, int l)
{
	for (int i = 0; i < l; ++i)
	{
		std::cout << pArray[i] << " ";
	}

	std::cout << std::endl;
}

double* merge(double* pLeft, double* pRight, int iLeft, int iRight)
{

	double *pRet = (double*) malloc(sizeof(double) * (iLeft + iRight));

	int k = 0;
	int l = 0;

	for (int i = 0; i < (iLeft + iRight); ++i)
	{
		if ((l == iRight) || (pLeft[k] < pRight[l]))
		{
			pRet[i] = pLeft[k++];
			continue;
		}

		if ((k == iLeft) || (pLeft[k] >= pRight[l]))
		{
			pRet[i] = pRight[l++];
			continue;
		}
	}

	free(pLeft);
	free(pRight);

	return pRet;

}

double* mergesort(double* a, int l)
{

	if (l == 1) return a;

	int iLeft = l/2;
	int iRight = l - iLeft;

	double* pLeft = (double*) malloc(sizeof(double) * iLeft);
	double* pRight = (double*) malloc(sizeof(double) * iRight);

	memcpy( pLeft, a, sizeof(double) * iLeft );
	memcpy( pRight, a + iLeft , sizeof(double) * iRight );

	printArray(pLeft, iLeft);
	printArray(pRight, iRight);

	free(a);

	return merge( mergesort(pLeft, iLeft), mergesort(pRight, iRight), iLeft, iRight );

}




int main(int argc, char** argv) {

	double* pArray = (double*) malloc(sizeof(double) * 6);

	pArray[0] = 5;
	pArray[1] = 4;
	pArray[2] = 3;
	pArray[3] = 2;
	pArray[4] = 1;
	pArray[5] = 0;

	//printArray(pArray, 6);

	//pArray = mergesort(pArray, 6);

	//printArray(pArray, 6);

	free(pArray);

	std::string sBAMFile = "accepted_hits.bam";
	std::string sBAMidxFile = "accepted_hits.bam.bai";

	std::vector<std::string>* pIgnores = new std::vector<std::string>();
	pIgnores->push_back("start_codon");
	pIgnores->push_back("stop_codon");
	pIgnores->push_back("CDS");


	std::string sGFFfile;

	if (argc > 1)
	{
		sGFFfile = std::string(argv[1]);
	} else {
		sGFFfile = "/home/markus/references/Saccharomyces_cerevisiae.R64-1-1.80.gtf";
		//sGFFfile = "/usr/local/storage/references/Homo_sapiens.GRCh38.81.gtf";
	}

	GffLoader* pLoader = new GffLoader(sGFFfile, pIgnores, "");

	delete pIgnores; // ce5.sorted.bam ce5.sorted.bam.bai

	XAMSplitStreamer* pStreamer = new XAMSplitStreamer(sBAMFile, sBAMidxFile);

	std::vector< std::pair< std::string, std::string > >* pExpandVector = new std::vector< std::pair< std::string, std::string > >();

	pExpandVector->push_back( std::pair< std::string, std::string >("transcript", "intron") );
	pExpandVector->push_back( std::pair< std::string, std::string >("gene", "intertrans") );
	pExpandVector->push_back( std::pair< std::string, std::string >("chromosome", "intergenic") );


	std::vector<std::string>* pAvailSeqNames = pLoader->getSeqNames();

	pStreamer->printSeqNames();

	std::string* pFlattenLevel = new std::string( "gene" );

	for (size_t i = 0; i < 1; ++i) // pAvailSeqNames->size()
	{

		std::string sSeqName = pAvailSeqNames->at(i);

		uint32_t iBAMHeaderIndex = pStreamer->getSeqNameID(&sSeqName);
		uint32_t iSeqLength = pStreamer->getSeqLength(&sSeqName);

		if (iBAMHeaderIndex == -1)
		{
			std::cerr << "Could not find sequences for " << sSeqName << std::endl;
			continue;
		}

		std::vector<GffEntry*>* pGffEntries = pLoader->getEntriesForSeqName(&sSeqName);

		GffEntry* pChromosome = new GffEntry( sSeqName, "splitter", "chromosome", 1, iSeqLength );
		pChromosome->addChildren( pGffEntries );
		pChromosome->sortChildren( NULL );

		pChromosome->getChildren()->size();

		std::cout << *(pChromosome->getFeature()) << " " << *pFlattenLevel << " " << pChromosome->getChildren()->size() << pChromosome->getLength() << std::endl;

		pChromosome->flatten( pFlattenLevel );

		pStreamer->process(iBAMHeaderIndex, pChromosome);

		delete pChromosome;

	}

	delete pFlattenLevel;

	delete pExpandVector;
	delete pLoader;
	delete pStreamer;

	return 0;
}
