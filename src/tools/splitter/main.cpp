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

#include <gtxloader/GffEntry.h>
#include <gtxloader/GffLoader.h>
#include <splits/XAMSplitStreamer.h>
#include <splits/FeatureManager.h>
#include <gtxloader/GffSplitLocator.h>

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


	std::string sGFFFile;
	std::string sSplitInfoFile;
	std::string sSplitInfoOutFile;

	if (argc != 4)
	{
		std::cerr << "You need to specify gene annotation and input/output splitinfo!" << std::endl;
		std::cerr << "call: splitter mygenome.gff star.splitinfo star.geneinfo" << std::endl;
		std::cerr << "You submitted " << argc << " arguments: " << std::endl;

		for (uint32_t i = 0; i < argc; ++i)
		{
			std::cerr << argv[i] << std::endl;
		}
	}

	if (argc == 4)
	{
		sGFFFile = std::string(argv[1]);
		sSplitInfoFile = std::string(argv[2]);
		sSplitInfoOutFile = std::string(argv[3]);


	} else {
		sGFFFile = std::string("/home/markus/references/Saccharomyces_cerevisiae.R64-1-1.80.gtf");
		sGFFFile = std::string("/usr/local/storage/references/Saccharomyces_cerevisiae.R64-1-1.80.gtf");

		sSplitInfoFile = "/usr/local/storage2/snyder/DNA_damage/100812_SPADE_FC62757_L3_pf/tophat.splitinfo2";
		sSplitInfoOutFile = "/usr/local/storage2/snyder/DNA_damage/100812_SPADE_FC62757_L3_pf/tophat.genesplit";
	}

	GffLoader* pLoader = NULL;//new GffLoader(sGFFFile, NULL);//new GffLoader(sGFFFile, pIgnores);
	delete pIgnores;

	pLoader->printStatistics( NULL );

	return 0;

	std::vector< std::pair< std::string, std::string > >* pExpandVector = new std::vector< std::pair< std::string, std::string > >();
	pExpandVector->push_back( std::pair< std::string, std::string >("transcript", "intron") );
	pExpandVector->push_back( std::pair< std::string, std::string >("gene", "intertrans") );
	pExpandVector->push_back( std::pair< std::string, std::string >("chromosome", "intergenic") );

	GffSplitLocator* pSplitLocator = new GffSplitLocator( sSplitInfoFile, sSplitInfoOutFile, pLoader );

	pSplitLocator->start(NULL);
	delete pSplitLocator;
	return 0;


	std::vector<std::string>* pAvailSeqNames = pLoader->getSeqNames();
	XAMSplitStreamer* pStreamer = new XAMSplitStreamer(sBAMFile, sBAMidxFile);
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

		GffEntry* pChromosome = pLoader->getChromosome( &sSeqName );

		std::cout << *(pChromosome->getFeature()) << " " << *pFlattenLevel << " " << pChromosome->getChildren()->size() << pChromosome->getLength() << std::endl;

		pChromosome->flatten( pFlattenLevel, NULL );

		pStreamer->process(iBAMHeaderIndex, pChromosome);

		delete pChromosome;

	}

	delete pFlattenLevel;

	delete pExpandVector;
	delete pLoader;
	delete pStreamer;

	return 0;
}
