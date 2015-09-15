/* 
 * File:   mainGTFTool.cpp
 * Author: joppich
 *
 * Created on September 11, 2015, 12:15 PM
 */

#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <iostream>

#include <omp.h>

#include "../GffLoader/GffEntry.h"
#include "../GffLoader/GffLoader.h"
#include "../Utils/CLParser.h"
#include "../Utils/Utils.h"

/*
 * 
 */
int main(int argc, char** argv) {


    CLParser* pConfig = new CLParser(argc, argv);
    std::map< std::string, std::string* >* pArguments = pConfig->getArguments();

    if (pArguments == NULL) {
        std::cerr << "Problem with configuration:" << std::endl;

        for (uint32_t i = 0; i < argc; ++i) {
            std::cerr << argv[i] << std::endl;
        }

        return 1;
    }

    std::string* pGTF = pConfig->getArgument("gtf");

    if (pGTF == NULL) {
        std::cerr << "Problem with configuration: missing argument gtf" << std::endl;

        return 2;
    }
    
    std::cerr << "GTF File: " << *pGTF << std::endl;
    
    std::string* pStatsTSV = pConfig->getArgument("stats");

    if (pStatsTSV == NULL) {
        std::cerr << "Problem with configuration: missing argument stats" << std::endl;

        return 2;
    }
    
    std::cerr << "Stats File: " << *pStatsTSV << std::endl;

    std::vector< std::string >* pStats = Utils<int,int>::readByLine( pStatsTSV );
    std::vector< GffLoader::sStatisticElement* >* pStatPairs = NULL;
    if (pStats != NULL)
    {
        pStatPairs = new std::vector< GffLoader::sStatisticElement* >();
        
        for (uint32_t i = 0; i < pStats->size(); ++i)
        {
            GffLoader::sStatisticElement* pElement = new GffLoader::sStatisticElement();
            
            std::string sLine = pStats->at(i);
                        
            std::vector<std::string> vElems = Utils<int,int>::split(sLine, '\t');
            
            pElement->sParent = vElems.at(0);
            pElement->sBase = vElems.at(1);
            pElement->iModifier = atoi(vElems.at(2).c_str());
            
            pStatPairs->push_back( pElement );
            
        }
        
    }
    
    
    
    
    
    
    GffLoader* pLoader = new GffLoader(*pGTF, NULL); //new GffLoader(sGFFFile, pIgnores);    
    pLoader->printStatistics( pStatPairs );
    
    if (pStats != NULL)
    {
        pStats->clear();
        
        delete pStats;
    }

    delete pLoader;
    delete pConfig;

    return 0;

}

