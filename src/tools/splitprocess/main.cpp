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

#include <gtxloader/GffEntry.h>
#include <gtxloader/GffLoader.h>
#include <gtxloader/GffSplitLocator.h>
#include <utils/CLParser.h>
#include <utils/Utils.h>

/*
 * 
 */
int main(int argc, char** argv) {


    CLParser* pConfig = new CLParser(argc, argv);

    GffLoader *pLoader = new GffLoader(pConfig);

    std::string* pInputFile = pConfig->getArgument("splits");
    std::string* pOutputFile = pConfig->getArgument("out");
    GffSplitLocator* pLocator = new GffSplitLocator(*pInputFile, *pOutputFile, pLoader);

    pLocator->start(NULL);


    delete pLoader;
    delete pConfig;

    return 0;

}

