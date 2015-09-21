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
#include <utils/CLParser.h>
#include <utils/Utils.h>

/*
 * 
 */
int main(int argc, char** argv) {


    CLParser* pConfig = new CLParser(argc, argv);

    GffLoader *pLoader = new GffLoader(pConfig); //new GffLoader(sGFFFile, pIgnores);

    pLoader->run();


    if (pConfig->isSet("stats") == true) {
        std::string *pStatsTSV = pConfig->getArgument("stats");

        std::cerr << "Stats File: " << *pStatsTSV << std::endl;

        std::vector<std::string> *pStats = Utils<int, int>::readByLine(pStatsTSV);
        std::vector<GffLoader::sStatisticElement *> *pStatPairs = NULL;
        if (pStats != NULL) {
            pStatPairs = new std::vector<GffLoader::sStatisticElement *>();

            for (uint32_t i = 0; i < pStats->size(); ++i) {
                GffLoader::sStatisticElement *pElement = new GffLoader::sStatisticElement();

                std::string sLine = pStats->at(i);

                std::vector<std::string> vElems = Utils<int, int>::split(sLine, '\t');

                pElement->sParent = vElems.at(0);
                pElement->sBase = vElems.at(1);
                pElement->iModifier = atoi(vElems.at(2).c_str());

                pStatPairs->push_back(pElement);

            }

        }

        pLoader->printStatistics(pStatPairs);

        if (pStats != NULL) {
            pStats->clear();

            delete pStats;
        }

    }


    if (pConfig->isSet("validate") == true)
    {
        pLoader->printValidation();
    }


    delete pLoader;
    delete pConfig;

    return 0;

}

