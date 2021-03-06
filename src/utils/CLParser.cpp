/* 
 * File:   CLParser.cpp
 * Author: joppich
 * 
 * Created on September 11, 2015, 9:53 AM
 */

#include <string>
#include <iostream>
#include <map>

#include "CLParser.h"

#include "Utils.h"

CLParser::CLParser(std::string sArgs)
{

    const std::string sInput = std::string(sArgs);
    std::vector<std::string> vArgs = Utils::split(sInput, ' ');
    vArgs.insert(vArgs.begin(), std::string("."));

    bool m_bSuccess = this->initialize(vArgs.size(), Utils::vec2char(vArgs));



}

CLParser::CLParser(int argc, char** argv) {

    bool m_bSuccess = this->initialize(argc, argv);

}

CLParser::CLParser(const CLParser& orig) {
}

CLParser::~CLParser() {
}

bool CLParser::initialize(int argc, char** argv)
{
    m_pThisExecutable = NULL;
    m_pCLArguments = NULL;

    if (argc > 0)
        m_pThisExecutable = new std::string(argv[0]);

    m_pCLArguments = new std::map< std::string, std::string* >();

    std::map< std::string, std::string* >::iterator oIt;
    bool bAwaitsInput = false;

    for (uint32_t i = 1; i < argc; ++i) {

        std::string* pArgument = new std::string(argv[i]);

        // does it start with - or -- ?
        if ((!bAwaitsInput) && (pArgument->at(0) == pArgument->at(1) == '-')) {

            std::string sArgument = pArgument->substr(2, pArgument->length());

            oIt = m_pCLArguments->insert(m_pCLArguments->end(), std::pair<std::string, std::string*>(sArgument, NULL));
            delete pArgument;

            bAwaitsInput = false;

            continue;
        }

        if ((!bAwaitsInput) && (pArgument->at(0) == '-')) {

            std::string sArgument = pArgument->substr(1, pArgument->length());

            oIt = m_pCLArguments->insert(m_pCLArguments->end(), std::pair<std::string, std::string*>(sArgument, NULL));
            delete pArgument;

            bAwaitsInput = true;
            continue;
        }

        if (!bAwaitsInput) {
            this->handleParsingError(pArgument);

            return false;
        }

        // bAwaitsInput == true

        (*oIt).second = pArgument;
        bAwaitsInput = false;


    }
}
