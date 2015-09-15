/*
 * Utils.h
 *
 *  Created on: Jul 13, 2015
 *      Author: joppich
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <inttypes.h>
#include <vector>
#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

template <typename T1, typename T2>
class Utils {
public:

    static uint32_t find(std::vector< std::pair<T1, T2> >* pVec, T1* pElem1, T2* pElem2) {

        if (pVec == NULL)
            return -1;

        if (pElem1 != NULL) {

            for (uint32_t i = 0; i < pVec->size(); ++i) {
                if (pVec->at(i).first == *pElem1)
                    return i;
            }

        } else {

            for (uint32_t i = 0; i < pVec->size(); ++i) {
                if (pVec->at(i).second == *pElem2)
                    return i;
            }

        }

        return -1;

    }

    static std::vector< std::string >* readByLine(std::string* pFilename) {

        std::ifstream oInputStream;
        std::string sLine;
        oInputStream.open(pFilename->c_str());

        if (oInputStream.eof() == true)
            return NULL;

        std::vector< std::string >* pVec = new std::vector< std::string >();

        while (!oInputStream.eof()) {

            std::getline(oInputStream, sLine);

            if (sLine.length() == 0)
                continue;

            pVec->push_back(sLine);

        }

        return pVec;
    }

    static std::vector<std::string> &split(const std::string &sString, char cDelim, std::vector<std::string> &vElems) {
        std::stringstream sStringStream( sString );
        std::string sItem;
        while (std::getline(sStringStream, sItem, cDelim)) {
            vElems.push_back(sItem);
        }
        return vElems;
    }

    static std::vector<std::string> split(const std::string &sString, char cDelim) {
        std::vector<std::string> vElems;
        Utils<T1,T2>::split(sString, cDelim, vElems);
        return vElems;
    }

};

#endif /* UTILS_H_ */