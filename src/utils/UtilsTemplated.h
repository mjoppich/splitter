//
// Created by joppich on 9/21/15.
//

#ifndef PROJECT_UTILSTEMPLATED_H
#define PROJECT_UTILSTEMPLATED_H

#include <string.h>
#include <vector>
#include <inttypes.h>

template<typename T1, typename T2>
class UtilsTemplated {
public:

    static uint32_t find(std::vector<std::pair<T1, T2> > *pVec, T1 *pElem1, T2 *pElem2) {

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

};

#endif //PROJECT_UTILSTEMPLATED_H
