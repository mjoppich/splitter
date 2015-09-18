/* 
 * File:   CLParser.h
 * Author: joppich
 *
 * Created on September 11, 2015, 9:53 AM
 */

#ifndef CLPARSER_H
#define	CLPARSER_H

#include <inttypes.h>

class CLParser {
public:
    CLParser(int argc, char** argv);
    CLParser(const CLParser& orig);
    virtual ~CLParser();

    std::map< std::string, std::string* >* getArguments() {
        return m_pCLArguments;
    }

    std::string* getArgument(std::string sArg) {

        std::map<std::string, std::string*>::iterator oIt = m_pCLArguments->find(sArg);

        if (oIt == m_pCLArguments->end())
            return NULL;

        return (*oIt).second;

    }

    bool isSet(std::string sArg) {
        std::map<std::string, std::string*>::iterator oIt = m_pCLArguments->find(sArg);

        if (oIt == m_pCLArguments->end())
            return false;
        
        return true;
    }

private:

    void handleParsingError(std::string* pArgument) {
        std::cerr << "Error: invalid input arguments at " << std::endl;

        this->deleteArguments();

        delete pArgument;

        m_pCLArguments = NULL;
    }

    void deleteArguments() {

        std::map< std::string, std::string*>::iterator oIt = m_pCLArguments->begin();

        for (; oIt != m_pCLArguments->end(); ++oIt) {
            delete (*oIt).second;
        }

        m_pCLArguments->empty();

        delete m_pCLArguments;

    }



    std::map< std::string, std::string* >* m_pCLArguments;
    std::string* m_pThisExecutable;

};

#endif	/* CLPARSER_H */

