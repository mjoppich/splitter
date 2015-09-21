#ifndef CLRUNNABLE_H
#define CLRUNNABLE_H

#include "CLParser.h"
#include <inttypes.h>

class CLRunnable
{

public:
    CLRunnable(CLParser* pParser)
    : m_pParser(pParser)
    {

    }

    virtual void run() = 0;

protected:

    CLParser* getConfig();

    CLParser* m_pParser;

    virtual bool checkConfig() = 0;
    virtual uint32_t prepareRun(CLParser* pParser) = 0;


};

#endif // CLRUNNABLE_H
