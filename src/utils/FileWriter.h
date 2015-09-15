/* 
 * File:   FileWriter.h
 * Author: joppich
 *
 * Created on August 11, 2015, 11:12 AM
 */

#ifndef FILEWRITER_H
#define	FILEWRITER_H

#include <pthread.h>
#include <ostream>
#include <fstream>
#include <inttypes.h>

class FileWriter {
public:
    FileWriter(std::string sFileName);

    virtual ~FileWriter();

    void writeDirect(std::string& sString) {
        pthread_mutex_lock(&mLock);

        (*m_pOut) << sString;
        //m_pOut->flush();

        pthread_mutex_unlock(&mLock);

        return;

    }

private:

    std::ofstream* m_pOut;
    bool m_bSelfCreated;
    pthread_mutex_t mLock;

    uint32_t m_iBufferSize;
    char* m_pBuffer;

};

#endif	/* FILEWRITER_H */

