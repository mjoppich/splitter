/* 
 * File:   FileWriter.cpp
 * Author: joppich
 * 
 * Created on August 11, 2015, 11:12 AM
 */

#include "FileWriter.h"

FileWriter::FileWriter(std::string sFileName) {
    m_bSelfCreated = true;

    m_pOut = new std::ofstream();

    m_pOut->open(sFileName);

    m_bSelfCreated = true;

    pthread_mutex_init(&mLock, NULL);

    m_iBufferSize = 65536;
    m_pBuffer = (char*) malloc(sizeof (char) * m_iBufferSize);
    m_pOut->rdbuf()->pubsetbuf(m_pBuffer, m_iBufferSize);
}

FileWriter::~FileWriter() {
    pthread_mutex_destroy(&mLock);
    m_pOut->close();
    delete m_pOut;
    free(m_pBuffer);

}

