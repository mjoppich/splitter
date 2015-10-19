/*
 * Copyright 2015 <copyright holder> <email>
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#ifndef GFFSPLITLOCATOR_H
#define GFFSPLITLOCATOR_H

#include "../utils/LineProcessor.h"
#include <utils/FileWriter.h>
#include "GffEntry.h"
#include "GffTranscript.h"
#include "GffLoader.h"
#include <string>
#include <vector>
#include <sstream>

class GffSplitLocator : public LineProcessor {
public:

    GffSplitLocator(std::string sFileName, std::string sOutName, GffLoader* pGffLoader)
    : LineProcessor(sFileName), m_pGffLoader(pGffLoader)
    {
        m_pLineElements = new std::vector<std::string>();
        m_pSearchLevel = new std::string("gene");
        m_pFileWriter = new FileWriter(sOutName);
        
        m_iLines = 0;
        
    }
    
    ~GffSplitLocator()
    {
        
        std::cout << "Finished " << m_iLines << " lines" << std::endl;
        
        delete m_pFileWriter;
    }

private:

    void process(std::string& sLine, void* pData);

    std::vector<std::string>* m_pLineElements;
    
    uint32_t m_iLines;
    std::string* m_pSearchLevel;
    GffLoader* m_pGffLoader;
    FileWriter* m_pFileWriter;
    
};

#endif // GFFSPLITLOCATOR_H
