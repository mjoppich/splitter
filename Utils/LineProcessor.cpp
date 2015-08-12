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

#include <string>
#include <fstream>
#include "../Utils/LineProcessor.h"

LineProcessor::LineProcessor(std::string sFileName)
{

  m_sFileName = sFileName;
  
  m_pInFile = new std::ifstream();
  m_pInFile->open( sFileName.c_str() );
  
}

void LineProcessor::start(void* pData)
{
  
  std::string sLine;
  
  // TODO for parallelisation read in chunks and then start thread for each chunk
  while (std::getline( *m_pInFile, sLine) )
  {
    
    this->process( sLine, pData);
  }
  
}
