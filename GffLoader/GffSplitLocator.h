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

#ifndef GFFLOCATOR_H
#define GFFLOCATOR_H

#include "../Utils/LineProcessor.h"
#include <string>
#include <vector>

class GffLocator : public LineProcessor
{
public:
GffLocator(std::string sFileName)
 : LineProcessor(sFileName)
 {
   m_pLineElements = new std::vector<std::string>();
 }

private:
  
void process(std::string& sLine, void* pData)
{
  this->split(sLine, '\t', m_pLineElements);
  
  
  
  m_pLineElements->clear();
}

  std::vector<std::string>* m_pLineElements;
  
};

#endif // GFFLOCATOR_H
