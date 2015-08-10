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

#ifndef LINEPROCESSOR_H
#define LINEPROCESSOR_H

#include <fstream>
#include <string>


class LineProcessor
{
public:
LineProcessor( std::string sFileName);

  void start(std::string& sLine, void* pData);

private:
      virtual process(std::string& sLine, void* pData );

  
  std::string m_sFileName;
  std::ifstream m_pInFile;
};

#endif // LINEPROCESSOR_H
