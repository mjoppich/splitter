/*
 * ThreadedHashMap.h
 *
 *  Created on: Jul 20, 2015
 *      Author: joppich
 */

#ifndef THREADEDHASHMAP_H_
#define THREADEDHASHMAP_H_

#include <inttypes.h>
#include <iostream>
#include <malloc.h>
#include <vector>
#include <omp.h>

template <class K, class O>
class ThreadedHashMap
{
public:
	ThreadedHashMap(uint32_t iBins)
	{

		m_iBins = iBins;
		m_pBins = (std::vector<O> **) malloc(sizeof(std::vector<O>**) * iBins);
		m_pLocks = (omp_lock_t*) malloc(sizeof(omp_lock_t) * iBins);
		m_pLoadFactors = (float*) malloc(sizeof(float) * iBins);

		for (uint32_t i = 0; i < m_iBins; ++i)
		{
			m_pBins[i] = new std::vector<O>();
			omp_init_lock (&(m_pLocks[i]));
		}

	}

	virtual float getAvgLoadFactor()
	{

		float fTotalSize = 0;

		for (uint32_t i = 0; i < m_iBins; ++i)
		{
			m_pLoadFactors[i] = m_pBins[i]->size();
			fTotalSize += m_pLoadFactors[i];
		}

		std::cerr << "Total Elements stored: " << fTotalSize << std::endl;

		float fAvgLoadFactor = fTotalSize;

		/*

		for (u32 i = 0; i < m_iBins; ++i)
		{
			m_pLoadFactors[i] = m_pLoadFactors[i];
			fAvgLoadFactor += m_pLoadFactors[i];
		}
		*/

		fAvgLoadFactor = fAvgLoadFactor / (float) m_iBins;

		return fAvgLoadFactor;
	}

    virtual uint64_t getElementCount()
    {
        uint64_t iElements = 0;

        for (uint32_t i = 0; i < m_iBins; ++i)
        {
            iElements += m_pBins[i]->size();
        }

        return iElements;
    }

	virtual ~ThreadedHashMap()
	{

		if (m_pBins != NULL)
		{
			for (uint32_t i = 0; i < m_iBins; ++i)
			{

				for (size_t j = 0; j < m_pBins[i]->size(); ++j)
				{
					delete m_pBins[i]->at(j);
				}

				delete m_pBins[i];
				omp_destroy_lock(&(m_pLocks[i]));

			}

			free(m_pLocks);
			free(m_pBins);
			free(m_pLoadFactors);
		}
	}

	O insert(K oKey, O oObject)
	{


		uint32_t iIndex = this->getIndexForKey(oKey);

		omp_set_lock (&(m_pLocks[iIndex]));

		m_pBins[iIndex]->push_back( oObject );
		O oReturnObject = m_pBins[iIndex]->back();

		omp_unset_lock(&(m_pLocks[iIndex]));

		return oReturnObject;

	}

	std::vector<O>* retrieveObjectsUnguarded(K oKey)
	{

		uint32_t iIndex = this->getIndexForKey(oKey);

		std::vector<O> *pReturn = findElements( m_pBins[iIndex], oKey);

		return pReturn;

	}
	std::vector<O>* retrieveObjects(K oKey)
	{

		uint32_t iIndex = this->getIndexForKey(oKey);

		omp_set_lock (&(m_pLocks[iIndex]));
		std::vector<O> *pReturn = findElements( m_pBins[iIndex], oKey);
		omp_unset_lock(&(m_pLocks[iIndex]));


		return pReturn;

	}

protected:
	virtual std::vector<O> *findElements(std::vector<O>* pLookIn, K oKey) = 0;
	virtual uint32_t getIndexForKey(K okey) = 0;

	std::vector<O>** m_pBins;
	omp_lock_t* m_pLocks;

	uint32_t m_iBins;
	float* m_pLoadFactors;

};

#endif /* THREADEDHASHMAP_H_ */
