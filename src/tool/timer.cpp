/************************************************************************
Copyright Paris13, CNRS

This software is a computer program whose purpose is to [describe
functionalities and technical features of your software].

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
************************************************************************/

/**
* \file    src/core/timer.cpp
* \brief   Timer
* \author  N. Lermé
* \version 1.0
* \date    09-2008
*
*/

#include "timer.h"


CTimer::CTimer()
{
	Init();
}

void CTimer::Init()
{
	#ifdef _WIN32
		QueryPerformanceFrequency(&m_TicksPerSecond);
		QueryPerformanceCounter(&m_AbsoluteTime);
	#else
		gettimeofday(&m_AbsoluteTime,NULL);
	#endif
}

double CTimer::Time() const
{
	double ElapsedTime;

	#ifdef _WIN32
		LARGE_INTEGER NewTime;
		QueryPerformanceCounter(&NewTime);
		ElapsedTime = (double(NewTime.QuadPart) - double(m_AbsoluteTime.QuadPart)) / double(m_TicksPerSecond.QuadPart);
	#else
		timeval NewTime;
		gettimeofday(&NewTime, NULL);
		ElapsedTime = (NewTime.tv_usec - m_AbsoluteTime.tv_usec) / 1000000. + NewTime.tv_sec - m_AbsoluteTime.tv_sec;
	#endif

	return ElapsedTime;
}
