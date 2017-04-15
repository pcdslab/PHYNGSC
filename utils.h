/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
*/

#ifndef _UTILS_H
#define _UTILS_H

#include "defs.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <string.h>

namespace utils
{

// --------------------------------------------------------------------------------------------
inline uint32 to_string(uchar* str, uint32 value)
{
	uint32 digits;
	uint32 power = 1;

	if (value == 0)
	{
		str[0] = '0';
		return 1;
	}

	for (digits = 0; digits < 10; ++digits)
	{
		if (value < power)
			break;
		power *= 10;
	}

	power /= 10;
	for (uint32 i = 0; power; ++i, power /= 10)
	{
		int32 d = value / power;
		str[i] = (uchar)('0' + d);
		value -= d * power;
	}

	return digits;
}

// --------------------------------------------------------------------------------------------
inline bool extend_string(uchar *&str, uint32 &size)
{
	uint32 new_size = size * 2;
	uchar *p = new uchar[new_size+1];

	if (!p)
		return false;

	std::copy(str, str+size, p);
	size = new_size;
	delete[] str;
	str = p;

	return true;
}

// --------------------------------------------------------------------------------------------
inline bool extend_string_to(uchar *&str, uint32 &size, uint32 new_size)
{
	if (new_size <= size)
		return true;

	uchar *p = new uchar[new_size+1];

	if (!p)
		return false;

	std::copy(str, str+size, p);

	size = new_size;
	delete[] str;
	str = p;

	return true;
}

// --------------------------------------------------------------------------------------------
inline uint32 int_log(uint32 x, uint32 base)
{
	uint32 r = 0;

	if (base == 0)
		return 1;
	if (base == 1)
		base++;

	for (uint32 tmp = base; tmp <= x; tmp *= base)
		++r;

	return r;
}

// --------------------------------------------------------------------------------------------
inline bool is_num(const uchar* str, uint32 len)
{
	for (uint32 i = 0; i < len; ++i)
		if (str[i] < '0' || str[i] > '9')
			return false;

	return (len == 1) || (len > 1 && str[0] != '0');
}

// --------------------------------------------------------------------------------------------
inline uint32 to_num(const uchar *str, uint32 len)
{
	uint32 r = 0;

	for (uint32 i = 0; i < len; ++i)
		r = r * 10 + (str[i] - '0');

	return r;
}

}

// --------------------------------------------------------------------------------------------
//
// --------------------------------------------------------------------------------------------
#if defined(_WIN32)
#if !defined(copy_n)	// TODO: check condition on different WIN32 compilers
#include <algorithm>
#include <iterator>

template <class _InputIter, class _Size, class _OutputIter>
void copy_n(_InputIter first, _Size count, _OutputIter result) 
{
	std::copy(first, first + count, result);
}
#endif

#else

#include <ext/algorithm>
using __gnu_cxx::copy_n;

#endif

#endif