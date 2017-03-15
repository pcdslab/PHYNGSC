/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
*/

#include "defs.h"
#include "bit_stream.h"
#include <algorithm>
#include <stdio.h>


// --------------------------------------------------------------------------------------------
BitStream::BitStream()
	:	IO_BUFFER_SIZE(0)
	,	file_pos(0)
	,	io_buffer(NULL)
	,	io_buffer_pos(0)
	,	io_buffer_size(0)
	,	word_buffer(0)
	,	word_buffer_pos(0)
	,	word_buffer_size(32)
{
	for (int32 i = 0; i < 32; ++i)
	{
		n_bit_mask[i] = (1u << i) - 1;
	}
}

// --------------------------------------------------------------------------------------------
BitStream::~BitStream()
{
	if (io_buffer)
		delete[] io_buffer;
}

// --------------------------------------------------------------------------------------------
bool BitStream::WriteBuffer()
{
	io_buffer_pos = 0;
	return true;
}

// --------------------------------------------------------------------------------------------
bool BitStream::Create(uint32 buffer_type)
{
	IO_BUFFER_SIZE = buffer_type == 0 ? 1 << 15 : 1 << 20 ;
	io_buffer = new uchar[IO_BUFFER_SIZE];
	io_buffer_size = IO_BUFFER_SIZE;

	word_buffer_size = 32;
	io_buffer_pos = 0;
	file_pos = 0;

	return true;	
}

// --------------------------------------------------------------------------------------------
bool BitStream::Close()
{
	if (io_buffer)
	{
		delete[] io_buffer;
		io_buffer = NULL;
	}
	io_buffer_pos = 0;
	return true;
}

