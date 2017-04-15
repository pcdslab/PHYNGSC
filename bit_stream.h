/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
*/

#ifndef _BIT_STREAM_H
#define _BIT_STREAM_H

#include "defs.h"
#include "huffman.h"
#include <vector>

// --------------------------------------------------------------------------------------------
//
// --------------------------------------------------------------------------------------------
class BitStream 
{
public:
	static const uint32 DEFAULT_IO_BUFFER_SIZE_W	=  1 << 23;

	BitStream();
	~BitStream();

	inline bool FlushFullWordBuffer();
	inline bool FlushPartialWordBuffer();
	inline bool FlushInputWordBuffer();

	inline bool PutBit(const bool word);
	inline bool PutBit(const uint32 word);
	inline bool Put2Bits(const uint32 word);
	inline bool PutBits(uint32 word, int32 n_bits);
	inline bool PutBytes(const uchar *data, int32 n_bytes);
	inline bool PutByte(const uchar byte);
	inline bool PutWord(const uint32 data);
	inline bool PutDWord(const uint64 data);
	inline bool Put2Bytes(const uint32 data);

	inline bool GetBit(uint32 &word);
	inline bool Get2Bits(uint32 &word);
	inline bool GetBits(uint32 &word, uint32 n_bits);
	inline bool GetBytes(uchar *data, uint32 n_bytes);
	inline bool GetByte(uint32 &byte);
	inline bool GetWord(uint32 &data);
	inline bool GetDWord(uint64 &data);
	inline bool Get2Bytes(uint32 &data);

	static inline uint32 BitLength(const uint64 x);

	// buffer_type: 0 -> for the processor buffer
	//              1 -> for the [title|DNA|Quality] buffer
	bool Create(uint32 buffer_type);
	bool Close();

	std::vector<uchar> GetIO_Buffer() const { return io_buffer; };
	int32 GetIO_Buffer_Pos() { return io_buffer_pos; };

	void SetIO_Buffer(uchar *read_buffer, uint32 r_buff_size, uint32 start_pos = 0)
	{ 
		io_buffer.assign(read_buffer+start_pos, read_buffer+r_buff_size); 
	};
	bool SetIO_Buffer_Pos(uint64 pos);

private:
	int32 IO_BUFFER_SIZE;

	uint64 file_pos;
	std::vector<uchar> io_buffer;
	int32 io_buffer_pos;
	uint32 word_buffer;
	int32 word_buffer_pos;
	int32 word_buffer_size;
	uint32 n_bit_mask[32];
};

// --------------------------------------------------------------------------------------------
bool BitStream::PutBit(const bool word)
{
	if (word_buffer_pos < word_buffer_size)
	{
		word_buffer <<= 1;
		if (word)
			word_buffer += 1;
		++word_buffer_pos;
	}
	else
	{
		PutWord(word_buffer);
		word_buffer_pos = 1;
		if (word)
			word_buffer = 1;
		else
			word_buffer = 0;
	}

	return true;			
}

// --------------------------------------------------------------------------------------------
bool BitStream::PutBit(const uint32 word)
{
	if (word_buffer_pos < word_buffer_size)
	{
		word_buffer <<= 1;
		word_buffer += word;
		++word_buffer_pos;
	}
	else
	{
		PutWord(word_buffer);
		word_buffer_pos = 1;
		word_buffer = word;
	}

	return true;			
}

// --------------------------------------------------------------------------------------------
bool BitStream::Put2Bits(const uint32 word)
{
	if (word_buffer_pos + 2 <= word_buffer_size)
	{
		word_buffer <<= 2;
		word_buffer += word;
		word_buffer_pos += 2;
	}
	else if (word_buffer_pos == word_buffer_size)
	{
		PutWord(word_buffer);
		word_buffer_pos = 2;
		word_buffer = word;
	}
	else
	{
		word_buffer <<= 1;
		word_buffer += word >> 1;
		PutWord(word_buffer);
		word_buffer = word & 1;
		word_buffer_pos = 1;
	}

	return true;			
}

// --------------------------------------------------------------------------------------------
bool BitStream::PutBits(uint32 word, int32 n_bits)
{
	my_assert(n_bits != 0);
	int32 rest_bits = word_buffer_size - word_buffer_pos;
	if (n_bits >= rest_bits)
	{
		n_bits -= rest_bits;
		word_buffer <<= rest_bits;
		word_buffer += word >> n_bits;
		word &= n_bit_mask[n_bits];
		word_buffer_pos = 0;
		PutWord(word_buffer);
		word_buffer = 0;
	}
	
	word_buffer <<= n_bits;
	word_buffer += word;
	word_buffer_pos += n_bits;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitStream::FlushFullWordBuffer()
{
	PutWord(word_buffer);
	
	word_buffer = 0;
	word_buffer_pos = 0;
	
	return true;
}

// --------------------------------------------------------------------------------------------
bool BitStream::FlushPartialWordBuffer()
{
	word_buffer <<= (32 - word_buffer_pos) & 7;

	if (word_buffer_pos > 24)
		PutByte(word_buffer >> 24);
	if (word_buffer_pos > 16)
		PutByte((word_buffer >> 16) & 0xFF);
	if (word_buffer_pos > 8)
		PutByte((word_buffer >> 8) & 0xFF);
	if (word_buffer_pos > 0)
		PutByte(word_buffer & 0xFF);
	
	word_buffer = 0;
	word_buffer_pos = 0;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitStream::FlushInputWordBuffer()
{
	word_buffer_pos = 0;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitStream::PutByte(const uchar byte)
{
	io_buffer.push_back(byte);
	io_buffer_pos++;
	file_pos++;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitStream::Put2Bytes(const uint32 data)
{
	PutByte((uchar) (data >> 8));
	PutByte((uchar) (data & 0xFF));

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitStream::PutBytes(const uchar *data, int32 n_bytes)
{
	for (int32 i = 0; i < n_bytes; ++i)
		io_buffer.push_back(data[i]);

	io_buffer_pos += n_bytes;
	file_pos += n_bytes;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitStream::PutDWord(const uint64 data)
{
	PutByte(data >> 56);
	PutByte((data >> 48) & 0xFF);
	PutByte((data >> 40) & 0xFF);
	PutByte((data >> 32) & 0xFF);
	PutByte((data >> 24) & 0xFF);
	PutByte((data >> 16) & 0xFF);
	PutByte((data >> 8) & 0xFF);
	PutByte(data & 0xFF);

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitStream::PutWord(const uint32 data)
{
	PutByte(data >> 24);
	PutByte((data >> 16) & 0xFF);
	PutByte((data >> 8) & 0xFF);
	PutByte(data & 0xFF);

	return true;
}

// --------------------------------------------------------------------------------------------
uint32 BitStream::BitLength(const uint64 x)
{
	for (uint32 i = 0; i < 32; ++i)
	{
		if (x < (1ull << i))
			return i;
	}

	return 64;
}

// --------------------------------------------------------------------------------------------
bool BitStream::GetBit(uint32 &word)
{
	if (word_buffer_pos == 0)
	{
		if (!GetByte(word_buffer))
			return false;
		word_buffer_pos = 7;
		word = word_buffer >> 7;
	}
	else
	{
		word = (word_buffer >> (--word_buffer_pos)) & 1;
	}

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitStream::Get2Bits(uint32 &word)
{
	if (word_buffer_pos >= 2)
	{
		word = (word_buffer >> (word_buffer_pos-2)) & 3;
		word_buffer_pos -= 2;
	}
	else if (word_buffer_pos == 0)
	{
		if (!GetByte(word_buffer))
			return false;
		word = word_buffer >> 6;
		word_buffer_pos = 6;
	}
	else
	{
		word = (word_buffer & 1) << 1;
		if (!GetByte(word_buffer))
			return false;
		word += word_buffer >> 7;
		word_buffer_pos = 7;
	}

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitStream::GetBits(uint32 &word, uint32 n_bits)
{
	my_assert(n_bits != 0);
	word = 0;
	while (n_bits)
	{
		if (word_buffer_pos == 0)
		{
			if (!GetByte(word_buffer))
				return false;
			word_buffer_pos = 8;
		}

		if ((int32) n_bits > word_buffer_pos)
		{
			word <<= word_buffer_pos;
			word += word_buffer & n_bit_mask[word_buffer_pos];
			n_bits -= word_buffer_pos;
			word_buffer_pos = 0;
		}
		else
		{
			word <<= n_bits;
			word_buffer_pos -= n_bits;
			word += (word_buffer >> word_buffer_pos) & n_bit_mask[n_bits];
			return true;
		}
	}

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitStream::GetByte(uint32 &byte)
{
	byte = io_buffer[io_buffer_pos++];
	file_pos++;
	return true;
}

// --------------------------------------------------------------------------------------------
bool BitStream::GetBytes(uchar *data, uint32 n_bytes)
{
	for (int32 i = 0; i < (int32)n_bytes; ++i)
		data[i] = io_buffer.at(io_buffer_pos+i);

	io_buffer_pos += n_bytes;
	file_pos      += n_bytes;
	
	return true;
}

// --------------------------------------------------------------------------------------------
bool BitStream::GetWord(uint32 &data)
{
	uint32 c;
	bool r;

	r = GetByte(c);
	data = c;
	r &= GetByte(c);
	data = (data << 8) + c;
	r &= GetByte(c);
	data = (data << 8) + c;
	r &= GetByte(c);
	data = (data << 8) + c;

	return r;
}

// --------------------------------------------------------------------------------------------
bool BitStream::GetDWord(uint64 &data)
{
	uint32 c;
	bool r;

	r = GetByte(c);
	data = c;
	r &= GetByte(c);
	data = (data << 8) + c;
	r &= GetByte(c);
	data = (data << 8) + c;
	r &= GetByte(c);
	data = (data << 8) + c;
	r &= GetByte(c);
	data = (data << 8) + c;
	r &= GetByte(c);
	data = (data << 8) + c;
	r &= GetByte(c);
	data = (data << 8) + c;
	r &= GetByte(c);
	data = (data << 8) + c;

	return r;
}

// --------------------------------------------------------------------------------------------
bool BitStream::Get2Bytes(uint32 &data)
{
	uint32 c;
	bool r;

	r = GetByte(c);
	data = c;
	r &= GetByte(c);
	data = (data << 8) + c;

	return r;
}

#endif