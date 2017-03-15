/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
*/

#ifndef _BIT_MEMORY_H
#define _BIT_MEMORY_H

#include "defs.h"
#include "utils.h"
#include <vector>
#include <algorithm>

// --------------------------------------------------------------------------------------------
//
// --------------------------------------------------------------------------------------------
class BitMemory 
{
public:
	inline BitMemory();
	inline BitMemory(const BitMemory &y); 
	inline ~BitMemory();

	inline bool FlushFullWordBuffer();
	inline bool FlushPartialWordBuffer();
	inline bool FlushInputWordBuffer();

	inline bool Open(uchar *p, int64 size, bool force_open = false);
	inline bool Create(int64 size = 1);
	inline bool Complete();
	inline bool Close();
	inline bool Restart();
	inline bool TakeOwnership();

	uint64 GetPos() const {return mem_buffer_pos;};
	uchar* GetMemoryPtr() const {return mem_buffer;};
	inline bool SetPos(int64 pos);

	inline bool PutBit(const uint32 word);
	inline bool Put2Bits(const uint32 word);
	inline bool PutBits(uint32 word, int32 n_bits);
	inline bool PutBytes(const uchar *data, int32 n_bytes);
	inline bool PutByte(const uchar byte);
	inline bool PutNBytes(const uchar byte, int32 n_bytes);
	inline bool PutWord(const uint32 data);
	inline bool PutDWord(const uint64 data);
	inline bool Put2Bytes(const uint32 data);

	inline bool GetBit(uint32 &word);
	inline bool Get2Bits(uint32 &word);
	inline bool GetBits(uint32 &word, uint32 n_bits);
	inline bool GetBytes(uchar *data, uint32 n_bytes);
	inline bool GetBool(bool &byte);
	inline bool GetByte(uint32 &byte);
	inline bool GetWord(uint32 &data);
	inline bool GetWord(int32 &data);
	inline bool GetDWord(uint64 &data);
	inline bool Get2Bytes(uint32 &data);

	inline bool SkipBytes(uint32 n_bytes);

	inline uint32 BitLength(const uint64 x);
	uint64 GetFilePos() { return mem_buffer_pos; };

private:
	uchar *mem_buffer;
	bool mem_buffer_ownership;
	int64 mem_buffer_pos;
	int64 mem_buffer_size;
	uint32 word_buffer;
	int32 word_buffer_pos;
	int32 word_buffer_size;
	MemoryModeEnum mode;
	uint32* n_bit_mask;

	template <typename _Ty, uint32 _BitNum>
	struct TBitMask
	{
		static _Ty* Get()
		{
			static _Ty bit_mask[_BitNum];
			for (_Ty i = 0; i < _BitNum; ++i)
			{
				bit_mask[i] = (1u << i) - 1;
			}
			return bit_mask;
		}
	};
};

// --------------------------------------------------------------------------------------------
BitMemory::BitMemory()
	:	mem_buffer(NULL)
	,	mem_buffer_ownership(false)
	,	mem_buffer_pos(0)
	,	mem_buffer_size(0)
	,	word_buffer(0)
	,	word_buffer_pos(0)
	,	word_buffer_size(32)
	,	mode(MEM_MODE_NONE)
{
	n_bit_mask = TBitMask<uint32, 32>::Get();
}

// --------------------------------------------------------------------------------------------
BitMemory::BitMemory(const BitMemory &y)
{
	// rebind pointers
	BitMemory &x = const_cast<BitMemory &>(y);

	mem_buffer = x.mem_buffer;
	x.mem_buffer = NULL;
	mem_buffer_ownership = x.mem_buffer_ownership;
	x.mem_buffer_ownership = false;

	mem_buffer_pos = x.mem_buffer_pos;
	x.mem_buffer_pos = 0;

	mem_buffer_size = x.mem_buffer_size;
	x.mem_buffer_size = 0;

	word_buffer = x.word_buffer;
	x.word_buffer = 0;

	word_buffer_pos = x.word_buffer_pos;
	x.word_buffer_pos = 0;

	word_buffer_size = x.word_buffer_size;
	x.word_buffer_size = 0;

	mode = x.mode;
	x.mode = MEM_MODE_NONE;

	n_bit_mask = TBitMask<uint32, 32>::Get();
}

// --------------------------------------------------------------------------------------------
BitMemory::~BitMemory()
{
	if (mem_buffer && mem_buffer_ownership)
		delete[] mem_buffer;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::PutBit(const uint32 word)
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
bool BitMemory::Put2Bits(const uint32 word)
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
bool BitMemory::PutBits(uint32 word, int32 n_bits)
{
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
	
	word_buffer     <<= n_bits;
	word_buffer     += word;
	word_buffer_pos += n_bits;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::FlushFullWordBuffer()
{
	PutWord(word_buffer);
	
	word_buffer     = 0;
	word_buffer_pos = 0;
	
	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::FlushPartialWordBuffer()
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
	
	word_buffer     = 0;
	word_buffer_pos = 0;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::FlushInputWordBuffer()
{
	word_buffer_pos = 0;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::PutByte(const uchar byte)
{
	if (mem_buffer_pos + 1 > mem_buffer_size)
	{
		mem_buffer_size = (uint64) ((mem_buffer_pos + 1) * 1.5);
		uchar *new_mem_buffer = new uchar[mem_buffer_size];
		if (mem_buffer)
		{
			copy_n(mem_buffer, mem_buffer_pos, new_mem_buffer);
			delete[] mem_buffer;
		}
		mem_buffer = new_mem_buffer;
	}

	mem_buffer[mem_buffer_pos++] = byte;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::Put2Bytes(const uint32 data)
{
	PutByte((uchar) (data >> 8));
	PutByte((uchar) (data & 0xFF));

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::PutBytes(const uchar *data, int32 n_bytes)
{
	if (mem_buffer_pos + n_bytes > mem_buffer_size)
	{
		mem_buffer_size = (uint64) ((mem_buffer_pos + n_bytes) * 1.5);
		uchar *new_mem_buffer = new uchar[mem_buffer_size];
		if (mem_buffer)
		{
			copy_n(mem_buffer, mem_buffer_pos, new_mem_buffer);
			delete[] mem_buffer;
		}
		mem_buffer = new_mem_buffer;
	}
	copy_n(data, n_bytes, mem_buffer+mem_buffer_pos);
	mem_buffer_pos += n_bytes;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::PutNBytes(const uchar byte, int32 n_bytes)
{
	if (mem_buffer_pos + n_bytes > mem_buffer_size)
	{
		mem_buffer_size = (uint64) ((mem_buffer_pos + n_bytes) * 1.5);
		uchar *new_mem_buffer = new uchar[mem_buffer_size];
		if (mem_buffer)
		{
			copy_n(mem_buffer, mem_buffer_pos, new_mem_buffer);
			delete[] mem_buffer;
		}
		mem_buffer = new_mem_buffer;
	}
	std::fill_n(mem_buffer+mem_buffer_pos, n_bytes, byte);
	mem_buffer_pos += n_bytes;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::PutDWord(const uint64 data)
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
bool BitMemory::PutWord(const uint32 data)
{
	PutByte(data >> 24);
	PutByte((data >> 16) & 0xFF);
	PutByte((data >> 8) & 0xFF);
	PutByte(data & 0xFF);

	return true;
}
// --------------------------------------------------------------------------------------------
uint32 BitMemory::BitLength(const uint64 x)
{
	for (uint32 i = 0; i < 32; ++i)
	{
		if (x < (1ull << i))
			return i;
	}
	return 64;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::GetBit(uint32 &word)
{
	if (word_buffer_pos == 0)
	{
		if (!GetByte(word_buffer))
			return false;
		word_buffer_pos = 7;
		word = word_buffer >> 7;
	}
	else
		word = (word_buffer >> (--word_buffer_pos)) & 1;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::Get2Bits(uint32 &word)
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
bool BitMemory::GetBits(uint32 &word, uint32 n_bits)
{
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
bool BitMemory::GetByte(uint32 &byte)
{
	if (mem_buffer_pos >= mem_buffer_size)
	{
		byte = (uint32)INVALID_DWORD;
		return false;
	}
	byte = mem_buffer[mem_buffer_pos++];

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::GetBool(bool &byte)
{
	if (mem_buffer_pos >= mem_buffer_size)
	{
		byte = false;
		return false;
	}
	byte = (mem_buffer[mem_buffer_pos++] == 1);

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::GetBytes(uchar *data, uint32 n_bytes)
{
	if (mem_buffer_pos + n_bytes >= mem_buffer_size)
	{
		std::fill(data, data+n_bytes, INVALID_BYTE);
		return false;
	}

	copy_n(mem_buffer+mem_buffer_pos, n_bytes, data);
	mem_buffer_pos += n_bytes;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::GetWord(uint32 &data)
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
bool BitMemory::GetWord(int32 &data)
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
bool BitMemory::GetDWord(uint64 &data)
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
bool BitMemory::Get2Bytes(uint32 &data)
{
	uint32 c;
	bool r;

	r = GetByte(c);
	data = c;
	r &= GetByte(c);
	data = (data << 8) + c;

	return r;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::SkipBytes(uint32 n_bytes)
{
	if (mem_buffer_pos + n_bytes >= mem_buffer_size)
		return false;

	mem_buffer_pos += n_bytes;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::TakeOwnership()
{
	if (mem_buffer && mem_buffer_ownership)
	{
		mem_buffer_ownership = false;	
		return true;
	}

	return false;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::Open(uchar *p, int64 size, bool force_open)
{
	if (!force_open)
	{
		if (mode != MEM_MODE_NONE)
			return false;
	}

	if (mem_buffer && mem_buffer_ownership)
		delete[] mem_buffer;

	mem_buffer_size = size;
/*	if (!mem_buffer_size)
		mem_buffer_size = 1;
	mem_buffer = new uchar[mem_buffer_size];
	copy_n(p, size, mem_buffer);*/
	mem_buffer = p;
	mem_buffer_ownership = false;

	if (!mem_buffer_size)
	{
		mem_buffer_size = 1;
		mem_buffer = new uchar[mem_buffer_size];
		mem_buffer_ownership = true;
		copy_n(p, size, mem_buffer);
	}

	mode = MEM_MODE_READ;
	
	word_buffer_size = 8;
	mem_buffer_pos   = 0;	

	return mode == MEM_MODE_READ;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::Restart()
{
	mode = MEM_MODE_READ;
	
	word_buffer_size = 8;
	mem_buffer_pos   = 0;
	word_buffer_pos  = 0;
	word_buffer		 = 0;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::Create(int64 size)
{
	if (mode != MEM_MODE_NONE)
		return false;

	if (mem_buffer && mem_buffer_ownership)
		delete[] mem_buffer;

	if (!size)
		size = 1;
	mem_buffer_size = size;
	mem_buffer = new uchar[mem_buffer_size];
	mem_buffer_ownership = true;

	mode = MEM_MODE_WRITE;
	mem_buffer_pos = 0;

	word_buffer_size = 32;

	return mode == MEM_MODE_WRITE;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::Close()
{
	if (mode != MEM_MODE_WRITE && mode != MEM_MODE_READ)
		return false;

	if (mode == MEM_MODE_WRITE)
	{
		if (word_buffer_pos)
			FlushPartialWordBuffer();
	}

	if (mem_buffer && mem_buffer_ownership)
		delete[] mem_buffer;

	mem_buffer = NULL;
	mem_buffer_ownership = false;
	mode = MEM_MODE_NONE;

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::Complete()
{
	if (mode != MEM_MODE_WRITE && mode != MEM_MODE_READ)
		return false;

	if (mode == MEM_MODE_WRITE)
	{
		if (word_buffer_pos)
			FlushPartialWordBuffer();
	}

	return true;
}

// --------------------------------------------------------------------------------------------
bool BitMemory::SetPos(int64 pos)
{
	if (mode != MEM_MODE_READ)
		return false;

	if ((int64) pos > mem_buffer_size)
		return false;
	mem_buffer_pos = pos;

	word_buffer_pos = 0;
	word_buffer     = 0;

	return true;
}

#endif
