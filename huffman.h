/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
*/

#ifndef _HUFFMAN_H
#define _HUFFMAN_H

#include "defs.h"
#include "bit_memory.h"

// --------------------------------------------------------------------------------------------
//
// --------------------------------------------------------------------------------------------
class HuffmanEncoder 
{
public:
	struct Code
	{
		uint32 code;
		uint32 len;
	};

	inline HuffmanEncoder(uint32 _size = 0);
	inline ~HuffmanEncoder();

	static void StoreTree(BitStream& bit_stream, HuffmanEncoder& tree);
	static void LoadTree(BitStream& bit_stream, HuffmanEncoder& tree);

	inline bool Insert(const uint32 frequency);
	Code* Complete(bool compact = true);
	inline int32 Decode(const uint32 bit);

	inline bool Restart(uint32 _size = 0);
	inline bool RestartDecompress(uint32 _size = 0, uint32 _root_id = 0);
	inline int32 DecodeFast(const uint32 bits);

	bool StoreTree(uchar *&mem, uint32 &len);
	bool LoadTree(uchar *mem, uint32 len);

	uint32 GetMinLen() const { return min_len; }
	uint32 GetBitsPerId() const { return bits_per_id; }
	uint32 GerSymbolsNum() const { return n_symbols; }

	const Code* GetCodes() const {return codes;}

private:
	struct Frequency 
	{
		uint32 symbol;
		uint32 frequency;

		bool operator<(const Frequency &x) const
		{
			return frequency > x.frequency || (frequency == x.frequency && symbol > x.symbol);
		}
	};

	struct Node
	{
		int32 left_child;
		int32 right_child;
	};

	uint32 size;
	uint32 n_symbols;
	uint32 min_len;

	int32 root_id;
	int32 cur_id;
	int32 tmp_id;

	uint32 bits_per_id;
	BitMemory *bit_memory;
	int32 *speedup_tree;
	Node *tree;
	Code *codes;
	Frequency *heap;

	inline void EncodeProcess(int32 node_id);
	inline int32 DecodeProcess(int32 node_id);

	void ComputeSpeedupTree();
};

// --------------------------------------------------------------------------------------------
HuffmanEncoder::HuffmanEncoder(uint32 _size)
	:	size(_size)
	,	n_symbols(0)
	,	min_len(1)
	,	root_id(0)
	,	cur_id(0)
	,	tmp_id(0)
	,	bits_per_id(0)
	,	bit_memory(NULL)
	,	speedup_tree(NULL)
{
	if (size)
	{
		tree  = new Node[2*size-1];
		codes = new Code[2*size-1];
		heap  = new Frequency[size];
	}
	else
	{
		tree  = NULL;
		codes = NULL;
		heap  = NULL;
	}
}

// --------------------------------------------------------------------------------------------
HuffmanEncoder::~HuffmanEncoder()
{
	if (tree)
		delete[] tree;
	if (codes)
		delete[] codes;
	if (heap)
		delete[] heap;

	if (speedup_tree)
		delete[] speedup_tree;

	if (bit_memory)
		delete bit_memory;
}

// --------------------------------------------------------------------------------------------
void HuffmanEncoder::EncodeProcess(int32 node_id)
{
	if (tree[node_id].left_child == -1)		// leaf
	{
		bit_memory->PutBit(1);
		bit_memory->PutBits(node_id, bits_per_id);
	}
	else
	{
		bit_memory->PutBit(0);
		EncodeProcess(tree[node_id].left_child);
		EncodeProcess(tree[node_id].right_child);
	}
}

// --------------------------------------------------------------------------------------------
int32 HuffmanEncoder::DecodeProcess(int32 node_id)
{
	my_assert(node_id >= 0);

	uint32 flag(0);
	uint32 tmp;

	bit_memory->GetBit(flag);

	if (!flag)
	{
		--tmp_id;
		tree[node_id].left_child  = DecodeProcess(tmp_id);
		tree[node_id].right_child = DecodeProcess(tmp_id);
		return node_id;
	}
	else
	{
		bit_memory->GetBits(tmp, bits_per_id);
		//		tree[tmp].left_child  = -1;
		//		tree[tmp].right_child = -1;
		return -(int32)tmp;
	}
}

// --------------------------------------------------------------------------------------------
bool HuffmanEncoder::Insert(const uint32 frequency)
{
	if (n_symbols == size)
		return false;

	heap[n_symbols].symbol    = n_symbols;
	heap[n_symbols].frequency = frequency;
	n_symbols++;

	return true;
}

// --------------------------------------------------------------------------------------------
int32 HuffmanEncoder::Decode(const uint32 bit) 
{
	if (cur_id <= 0)
		cur_id = root_id;
	if (bit)
		cur_id = tree[cur_id].right_child;
	else
		cur_id = tree[cur_id].left_child;

	if (cur_id <= 0)
		return -cur_id;				// Symbol found
	
	return -1;		
}

// --------------------------------------------------------------------------------------------
inline int32 HuffmanEncoder::DecodeFast(const uint32 bits)
{
	cur_id = speedup_tree[bits];
//	if (cur_id < n_symbols)
	if (cur_id <= 0)
		return -cur_id;				// Symbol found
	
	return -1;					// Not found yet
}


// --------------------------------------------------------------------------------------------
bool HuffmanEncoder::Restart(uint32 _size)
{
	if (size != _size)
	{
		if (tree)
			delete[] tree;
		if (codes)
			delete[] codes;
		if (heap)
			delete[] heap;
		if (speedup_tree)
			delete[] speedup_tree;

		if (size)
		{
			tree  = new Node[2*size-1];
			codes = new Code[2*size-1];
			heap  = new Frequency[size];
		}
		else
		{
			tree  = NULL;
			codes = NULL;
			heap  = NULL;
		}
		speedup_tree = NULL;
	}

	n_symbols = 0;

	return true;
}

// --------------------------------------------------------------------------------------------
bool HuffmanEncoder::RestartDecompress(uint32 _size, uint32 _root_id)
{
	if (tree)
		delete[] tree;
	if (codes)
		delete[] codes;
	if (heap)
		delete[] heap;
	if (speedup_tree)
		delete[] speedup_tree;

	size = _size;
	if (size)
	{
		uint32 tsize = _root_id - _size + 2;
		tree  = new Node[tsize];
		std::fill((uchar*)tree, (uchar*)tree + tsize*sizeof(Node), 0);
		//		codes = new Code[2*size-1];
		codes = NULL;
		//		heap  = new Frequency[size];
		heap = NULL;
	}
	else
	{
		tree  = NULL;
		codes = NULL;
		heap  = NULL;
	}
	n_symbols = size;

	speedup_tree = NULL;

	return true;
}


#endif