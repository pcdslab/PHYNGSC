/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
*/

#include "defs.h"
#include "huffman.h"
#include "utils.h"
#include "bit_stream.h"
#include <algorithm>


// --------------------------------------------------------------------------------------------
HuffmanEncoder::Code* HuffmanEncoder::Complete(bool compact)
{
	if (!n_symbols)
		return NULL;

	// Make heap of symbols
	std::make_heap(heap, heap+n_symbols);
	
	// Prepare leaves of the tree
	for (uint32 i = 0; i < n_symbols; ++i)
	{
		codes[i].code = 0;
		codes[i].len  = 0;
		tree[i].left_child = -1;
		tree[i].right_child = -1;
	}
	for (uint32 i = n_symbols; i < 2*n_symbols-1; ++i)
	{
		codes[i].code = 0;
		codes[i].len  = 0;
	}

	// Build tree
	int32 heap_size = n_symbols;

	// Remove symbols with 0 frequency
	if (compact)
	{
		while (heap_size > 2 && heap[0].frequency == 0)
		{
			std::pop_heap(heap, heap+heap_size--);
		}
	}

	int32 present_symbols = heap_size;

	if (!present_symbols)
		return codes;

	for (int32 i = 0; i < present_symbols-1; ++i)
	{
		Frequency left = heap[0];
		std::pop_heap(heap, heap+heap_size--);
		Frequency right = heap[0];
		std::pop_heap(heap, heap+heap_size--);
		
		heap[heap_size].symbol = n_symbols+i;
		heap[heap_size].frequency = left.frequency + right.frequency;
		std::push_heap(heap, heap+ ++heap_size);

		tree[n_symbols+i].left_child  = left.symbol;
		tree[n_symbols+i].right_child = right.symbol;
	}

	// Compute codes
	for (uint32 i = n_symbols+present_symbols-2; i >= n_symbols; --i)
	{
		codes[tree[i].left_child].len   = codes[i].len+1;
		codes[tree[i].left_child].code  = (codes[i].code << 1);
		codes[tree[i].right_child].len  = codes[i].len+1;
		codes[tree[i].right_child].code = (codes[i].code << 1) | 1;
	}

	root_id = n_symbols + present_symbols - 2;
	cur_id = root_id;

	return codes;
} 

// --------------------------------------------------------------------------------------------
bool HuffmanEncoder::StoreTree(uchar *&mem, uint32 &len)
{
	if (bit_memory)
		delete bit_memory;

	bit_memory = new BitMemory();

	const uint32 tree_size = n_symbols;		// in prev version 'size' was used
	bits_per_id = utils::int_log(tree_size, 2);
	if (tree_size & (tree_size-1))			// size is not power of 2
		bits_per_id++;

	min_len = n_symbols;					 // was 32 here and failed for n_symbols == 1
	for (uint32 i = 0; i < n_symbols; ++i)
	{
		if (codes[i].len < min_len && codes[i].len > 0)
			min_len = codes[i].len;
	}

	int32 node_id = root_id;
	bit_memory->PutWord(root_id);
	bit_memory->PutWord(n_symbols);
	bit_memory->PutByte((uchar) min_len);
	EncodeProcess(node_id);
	bit_memory->FlushPartialWordBuffer();
	
	mem = bit_memory->GetMemoryPtr();
	len = (int32) bit_memory->GetPos();

	return true;
}


// --------------------------------------------------------------------------------------------
bool HuffmanEncoder::LoadTree(uchar *mem, uint32 len)
{
	if (bit_memory)
		delete bit_memory;

	bit_memory = new BitMemory();
	bit_memory->Open(mem, len);

	uint32 tmp;  

	bit_memory->GetWord(tmp);
	root_id = tmp;
	bit_memory->GetWord(tmp);
	my_assert(tmp < 1 << 20);
	n_symbols = tmp;

	tmp_id = root_id;
	cur_id = root_id;

	bit_memory->GetByte(min_len);
	RestartDecompress(n_symbols, root_id);
	bits_per_id = utils::int_log(size, 2);

	if (size & (size-1))			// size is not power of 2
		bits_per_id++;

//	n_symbols = tmp;
	int32 node_id = root_id - n_symbols + 1;
	root_id = node_id;
	tmp_id = root_id;
	cur_id = root_id;
	DecodeProcess(node_id);

	if (!min_len)
		min_len = 1;
	ComputeSpeedupTree();

	delete bit_memory;
	bit_memory = NULL;

	return true;
}

// --------------------------------------------------------------------------------------------
void HuffmanEncoder::ComputeSpeedupTree()
{
	if (!min_len)
		return;

	if (speedup_tree)
		delete[] speedup_tree;
	speedup_tree = new int32[(uint32) (1 << min_len)];

	for (int32 i = 0; i < (1 << min_len); ++i)
	{
		cur_id = root_id;
		for (int32 j = min_len-1; j >= 0; --j)
		{
			Decode(i & (1 << j));
		}
		speedup_tree[i] = cur_id;
	}

	cur_id = root_id;
	tmp_id = root_id;
}


// --------------------------------------------------------------------------------------------
void HuffmanEncoder::StoreTree(BitStream& bit_stream, HuffmanEncoder& tree)
{
	static uchar* mem_buf = NULL;
	static uint32 mem_size;

	tree.StoreTree(mem_buf, mem_size);
	my_assert(mem_size > 0);
	my_assert(mem_buf);

	bit_stream.FlushPartialWordBuffer();
	bit_stream.PutWord(mem_size);
	bit_stream.PutBytes(mem_buf, mem_size);

	delete [] mem_buf;
}
