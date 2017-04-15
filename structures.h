/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Sebastian Deorowicz and Szymon Grabowski
  Contributions: Lucas Roguski
  
  Version: 1.00
*/

#ifndef _STRUCTURES_H
#define _STRUCTURES_H

#include "defs.h"
#include "huffman.h"
#include "utils.h"
#include <map>
#include <vector>

// --------------------------------------------------------------------------------------------
// Field structure for the analysis and storage of the title portion of a FASTQ record
// --------------------------------------------------------------------------------------------
struct Field 
{
	static const uint32 HUF_GLOBAL_SIZE = 512;
	static const uint32 HUF_LOCAL_SIZE = 256;

	struct BlockDesc 
	{
		bool is_block_constant;
		bool is_block_value_constant;
		bool is_block_delta_constant;
		uint32 block_delta_constant;
		inline BlockDesc(bool _is_block_constant = true, bool _is_block_value_constant = true, 
			bool _is_block_delta_constant = true, uint32 _block_delta_constant = 0);
	};

	uint32 len;
	uint32 min_len;
	uint32 max_len;
	uchar sep;
	bool is_constant;
	bool is_len_constant;
	bool is_numeric;
	int32 min_value;
	int32 max_value;
	int32 min_delta;
	int32 max_delta;
	uint32 no_of_bits_per_num;
	uint32 no_of_bits_per_value;
	uint32 no_of_bits_per_len;
	bool is_delta_coding;

	uint32 block_str_start;
	uint32 block_str_len;
	int32 block_value;
	int32 block_delta;

	std::vector<BlockDesc> block_desc;

	uchar *data;
	uint32 *global_stats;
	uint32 **stats;
	uint32 *raw_stats;
	bool *Ham_mask;
	HuffmanEncoder *Huffman_global;
	std::vector<HuffmanEncoder*> Huffman_local;

	std::map<int32, int32> num_values;
	std::map<int32, int32> delta_values;
	std::vector<std::map<uint32, uint32> > chars;

	inline Field();
	inline Field(const Field &x);
	inline ~Field();

	bool IsNum() const { return utils::is_num(data, len); }
	uint32 ToNum() const { return utils::to_num(data, len); }
	uint32 ToString() const { return utils::to_string(data, len); }
};


// --------------------------------------------------------------------------------------------
Field::BlockDesc::BlockDesc(bool _is_block_constant, bool _is_block_value_constant, 
							bool _is_block_delta_constant, uint32 _block_delta_constant)
	:	is_block_constant(_is_block_constant)
	,	is_block_value_constant(_is_block_value_constant)
	,	is_block_delta_constant(_is_block_delta_constant)
	,	block_delta_constant(_block_delta_constant)
{

}

// --------------------------------------------------------------------------------------------
Field::Field()
	:	len(0)
	,	min_len(0)
	,	max_len(0)
	,	sep(0)
	,	is_constant(false)
	,	is_len_constant(false)
	,	is_numeric(false)
	,	min_value(1)
	,	max_value(-1)
	,	min_delta(1)
	,	max_delta(-1)
	,	no_of_bits_per_num(0)
	,	no_of_bits_per_value(0)
	,	no_of_bits_per_len(0)
	,	is_delta_coding(false)
	,	block_str_start(0)
	,	block_str_len(0)
	,	block_value(0)
	,	block_delta(0)
	,	data(NULL)
	,	global_stats(NULL)
	,	stats(NULL)
	,	raw_stats(NULL)
	,	Ham_mask(NULL)
	,	Huffman_global(NULL)
{
	
}

// --------------------------------------------------------------------------------------------
Field::Field(const Field &y)
	:	len(y.len)
	,	min_len(y.min_len)
	,	max_len(y.max_len)
	,	sep(y.sep)
	,	is_constant(y.is_constant)
	,	is_len_constant(y.is_len_constant)
	,	is_numeric(y.is_numeric)
	,	min_value(y.min_value)
	,	max_value(y.max_value)
	,	min_delta(y.min_delta)
	,	max_delta(y.max_delta)
	,	no_of_bits_per_num(y.no_of_bits_per_num)
	,	no_of_bits_per_value(y.no_of_bits_per_value)
	,	no_of_bits_per_len(y.no_of_bits_per_len)
	,	is_delta_coding(y.is_delta_coding)
	,	block_str_start(y.block_str_start)
	,	block_str_len(y.block_str_len)
	,	block_value(y.block_value)
	,	block_delta(y.block_delta)
	,	block_desc(y.block_desc)
{
	Field &x = const_cast<Field &>(y);

	if (x.data && len)
	{
		data   = x.data;
		x.data = NULL;
	}
	else
	{
		data = NULL;
	}

	if (x.Ham_mask && len)
	{
		Ham_mask   = x.Ham_mask;
		x.Ham_mask = NULL;
	}
	else
	{
		Ham_mask = NULL;
	}

	if (x.stats)
	{
		stats       = x.stats;
		raw_stats   = x.raw_stats;
		x.stats     = NULL;
		x.raw_stats = NULL;
	}
	else
	{
		stats     = NULL;
		raw_stats = NULL;
	}

	if (x.global_stats)
	{
		global_stats = x.global_stats;
		x.global_stats = NULL;
	}
	else
	{
		global_stats = NULL;
	}

	if (x.Huffman_global)
	{
		Huffman_global = x.Huffman_global;
		x.Huffman_global = new HuffmanEncoder(Field::HUF_GLOBAL_SIZE);
	}
	else
	{
		Huffman_global = NULL;
	}
}

// --------------------------------------------------------------------------------------------
Field::~Field()
{
	if (data)
		delete[] data;
	if (Ham_mask)
		delete[] Ham_mask;

	if (stats)
		delete[] stats;
	if (global_stats)
		delete[] global_stats;
	if (raw_stats)
		delete[] raw_stats;

	if (Huffman_global)
		delete Huffman_global;

	for (uint32 i = 0; i < Huffman_local.size(); ++i)
	{
		if (Huffman_local[i])
			delete Huffman_local[i];
	}
	Huffman_local.resize(0);
}

// --------------------------------------------------------------------------------------------
// Record structure to identify title and fields of title, sequence and quality positions
// --------------------------------------------------------------------------------------------
struct Record 
{
	uint32 title_end_pos;
	uint32 seq_end_pos;
	int32 seq_len;
	int32 qua_len;
	std::vector<uint32> fields_end_pos;

	std::vector<uchar> title;
	std::vector<uchar> dna_seq;
	std::vector<uchar> quality;
	uchar plus = '+';

	inline Record();
	inline bool AppendData(uchar *str, uint32 len, uchar rec_line);
	inline bool AppendChar(uchar c, uchar rec_line);
};

// --------------------------------------------------------------------------------------------
Record::Record()
	:	title_end_pos(0)
	,	seq_end_pos(0)
	,	seq_len(0)
	,	qua_len(0)
{
	
}

// --------------------------------------------------------------------------------------------
bool Record::AppendData(uchar *str, uint32 len, uchar rec_line)
{
	switch (rec_line)
	{
		case 'T':
			for (uint32 i = 0; i < len; ++i)
				title.push_back(str[i]);
			break;
	  	case 'D':
	  		for (uint32 i = 0; i < len; ++i)
				dna_seq.push_back(str[i]);
			break;
		case 'Q':
	  		for (uint32 i = 0; i < len; ++i)
				quality.push_back(str[i]);
			break;
	}

	return true;
}

// --------------------------------------------------------------------------------------------
bool Record::AppendChar(uchar c, uchar rec_line)
{
	switch (rec_line)
	{
		case 'T': title.push_back(c); break;
	  	case 'D': dna_seq.push_back(c); break;
		case 'Q': quality.push_back(c); break;
	}

	return true;
}

// --------------------------------------------------------------------------------------------
// Processor's compression information structure
// --------------------------------------------------------------------------------------------
struct ProcessCompressionInfo
{
	int32 n_blocks;
  	int32 n_subblocks;
  	uint32 last_block_size;
  	int32 wr_overlap;
};

// --------------------------------------------------------------------------------------------
// Block header structure to store header information
// --------------------------------------------------------------------------------------------
struct BlockHeader
{
  	int32 WRID;				 // Working Region ID
  	int32 BHS;				 // Block Header Size
  	int32 BEWR;     		 // Bits Encoding Working Region ID   
  	int32 BESO;     		 // Bits Encoding SubBlock Offsets
  	int32 NOSB; 			 // Number Of SubBlocks.
  	short BCSS;  			 // Block Contains Splitted SubBlocks?
  	std::vector<int32> SBOL; // SubBlock's Offset list

  	inline int32 HeaderSize();
};

int32 BlockHeader::HeaderSize()
{
	// 6 bits to encode NOSB and 12 bits to encode BHS
    double size = BEWR + 6 + 12;
    // 5 bits and 2 bits to encode BESO and BCSS
    size += (BESO * SBOL.size()) + 5 + 2;

    // Return size of header in bytes
    return ceil(size / 8);

}

// --------------------------------------------------------------------------------------------
// Structure to store footer information
// --------------------------------------------------------------------------------------------
struct Footer
{
  int32 BEPS;
  int32 BEFS;           
  int32 BEBS;    
  int32 BESS;     
  int32 BELB;
  int32 BEOV;   
  short LBES;    
  int32 PS;         
  int64 FS;     
  int32 BS;       
  int32 SS;         
  std::vector<int32>  OV;        
  std::vector<int32>  CBO;      
  std::vector<uint32> LBS;
  std::vector<uint32> ABS;
};

#endif