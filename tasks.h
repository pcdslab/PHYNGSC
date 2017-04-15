/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
*/

#ifndef _TASKS_H
#define _TASKS_H

#include <vector>
#include <map>
#include "defs.h"
#include "huffman.h"
#include "structures.h"

// --------------------------------------------------------------------------------------------
void AnalyzeTitleFields(uchar *read_buffer, uint32 field_0_rec_0_pos, int32 c_field, 
  std::vector<Field> &fields, Record *&records, int64 no_records);

// --------------------------------------------------------------------------------------------
void AnalyzeDNA(uint32 &fastq_flags, std::vector<uint32> &dna_occ, std::vector<uchar> &symbols,
  uint32 *&sym_stats);

// --------------------------------------------------------------------------------------------
void AnalyzeQuality(uchar *read_buffer, Record *&records, int64 no_records, uint32 max_quality_length, 
  std::vector<uchar> &qualities, uchar *qua_code, uint32 **&quality_stats, uint32 *&raw_quality_stats);

// --------------------------------------------------------------------------------------------
void StoreTitle(BitStream &title_bit_stream, std::vector<Field> &fields, uchar *read_buffer,
  Record *&records, uint32 field_0_rec_0_pos, int64 no_records);

// --------------------------------------------------------------------------------------------
void StoreDNA(BitStream &dna_bit_stream, uchar *read_buffer, Record *&records, int64 no_records, 
  uint32 fastq_flags, std::vector<uchar> &symbols, uint32 *&sym_stats, uchar *sym_code);

// --------------------------------------------------------------------------------------------
void StoreQuality(BitStream &quality_bit_stream, uchar *read_buffer, Record *&records, int64 no_records, 
  uint32 max_quality_length, uchar *qua_code, std::vector<uchar> &qualities, uint32 **&quality_stats, 
  HuffmanEncoder::Code **&qua_huf_codes, HuffmanEncoder::Code *&raw_qua_huf_codes);

// --------------------------------------------------------------------------------------------
void FetchTitleHeader(BitStream &sb_bit_stream, std::vector<Field> &fields);

// --------------------------------------------------------------------------------------------
void FetchTitleBody(BitStream &sb_bit_stream, std::vector<Field> &fields, Record *&records, 
  int64 no_records, uint32 fastq_flags);

// --------------------------------------------------------------------------------------------
void FetchDNA(BitStream &sb_bit_stream, std::vector<uchar> &symbols, uint32 no_symbols, uint32 fastq_flags,
  Record *&records, int64 no_records, std::vector<uint32> &no_ambiguity);

// --------------------------------------------------------------------------------------------
void FetchQuality(BitStream &sb_bit_stream, std::vector<uchar> &qualities, uint32 no_qualities, 
  Record *&records, int64 no_records, uint32 fastq_flags, std::vector<uint32> &no_ambiguity,
  uint32 max_quality_length);

// --------------------------------------------------------------------------------------------
void MakeFooter(BitStream &footer_bit_stream, int32 g_size, int64 FASTQ_size, int32 n_blocks, 
  int32 n_subblocks, int32 *all_wr_overlaps, std::map<double,int32> &blocks_order, uint32 *lb_sizes);

// --------------------------------------------------------------------------------------------
void MakeHeader(BitStream &header_bit_stream, BlockHeader &p_block_header);

// --------------------------------------------------------------------------------------------
void ReadFooter(BitStream &footer_bit_stream, Footer &footer) ;

#endif