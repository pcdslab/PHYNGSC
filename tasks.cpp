/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
*/

#include <vector>
#include <map>
#include <math.h>
#include <algorithm>
#include "defs.h"
#include "utils.h"
#include "structures.h"
#include "bit_stream.h"
#include "huffman.h"

// --------------------------------------------------------------------------------------------
void AnalyzeTitleFields(uchar *read_buffer, uint32 field_0_rec_0_pos, int32 c_field, 
  std::vector<Field> &fields, Record *&records, int64 no_records)
{
  const uint32 MAX_FIELD_STAT_LEN  = 128;
  const uint32 DEFAULT_B_SIZE      = 32;
  const uint32 MAX_NUM_VAL_HUF     = 512;

  int32 prev_value = 0;

  for (uint32 i = 0; i < no_records; ++i)
  {
    uint32 start_pos = 0; 
    uint32 end_pos = records[i].fields_end_pos[c_field];

    if (c_field) // Other than the field 0 assing
    {
      start_pos = records[i].fields_end_pos[c_field-1] + 1;
    }
    else
    {
      if (i)
      {
        int32 rec_const_lines = (records[i-1].seq_end_pos - records[i-1].title_end_pos) * 2;
        start_pos = records[i-1].fields_end_pos[fields.size()-1] + rec_const_lines + 3;
      }
      else
      {
        start_pos = field_0_rec_0_pos;
      }
    }

    if (end_pos - start_pos > fields[c_field].max_len)
    {
      fields[c_field].max_len = end_pos - start_pos;
    }
    else if (end_pos - start_pos < fields[c_field].min_len)
    {
      fields[c_field].min_len = end_pos - start_pos;
    }

    // Check whether field in a block is constant
    if (i % DEFAULT_B_SIZE > 0){
      if (fields[c_field].block_str_len != end_pos-start_pos)
      {
        fields[c_field].block_desc.back().is_block_constant = false;
      }
      else if (fields[c_field].block_desc.back().is_block_constant)
      {
        fields[c_field].block_desc.back().is_block_constant = std::equal(read_buffer+start_pos, read_buffer+end_pos,
          read_buffer+fields[c_field].block_str_start);
      }
    }
    else
    {
      fields[c_field].block_str_start = start_pos;
      fields[c_field].block_str_len   = end_pos-start_pos;
      fields[c_field].block_desc.push_back(Field::BlockDesc(true, true, true, 0));
    }

    // Consider instead of using resize(), the use of push_back()
    // so you don't have memory allocations within a thread
    fields[c_field].chars.resize(fields[c_field].max_len);
    uint32 chars_len = MIN(MAX_FIELD_STAT_LEN, end_pos-start_pos);

    for (uint32 x = 0; x < chars_len; ++x)
    {
      fields[c_field].chars[x][read_buffer[start_pos+x]]++;
    }
    for (uint32 x = MAX_FIELD_STAT_LEN; x < end_pos-start_pos; ++x)
    {
      fields[c_field].chars[MAX_FIELD_STAT_LEN][read_buffer[start_pos+x]]++;
    }
    
    if (fields[c_field].is_constant)
    {
      if (end_pos - start_pos != fields[c_field].len)
      {
        fields[c_field].is_constant = false;
      }
      else
      {
        fields[c_field].is_constant = std::equal(fields[c_field].data, fields[c_field].data + fields[c_field].len, read_buffer+start_pos);
      }
    }

    if (fields[c_field].is_len_constant)
    {
      fields[c_field].is_len_constant = fields[c_field].len == end_pos-start_pos;
    }
    
    if (fields[c_field].is_numeric)
    {
      fields[c_field].is_numeric = utils::is_num(read_buffer+start_pos, end_pos-start_pos);
      if (fields[c_field].is_numeric)
      {
        int32 value = utils::to_num(read_buffer+start_pos, end_pos-start_pos);
        if (value < fields[c_field].min_value)
        {
          fields[c_field].min_value = value;
        }
        else if (value > fields[c_field].max_value)
        {
          fields[c_field].max_value = value;
        }

        if (i % DEFAULT_B_SIZE > 0)
        {
          if (fields[c_field].block_value != value)
            fields[c_field].block_desc.back().is_block_value_constant = false;

          // when there are more than 2 records on the new block
          if (i % DEFAULT_B_SIZE > 1)
          {
            if (value - prev_value != fields[c_field].block_delta)
              fields[c_field].block_desc.back().is_block_delta_constant = false;
          }
          else // j == 1, when there are 2 records on the new block
          {
            fields[c_field].block_delta = value - prev_value;
            fields[c_field].block_desc.back().block_delta_constant = value - prev_value;
          }
        }
        else
        {
          fields[c_field].block_value = value;
        }

        if (i >= 1)
        {
          if (i > 1)
          {
            if (value - prev_value > fields[c_field].max_delta)
            {
              fields[c_field].max_delta = value - prev_value;
            }
            if (value - prev_value < fields[c_field].min_delta)
            {
              fields[c_field].min_delta = value - prev_value;
            }
          }
          else // i == 1
          {
            fields[c_field].max_delta = value - prev_value;
            fields[c_field].min_delta = value - prev_value;
          }

          fields[c_field].delta_values[value - prev_value]++;
          if (fields[c_field].delta_values.size() > MAX_NUM_VAL_HUF)
          {
            fields[c_field].delta_values.clear();
          }
        }

        if (fields[c_field].num_values.size())
        {
          fields[c_field].num_values[value]++;
          if (fields[c_field].num_values.size() > MAX_NUM_VAL_HUF)
          {
            fields[c_field].num_values.clear();
          }
        }

        prev_value = value;
      }
    }
    if (!fields[c_field].is_constant)
    {
      for (uint32 p = 0; p < end_pos - start_pos && p < fields[c_field].len; ++p)
      {
        fields[c_field].Ham_mask[p] &= fields[c_field].data[p] == read_buffer[p+start_pos];
      }
    }
  }

    // Find better encoding of numeric values
  if (!fields[c_field].is_numeric)
  {
    if (!fields[c_field].is_constant)
    {
      fields[c_field].chars.resize(MIN(fields[c_field].max_len, MAX_FIELD_STAT_LEN+1));
      fields[c_field].no_of_bits_per_len = BitStream::BitLength(fields[c_field].max_len - fields[c_field].min_len);
    }
  } 
  else
  {
    int32 diff;
    if (fields[c_field].max_value - fields[c_field].min_value < fields[c_field].max_delta - fields[c_field].min_delta)
    {
      fields[c_field].is_delta_coding = false;
      diff = fields[c_field].max_value - fields[c_field].min_value;
    }
    else
    {
      fields[c_field].is_delta_coding = true;
      diff = fields[c_field].max_delta - fields[c_field].min_delta;
    }

    fields[c_field].no_of_bits_per_num = BitStream::BitLength(diff);
    diff = fields[c_field].max_value - fields[c_field].min_value;
    fields[c_field].no_of_bits_per_value = BitStream::BitLength(diff);
  }
}

// --------------------------------------------------------------------------------------------
void AnalyzeDNA(uint32 &fastq_flags, std::vector<uint32> &dna_occ, std::vector<uchar> &symbols,
  uint32 *&sym_stats)
{
  if (sym_stats)
  {
    delete[] sym_stats;
  }
  sym_stats = new uint32[symbols.size()];
  for (uint32 i = 0; i < symbols.size(); ++i)
  {
    sym_stats[i] = dna_occ[symbols[i]];
  }

  if (symbols.size() > 4)
  {
    fastq_flags &= ~FLAG_DNA_PLAIN;
  }
  else
  {
    std::vector<uint32> sym_tmp(dna_occ.begin(), dna_occ.end());
    std::sort(sym_tmp.begin(), sym_tmp.end());

    if (sym_tmp[0] > sym_tmp[2] + sym_tmp[3])
    {
      fastq_flags &= ~FLAG_DNA_PLAIN;
    }
    else
    {
      fastq_flags |= FLAG_DNA_PLAIN;
    }
  }
}

// --------------------------------------------------------------------------------------------
void AnalyzeQuality(uchar *read_buffer, Record *&records, int64 no_records, uint32 max_quality_length, 
  std::vector<uchar> &qualities, uchar *qua_code, uint32 **&quality_stats, uint32 *&raw_quality_stats)
{
  quality_stats = new uint32*[max_quality_length+1];
  raw_quality_stats = new uint32[(max_quality_length+1)*qualities.size()]();
  
  for (uint32 i = 0; i <= max_quality_length; ++i)
  {
    quality_stats[i] = raw_quality_stats+i*qualities.size();
  }
  for (uint32 i = 0; i < (max_quality_length+1)*qualities.size(); ++i)
  {
    raw_quality_stats[i] = 0;
  }
  
  for (uint32 i = 0; i < no_records; ++i)
  {
    int64 qua_start = records[i].seq_end_pos + 3; 
    int64 qua_end = qua_start + (records[i].seq_end_pos - records[i].title_end_pos) - 1;

    for (uint32 j = 0; qua_start < qua_end && read_buffer[qua_start] != '\0'; ++j, ++qua_start)
    {
      quality_stats[j+1][qua_code[read_buffer[qua_start]]]++;
      quality_stats[0][qua_code[read_buffer[qua_start]]]++;
    }
  }    
}

// --------------------------------------------------------------------------------------------
void StoreTitleHeader(BitStream &h_title_bit_stream, std::vector<Field> &fields)
{
  const uint32 MAX_NUM_VAL_HUF     = 512;
  const uint32 MAX_FIELD_STAT_LEN  = 128;

  int32 n_fields = fields.size();
  h_title_bit_stream.PutWord(n_fields);

  for (int32 i = 0; i < n_fields; ++i)
  {
    h_title_bit_stream.PutByte(fields[i].sep);
    h_title_bit_stream.PutByte(fields[i].is_constant);
    if (fields[i].is_constant)
    {
      h_title_bit_stream.PutWord(fields[i].len);
      h_title_bit_stream.PutBytes(fields[i].data, fields[i].len);
      continue;
    }

    h_title_bit_stream.PutByte(fields[i].is_numeric);
    if (fields[i].is_numeric)
    {
      h_title_bit_stream.PutWord(fields[i].min_value);
      h_title_bit_stream.PutWord(fields[i].max_value);
      h_title_bit_stream.PutWord(fields[i].min_delta);
      h_title_bit_stream.PutWord(fields[i].max_delta);

      int32 diff, base;
      std::map<int32, int32> &map_stats = fields[i].num_values;
      if (fields[i].max_value-fields[i].min_value < fields[i].max_delta-fields[i].min_delta)
      {
        diff = fields[i].max_value - fields[i].min_value;
        base = fields[i].min_value;
      }
      else
      {
        diff = fields[i].max_delta - fields[i].min_delta;
        base = fields[i].min_delta;
        map_stats = fields[i].delta_values;
      }

      diff++;
      if (diff <= (int32)MAX_NUM_VAL_HUF && map_stats.size())     // few values, so use Huffman for them
      {
        HuffmanEncoder* huf = fields[i].Huffman_global = new HuffmanEncoder(Field::HUF_GLOBAL_SIZE);
        for (uint32 j = 0; j < (uint32)diff; ++j)
        {
          huf->Insert(map_stats[base+j]);
        }
        huf->Complete();
        HuffmanEncoder::StoreTree(h_title_bit_stream, *huf);
      }

      continue;
    }

    h_title_bit_stream.PutByte(fields[i].is_len_constant);
    h_title_bit_stream.PutWord(fields[i].len);
    h_title_bit_stream.PutWord(fields[i].max_len);
    h_title_bit_stream.PutWord(fields[i].min_len);
    h_title_bit_stream.PutBytes(fields[i].data, fields[i].len);

    for (uint32 j = 0; j < fields[i].len; ++j)
    {
      h_title_bit_stream.PutBit(fields[i].Ham_mask[j]);
    }
    fields[i].Huffman_local.resize(MIN(fields[i].max_len+1, MAX_FIELD_STAT_LEN+1));

    for (uint32 j = 0; j < MIN(fields[i].max_len, MAX_FIELD_STAT_LEN); ++j)
    {
      fields[i].Huffman_local[j] = NULL;
      if (j >= fields[i].len || !fields[i].Ham_mask[j])
      {
        HuffmanEncoder* huf = fields[i].Huffman_local[j] = new HuffmanEncoder(Field::HUF_LOCAL_SIZE);
        for (uint32 k = 0; k < Field::HUF_LOCAL_SIZE; ++k)
        {
          huf->Insert(fields[i].chars[j][k]);
        }
        huf->Complete(true);
        HuffmanEncoder::StoreTree(h_title_bit_stream, *huf);
      }
    } 
    if (fields[i].max_len >= MAX_FIELD_STAT_LEN)
    {
      HuffmanEncoder* huf = fields[i].Huffman_local[MAX_FIELD_STAT_LEN] = new HuffmanEncoder(Field::HUF_LOCAL_SIZE);
      for (uint32 k = 0; k < Field::HUF_LOCAL_SIZE; ++k)
      {
        huf->Insert(fields[i].chars[(uchar) MAX_FIELD_STAT_LEN][k]);
      }
      huf->Complete(true);
      HuffmanEncoder::StoreTree(h_title_bit_stream, *huf);
    }

    h_title_bit_stream.FlushPartialWordBuffer();
  }
}

// --------------------------------------------------------------------------------------------
void StoreTitleBody(BitStream &b_title_bit_stream, std::vector<Field> &fields, uchar *read_buffer,
  Record *&records, uint32 field_0_rec_0_pos, int64 no_records)
{
  const int32 MAX_FIELD_STAT_LEN  = 128;
  const int32 DEFAULT_B_SIZE      = 32;

  int32 n_blocks = (no_records + DEFAULT_B_SIZE -1) / DEFAULT_B_SIZE;
  int32 n_fields = fields.size();
  std::vector<int32> prev_value(n_fields);
  int32 rec_count = 0;
  int32 first_rec_in_block;

  for (int32 block_no = 0; block_no < n_blocks; ++block_no)
  {
    first_rec_in_block = rec_count;
    rec_count += DEFAULT_B_SIZE;

    if (rec_count > no_records)
      rec_count = no_records;

    for (int32 i = 0; i < n_fields; ++i)
    {
      if (fields[i].is_constant)
        continue;

      Field::BlockDesc& block_desc = fields[i].block_desc[block_no];
      prev_value[i] = 0;
      if (!fields[i].is_numeric)
      {
        b_title_bit_stream.PutBit(block_desc.is_block_constant);
      }

      if (fields[i].is_numeric)
      {
        block_desc.is_block_delta_constant &= (int32)block_desc.block_delta_constant == fields[i].min_delta;
        if (fields[i].is_delta_coding)
        {
          b_title_bit_stream.PutBit(block_desc.is_block_delta_constant);
        }
        else
        {
          b_title_bit_stream.PutBit(block_desc.is_block_value_constant);
        }
      }   
    }

    for (int32 i = first_rec_in_block; i < rec_count; ++i)
    {
      for (int32 c_field = 0; c_field < n_fields; ++c_field)
      {
        Field &cur_field = fields[c_field];
        uint32 start_pos = 0;
        uint32 end_pos = records[i].fields_end_pos[c_field];
        
        if (c_field) // Other than the field 0 assing
          start_pos = records[i].fields_end_pos[c_field-1] + 1;
        else
        {
          if (i)
          {
            int32 rec_const_lines = (records[i-1].seq_end_pos - records[i-1].title_end_pos) * 2;
            start_pos = records[i-1].fields_end_pos[fields.size()-1] + rec_const_lines + 3;
          }
          else
            start_pos = field_0_rec_0_pos;
        }

        if (cur_field.is_constant)
          continue;

        if (cur_field.is_numeric)
        {
          int32 value = utils::to_num(read_buffer+start_pos, end_pos-start_pos);
          if ((i % DEFAULT_B_SIZE) == 0) // Is the first record of Block?
          {
            b_title_bit_stream.PutBits(value-cur_field.min_value, cur_field.no_of_bits_per_value);
          }
          else if ((cur_field.is_delta_coding && !cur_field.block_desc[block_no].is_block_delta_constant) ||
            (!cur_field.is_delta_coding && !cur_field.block_desc[block_no].is_block_value_constant))
          {
            int32 to_store;
            if (cur_field.is_delta_coding)
            {
              to_store = value - prev_value[c_field] - cur_field.min_delta;
            }
            else
            {
              to_store = value - cur_field.min_value;
            }

            if (cur_field.Huffman_global)
            {
              const HuffmanEncoder::Code* codes = cur_field.Huffman_global->GetCodes();
              b_title_bit_stream.PutBits(codes[to_store].code, codes[to_store].len);
            }
            else
            {
              b_title_bit_stream.PutBits(to_store, cur_field.no_of_bits_per_num);
            }
          }

          prev_value[c_field] = value;
          continue;
        }

        if ((i % DEFAULT_B_SIZE) > 0 && cur_field.block_desc[block_no].is_block_constant)
        {
          continue;
        }

        if (!cur_field.is_len_constant)
        {
          b_title_bit_stream.PutBits(end_pos-start_pos - cur_field.min_len, cur_field.no_of_bits_per_len);
        }
        
        for (uint32 j = 0; j < end_pos-start_pos; ++j)
        {
          if (j >= cur_field.len || !cur_field.Ham_mask[j])
          {
            uchar c = read_buffer[start_pos+j];
            const HuffmanEncoder::Code* codes = cur_field.Huffman_local[MIN(j, MAX_FIELD_STAT_LEN)]->GetCodes();
            b_title_bit_stream.PutBits(codes[c].code, codes[c].len);
          }
        }
      }
    }
    b_title_bit_stream.FlushPartialWordBuffer();
  }
}

// --------------------------------------------------------------------------------------------
void StoreDNA(BitStream &dna_bit_stream, uchar *read_buffer, Record *&records, int64 no_records, 
  uint32 fastq_flags, std::vector<uchar> &symbols, uint32 *&sym_stats, uchar *sym_code)
{
  HuffmanEncoder::Code *sym_huf_codes = NULL;
  HuffmanEncoder *Huffman_sym = NULL;

  for (uint32 i = 0; i < symbols.size(); ++i)
  {
    dna_bit_stream.PutByte(symbols[i]);
  }
  dna_bit_stream.FlushPartialWordBuffer();

  if ((fastq_flags & FLAG_DNA_PLAIN) == 0)
  {
    sym_huf_codes = new HuffmanEncoder::Code[symbols.size()];
    Huffman_sym   = new HuffmanEncoder((uint32) symbols.size());

    for (uint32 i = 0; i < symbols.size(); ++i)
    {
      Huffman_sym->Insert(sym_stats[i]);
    }

    HuffmanEncoder::Code* codes = Huffman_sym->Complete();
    for (uint32 i = 0; i < symbols.size(); ++i)
    {
      sym_huf_codes[i] = codes[i];
    }

    HuffmanEncoder::StoreTree(dna_bit_stream, *Huffman_sym);
  }

  for (int64 i = 0; i < no_records; ++i)
  {
    uint32 seq_start = records[i].title_end_pos + 1; 
    uint32 seq_end   = records[i].seq_end_pos;

    for (uint32 j = seq_start; read_buffer[j] != '\0' && j < seq_end; ++j)
    {
      uchar c = read_buffer[j];
      if ((fastq_flags & FLAG_DNA_PLAIN) != 0) 
        dna_bit_stream.Put2Bits(sym_code[c]);
      else 
        dna_bit_stream.PutBits(sym_huf_codes[sym_code[c]].code, sym_huf_codes[sym_code[c]].len);
    }
  }

  dna_bit_stream.FlushPartialWordBuffer();

  // Release Huffman trees - DNA symbols
  if (Huffman_sym)
  {
    delete Huffman_sym;
    Huffman_sym = NULL;
    delete[] sym_huf_codes;
    sym_huf_codes = NULL;
  }
}

// --------------------------------------------------------------------------------------------
void StoreQuality(BitStream &quality_bit_stream, uchar *read_buffer, Record *&records, int64 no_records, 
  uint32 max_quality_length, uchar *qua_code, std::vector<uchar> &qualities, uint32 **&quality_stats, 
  HuffmanEncoder::Code **&qua_huf_codes, HuffmanEncoder::Code *&raw_qua_huf_codes)
{
  // Store qualities
  for (uint32 i = 0; i < qualities.size(); ++i)
  {
    quality_bit_stream.PutByte(qualities[i]);
  }
  quality_bit_stream.FlushPartialWordBuffer();

  qua_huf_codes = new HuffmanEncoder::Code*[max_quality_length+1];
  raw_qua_huf_codes = new HuffmanEncoder::Code[(max_quality_length+1)*qualities.size()];
  for (uint32 i = 0; i <= max_quality_length; ++i)
  {
    qua_huf_codes[i] = &raw_qua_huf_codes[i*qualities.size()];
  }

  for (uint32 i = 0; i <= max_quality_length; ++i)
  {
    HuffmanEncoder *huf = new HuffmanEncoder((uint32) qualities.size());

    for (uint32 j = 0; j < qualities.size(); ++j)
    {
      huf->Insert(quality_stats[i][j]);
    }

    HuffmanEncoder::Code* codes = huf->Complete();
    for (uint32 j = 0; j < qualities.size(); ++j)
    {
      qua_huf_codes[i][j] = codes[j];
    }
    HuffmanEncoder::StoreTree(quality_bit_stream, *huf);
  }

  quality_bit_stream.FlushPartialWordBuffer();
  
  for (int64 i = 0; i < no_records; ++i)
  {
    uint32 qua_start = records[i].seq_end_pos + 3; 
    uint32 qua_end = qua_start + (records[i].seq_end_pos - records[i].title_end_pos) - 1;

    for (uint32 j = qua_start, k = 0; j < qua_end && read_buffer[j] != '\0'; ++j, ++k)
    {
      int32 qua = qua_code[read_buffer[j]];
      quality_bit_stream.PutBits(qua_huf_codes[k+1][qua].code, qua_huf_codes[k+1][qua].len);
    }
  }

  quality_bit_stream.FlushPartialWordBuffer();
}

// --------------------------------------------------------------------------------------------
void MakeFooter(BitStream &footer_bit_stream, int32 g_size, uint32 FASTQ_size, int32 n_blocks, int32 n_subblocks, 
  int32 *all_wr_overlaps, std::map<double,int32> &blocks_order, uint32 *lb_sizes)  
{
  // Bits Encoding Processes Size
  int32 BEPS = floor(log2(g_size)) + 1;
  // Bits Encoding File Size
  int32 BEFS = floor(log2(FASTQ_size)) + 1;           
  // Bits Encoding Block Size
  int32 BEBS = floor(log2(n_blocks)) + 1;    
  // Bits Encoding SubBlock Size
  int32 BESS = floor(log2(n_subblocks)) + 1;     
  // Bits Encoding Last Blocks Size
  int32 BELB = floor(log2(*std::max_element(lb_sizes, lb_sizes+g_size))) + 1;
  // Bits Encoding Overlaps
  int32 BEOV = floor(log2(*std::max_element(all_wr_overlaps, all_wr_overlaps+g_size))) + 1;  
  // Are Last Blocks of Equal Size?
  const uint32 LBES = *std::max_element(lb_sizes, lb_sizes+g_size) != *std::min_element(lb_sizes, lb_sizes+g_size) ? 0 : 1; 
  // Create Footer buffer
  footer_bit_stream.Create(0);

  footer_bit_stream.PutBits(BEPS,4);
  footer_bit_stream.PutBits(BEFS,6);
  footer_bit_stream.PutBits(BEBS,4);
  footer_bit_stream.PutBits(BESS,4);
  footer_bit_stream.PutBits(BELB,5);
  footer_bit_stream.PutBits(BEOV,4);
  footer_bit_stream.PutBit(LBES);

  footer_bit_stream.PutBits(g_size, BEPS);
  footer_bit_stream.PutBits(FASTQ_size, BEFS);
  footer_bit_stream.PutBits(n_blocks, BEBS);
  footer_bit_stream.PutBits(n_subblocks, BESS);

  // Bitpack overlaps
  for (int32 i = 1; i < g_size; ++i)
  {
    footer_bit_stream.PutBits(all_wr_overlaps[i], BEOV);
  }

  // Recalculate to know the minimun amount of bits to encode the Processes ranks
  // from P0 to P(p-1)
  BEPS = ceil(log2(g_size));
  // Bitpack block order
  for (std::map<double, int>::iterator it = blocks_order.begin(); it != blocks_order.end(); ++it)
  {
      footer_bit_stream.PutBits(it->second, BEPS);
  }
  // Bitpack size of last blocks
  if (LBES == 0)
  {
    for (int32 i = 0; i < g_size; ++i)
    {
      footer_bit_stream.PutBits(lb_sizes[i], BELB);
    }
  }

  // Round up to the next complete byte (< than 8 bits wasted)
  footer_bit_stream.FlushPartialWordBuffer();

  // Add footer size at the end (2bytes)
  footer_bit_stream.Put2Bytes(footer_bit_stream.GetIO_Buffer_Pos());
}

// --------------------------------------------------------------------------------------------
void MakeHeader(BitStream &header_bit_stream, BlockHeader &p_block_header)  
{
  // Create header buffer
  header_bit_stream.Create(0);
  // Number of subblocks in this block
  int32 n_sb = p_block_header.sb_offset_list.size();

  header_bit_stream.PutBits(p_block_header.wr_id, p_block_header.BEWR);
  header_bit_stream.PutBits(n_sb, 6);
  header_bit_stream.PutBits(p_block_header.BESO, 5);
  header_bit_stream.PutBits(p_block_header.split_sb, 2);

  // Bitpack subblock's offsets
  for (int32 i = 0; i < n_sb; ++i)
  {
    header_bit_stream.PutBits(p_block_header.sb_offset_list.at(i), p_block_header.BESO);
  }

  // Round up to the next complete byte (< than 8 bits wasted)
  header_bit_stream.FlushPartialWordBuffer();
}
