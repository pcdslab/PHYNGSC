/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
*/

#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <map>
#include <math.h>
#include <algorithm>
#include "defs.h"
#include "utils.h"
#include "structures.h"
#include "bit_stream.h"
#include "huffman.h"
#include "tasks.h"

// --------------------------------------------------------------------------------------------
// main()
// --------------------------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  int32 g_size, p_rank, provided;
  MPI_Status status;
  int32 p_subblock_count = 0;

  MPI_File input_FASTQ, output_NGSC;
  int32 err_read, err_write;

  double p_timer_start, p_timer_end;
  
  MPI_Offset FASTQ_size;

  MPI_Offset p_working_region, r_buffer_size = READ_BUFFER_SIZE, w_buffer_size = WRITE_BUFFER_SIZE;
  MPI_Offset p_bytes_read = 0, p_bytes_written = 0, p_bytes_to_copy = 0;
  MPI_Offset p_wr_start, p_wr_end, r_buffer_curr_pos = 0;
  uchar *read_buffer, *write_buffer, *copy_buffer;

  int32 overlap = 500, wr_ov_used = 0;
  int64 no_records = 0, rec_start_pos = 0;
  int32 no_threads = 0; 
  int64 records_per_th = 100000;

  ProcessCompressionInfo p_compress_info, *all_info = NULL;
  BlockHeader p_block_header;
  std::vector<double> timestamps;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &g_size);

  if (argc != 4) {
    if (p_rank == 0){
      fprintf(stderr, "\n[E] ERROR: Incorrect number of arguments. Usage:\n");
      fprintf(stderr, "\n           mpiexec -np p ./phyNGSC input_filename.fastq output_filename.ngsc num_threads.\n");
      fprintf(stderr, "                 mpiexec:              Command to run MPI applications.\n");
      fprintf(stderr, "                 -np p:                Number of MPI processes to be used, where p is a number greater than 1.\n");
      fprintf(stderr, "                 ./phyNGSC:            phyNGSC application.\n");
      fprintf(stderr, "                 input_filename.fastq: Name of FASTQ file.\n");
      fprintf(stderr, "                 output_filename.ngsc: Name of NGSC file (created).\n");
      fprintf(stderr, "                 num_threads:          Number of threads to be used per MPI process (must be greater than 1).\n");
    }
    MPI_Finalize();
    exit(1);
  }

  // Open input/output files
  err_read  = MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &input_FASTQ);
  err_write = MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &output_NGSC);
  no_threads = atoi(argv[3]);

  if (no_threads)
    records_per_th /= no_threads;
  else
  {
    if (p_rank == 0)
      fprintf(stderr, "\n[E] ERROR: Number of threads is less than 1.\n");
    MPI_Finalize();
    exit(1);
  }

  if (g_size < 2)
  {
    if (p_rank == 0)
      fprintf(stderr, "\n[E] ERROR: Number of MPI Processes (-np) is less than 2.\n");
    MPI_Finalize();
    exit(1);
  }

  if (err_read || err_write) {

    if (p_rank == 0)
      fprintf(stderr, "\n[E] ERROR: MPI failed to open file <<%s>>.[%d][%d]\n", err_read ? argv[1] : argv[2], err_read, err_write);
    MPI_Finalize();
    exit(2);
  }

  if(p_rank == 0)
    printf("\n[I] INFO: Processing <<%s>> with %d MPI processes and %d threads per process.\n", argv[1], g_size, atoi(argv[3]));

  // Start timer
  p_timer_start = MPI_Wtime();

  MPI_File_get_size(input_FASTQ, &FASTQ_size);
  p_working_region = FASTQ_size/g_size;
  p_wr_start = p_rank * p_working_region;
  p_wr_end = (p_rank != g_size - 1) ? (p_wr_start + p_working_region + overlap - 1) : FASTQ_size - 1;

  // If the p_working_region is smaller than 8Mb, then read that amount
  if (p_working_region < r_buffer_size)
  {
    r_buffer_size = p_wr_end - p_wr_start + 1;
    if (p_rank == g_size - 1)
      overlap = 0;
  }

  write_buffer = (uchar*) malloc(w_buffer_size * sizeof(uchar));
  read_buffer  = (uchar*) malloc(r_buffer_size * sizeof(uchar));
  MPI_File_read_at(input_FASTQ, p_wr_start, read_buffer, r_buffer_size, MPI_CHAR, MPI_STATUS_IGNORE);

  // Each p_rank (except for 0) locate the start of the first complete record in their working region
  int32 counter = 0, firstAtSign = 0;
  if (p_rank != 0) 
  {
    char checkChar = '@';

    while(true){
      if (read_buffer[counter] == checkChar)
      {
        if(checkChar == '@')
        {
          firstAtSign = counter;
          checkChar = '\n';
        }
        if (read_buffer[counter] == '\n')
        {
          if (read_buffer[counter+1] == '@')
            r_buffer_curr_pos = counter + 1;
          else
            r_buffer_curr_pos = firstAtSign; 
          
          break;
        }
      }
      counter++;        
    }
  }

  rec_start_pos = r_buffer_curr_pos;
  // Get the actual overlap used to start a complete rec, to build footer
  wr_ov_used = r_buffer_curr_pos;
  // Add the working region ID to the headers of this processor
  p_block_header.wr_id = p_rank;
  p_block_header.BEWR = ceil(log2(g_size));
  p_block_header.split_sb = 0;


  // Begin processing FASTQ file with each p_rank's corresponding working region
  // --------------------------------------------------------------------------------------------
  while(p_bytes_read < p_working_region)
  {
    std::vector<Field> fields;
    uint32 n_fields = 0;
    std::vector<uint32> dna_occ(256,0);
    uint32 fastq_flags;
    std::vector<uchar> symbols;
    std::vector<uchar> qualities;
    uchar sym_code[256]={0};
    uchar qua_code[256]={0};
    uint32 *sym_stats = NULL;
    uint32 **quality_stats = NULL;
    uint32 *raw_quality_stats = NULL;
    HuffmanEncoder::Code **qua_huf_codes = NULL;
    HuffmanEncoder::Code *raw_qua_huf_codes = NULL;

    uchar trans_amb_codes[256];
    bool alphabet_qua[256] = {false};

    std::fill_n(trans_amb_codes, 256, 0);
    trans_amb_codes['Y'] = 2;
    trans_amb_codes['R'] = 3;
    trans_amb_codes['W'] = 4;
    trans_amb_codes['S'] = 5;
    trans_amb_codes['K'] = 6;
    trans_amb_codes['M'] = 7;
    trans_amb_codes['D'] = 8;
    trans_amb_codes['V'] = 9;
    trans_amb_codes['H'] = 10;
    trans_amb_codes['B'] = 11;
    trans_amb_codes['N'] = 12;
    trans_amb_codes['X'] = 13;
    trans_amb_codes['U'] = 14;
    trans_amb_codes['.'] = 15;
    trans_amb_codes['-'] = 16;
    trans_amb_codes['A'] = 1;
    trans_amb_codes['C'] = 1;
    trans_amb_codes['T'] = 1;
    trans_amb_codes['G'] = 1;

    const char *c_separators = " ._,=:/-#\n";
    const std::vector<uchar> separators(c_separators, c_separators+strlen(c_separators));

    int64 *tmp_title_pos = new int64[no_threads*records_per_th];
    int64 *tmp_seq_pos   = new int64[no_threads*records_per_th];
    int64 *th_rec_count = new int64[no_threads];

    ++p_subblock_count;

    if(!tmp_title_pos || !tmp_seq_pos)
    {
      printf("\n[E] ERROR: p_Rank %d can't allocate memory for tmp_title_pos or tmp_seq_pos.\n", p_rank);
      printf("           Calling MPI_Finalize() and exit(3). 3 = mem alloc returned NULL)\n");
      MPI_Finalize();
      exit(3);
    }

    BitStream p_write_buff_bit_stream;
    BitStream h_title_bit_stream, b_title_bit_stream;
    BitStream dna_bit_stream;
    BitStream quality_bit_stream;

    // Variables to analize quality and dna seq
    uint32 max_quality_length = 0;
    uint32 min_quality_length = (uint32) -1;
    uint32 max_sequence_length = 0;
    uint32 min_sequence_length = (uint32) -1;
    uint32 global_max_sequence_length = 0;

    if(p_bytes_read != 0)
    {
      read_buffer = (uchar*) malloc(r_buffer_size * sizeof(uchar));

      if(!read_buffer)
      {
        printf("\n[E] ERROR: p_Rank %d can't allocate memory for the read_buffer.\n", p_rank);
        printf("           Calling MPI_Finalize() and exit(3). 3 = mem alloc returned NULL)\n");
        MPI_Finalize();
        exit(3);
      }

      MPI_File_read_at(input_FASTQ, r_buffer_curr_pos, read_buffer, r_buffer_size, MPI_CHAR, MPI_STATUS_IGNORE);
    }

    // OpenMP Parallel Region: Threads locate the start of each record in read_buffer
    // --------------------------------------------------------------------------------------------
    #pragma omp parallel num_threads(no_threads) reduction(+ : no_records)
    {
      // Obtain and print thread id
      int32 th_id = omp_get_thread_num();
      int64 i, th_i_start , th_i_end;
      bool check = true;
      int64 th_rec_no = th_id * records_per_th;

      th_i_start = th_id * (r_buffer_size) / no_threads;
      th_i_end = (th_id + 1) * (r_buffer_size) / no_threads;

      if(th_id == 0)
        th_i_start = rec_start_pos;

      for (i = th_i_start; i < th_i_end; ++i)
      {
        // Check for the 1st title of thread th_id
        if (check)
        {
          int64 th_1st_title_ends[5];
          int32 pos = 0;

          while(check)
          {
            if (read_buffer[i] == '\n')
            {
              th_1st_title_ends[pos++] = i;
              if (read_buffer[i + 1] == '+' && read_buffer[i + 2] == '\n')
              {
                if (pos > 1)
                {
                  tmp_title_pos[th_rec_no] = th_1st_title_ends[pos - 2];
                  break;
                }
                i += 2;
                th_1st_title_ends[pos++] = i;
              }
            }
            i++;
          }

          tmp_seq_pos[th_rec_no] = i;
          i += 2 + (i -  tmp_title_pos[th_rec_no]) + 1;

          ++th_rec_no;
          check = false;  
        }
        

        if (read_buffer[i] != '\n')
          continue;

        tmp_title_pos[th_rec_no] = i;
        
        while(read_buffer[++i] != '\n')
          continue;

        tmp_seq_pos[th_rec_no] = i;
        i +=  (tmp_seq_pos[th_rec_no] - tmp_title_pos[th_rec_no]) + 3;
        ++th_rec_no;

        if (th_id == (no_threads - 1) && i >= (r_buffer_size - overlap)){
          i = th_i_end ;
          continue;
        }

        // Check that threads keep using their memory space
        if (th_rec_no > ((th_id * records_per_th) + records_per_th))
        {
          printf("\n[!] WARNING: Thread %d in p_Rank %d can't allocate more records.\n", th_id, p_rank);
          printf("             Thread will exit loop with only %lld records\n", th_rec_no);
          i = th_i_end ;
        }
      }
      th_rec_count[th_id] = th_rec_no - (th_id * records_per_th);
      // Total number of records in the team of threads
      no_records  += th_rec_count[th_id];
    }

    Record *records = new Record[no_records];
    if(!records)
    {
      printf("\n[E] ERROR: p_Rank %d can't allocate memory for %lld records.\n", p_rank, no_records);
      printf("           Calling MPI_Finalize() and exit(3). 3 = mem alloc returned NULL)\n");
      MPI_Finalize();
      exit(3);
    }

    int64 f_start = rec_start_pos;

    // Initializing vector Fields to analize record titles
    for (uint32 i = rec_start_pos; i <= tmp_title_pos[0]; ++i)
    {
      if (!count(separators.begin(), separators.end(), read_buffer[i]))
        continue;

      fields.push_back(Field());

      fields[n_fields].data            = new uchar[i-f_start+1];
      std::copy(read_buffer+f_start, read_buffer+i, fields[n_fields].data);
      fields[n_fields].data[i-f_start] = '\0';
      fields[n_fields].len             = i - f_start;
      fields[n_fields].max_len         = fields[n_fields].len;
      fields[n_fields].min_len         = fields[n_fields].len;
      fields[n_fields].sep             = read_buffer[i];
      fields[n_fields].is_constant     = true;
      fields[n_fields].is_len_constant = true;
      fields[n_fields].is_numeric      = fields[n_fields].IsNum();
      fields[n_fields].Ham_mask        = new bool[fields[n_fields].len];

      if (fields[n_fields].is_numeric)
      {
        fields[n_fields].min_value = fields[n_fields].ToNum();
        fields[n_fields].max_value = fields[n_fields].min_value;
        fields[n_fields].num_values[fields[n_fields].min_value]++;
      }

      for (uint32 k = 0; k < fields[n_fields].len; ++k)
      {
        fields[n_fields].Ham_mask[k] = true;
      }
      fields[n_fields].block_desc.clear();

      f_start = i+1;
      n_fields++;
    }

    // OpenMP Parallel Region: Threads populate the Records structure
    // --------------------------------------------------------------------------------------------
    #pragma omp parallel num_threads(no_threads)
    {
      int32 th_id       = omp_get_thread_num();
      int64 th_cur_rec  = 0;
      int64 th_f_start  = rec_start_pos;
      int64 cur_pos     = th_id * records_per_th;

      for (int32 r = th_id; r > 0; --r)
        th_cur_rec += th_rec_count[r - 1];

      for (int64 r = th_cur_rec; r < (th_cur_rec + th_rec_count[th_id]); ++r)
      {
        records[r].title_end_pos  = tmp_title_pos[cur_pos];
        records[r].seq_end_pos    = tmp_seq_pos[cur_pos];
        ++cur_pos;

        // If other than the first record of the first thread
        if (r || th_id)
        {
          th_f_start = records[r].title_end_pos;
          while(read_buffer[th_f_start - 1] != '\n')
            th_f_start--;
        }

        my_assert(read_buffer[th_f_start] != '@');

        for (int64 f = th_f_start; f <= records[r].title_end_pos; ++f)
        {
          if (!count(separators.begin(), separators.end(), read_buffer[f]))
            continue;

          records[r].fields_end_pos.push_back(f);
        }

        if(records[r].fields_end_pos.size() != n_fields)
        {
          printf("\n[!] WARNING: Thread %d in p_Rank %d has diff number of fields.\n", th_id, p_rank);
          printf("             Thread has %lu/%u fields.\n", records[r].fields_end_pos.size(), n_fields);
        }
      }
    }

    // Release temporary variables memory
    delete[] tmp_title_pos;
    delete[] tmp_seq_pos;
    delete[] th_rec_count;
    tmp_title_pos = NULL;
    tmp_seq_pos   = NULL;
    th_rec_count  = NULL;

    fastq_flags = FLAG_PLUS_ONLY | FLAG_DNA_PLAIN | FLAG_CONST_NUM_FIELDS;
    fastq_flags &= ~(FLAG_LINE_BREAKS | FLAG_VARIABLE_LENGTH | FLAG_TRY_LZ | FLAG_USE_DELTA | 
    FLAG_DELTA_CONSTANT | FLAG_DELTA_NO_BEGIN_NUC);

    // OpenMP Parallel Region: Threads analize, process and compress records in buffer.
    // --------------------------------------------------------------------------------------------
    #pragma omp parallel shared(records, read_buffer, fields, rec_start_pos, dna_occ) num_threads(no_threads)
    {
      std::vector<uint32> th_dna_occ(256, 0);
      bool th_alphabet_qua[256] = {false};
      uint32 local_max_qua_len  = 0, local_max_seq_len = 0;
      uint32 local_min_qua_len  = (uint32) -1, local_min_seq_len = (uint32) -1;

      bool is_delta             = false;
      bool is_delta_constant    = true;
      bool has_no_begin_nuc     = false;
      bool is_first_th_rec      = true;

      uchar sequence_start      = 0;
      uchar quality_start       = 0;

      // Analize fields data in title portion of records
      #pragma omp for nowait
      for (uint32 i = 0; i < n_fields; ++i)
      {
        AnalyzeTitleFields(read_buffer, rec_start_pos, i, fields, records, no_records);
      }
      
      // Analize ambiguity in DNA and transfer to Quality
      #pragma omp for nowait
      for (int64 i = 0; i < no_records; ++i)
      { 
        uint32 j, k, q, tmp;
        uint32 dna_seq_length = records[i].seq_end_pos - records[i].title_end_pos;
        uint32 seq_start = records[i].title_end_pos + 1;
        uint32 qua_start = records[i].seq_end_pos + 3; 
        uint32 seq_end = records[i].seq_end_pos;
        uint32 qua_end = qua_start + dna_seq_length - 1;
        uint32 qua_len = dna_seq_length, seq_len = 0;
        
        if (is_first_th_rec)
        {
          sequence_start = read_buffer[seq_start];
          quality_start  = read_buffer[qua_start];

          is_delta = read_buffer[seq_start+1] >= '0' && read_buffer[seq_start+1] <= '3';

          if (sequence_start <= '3' && sequence_start >= '0')
          {
            has_no_begin_nuc = true;
            is_delta_constant = false;
          }

          is_first_th_rec = false;
        }
        
        if (is_delta && !has_no_begin_nuc && is_delta_constant)
        {
          is_delta_constant &= read_buffer[seq_start] == sequence_start;
          is_delta_constant &= read_buffer[qua_start] == quality_start;  
        }

        if(is_delta)
        {
          static const char delta_A[] = {'N', 'N', 'A', 'C', 'G', 'T'};
          static const char delta_C[] = {'N', 'N', 'C', 'A', 'T', 'G'};
          static const char delta_G[] = {'N', 'N', 'G', 'T', 'A', 'C'};
          static const char delta_T[] = {'N', 'N', 'T', 'G', 'C', 'A'};

          const uint32 translation = (uint32)is_delta_constant;
          const char* last_matrix = delta_A;

          uchar symbol;
          if (!has_no_begin_nuc)
          {
            symbol = read_buffer[seq_start]; 
            th_dna_occ[symbol]++;
          }
          else
          {
            symbol = 'A';
            my_assert(is_delta_constant == false);
          }
          
          uint32 n = 1 - (uint32)has_no_begin_nuc;
          for ( ; n < dna_seq_length; ++n)
          {
            switch (symbol)
            {
              case 'A': last_matrix = delta_A; break;
              case 'C': last_matrix = delta_C; break;
              case 'G': last_matrix = delta_G; break;
              case 'T': last_matrix = delta_T; break;
              case 'N': 
              default: break;
            }
            symbol = last_matrix[read_buffer[seq_start+n]-'.'];
            my_assert(symbol =='A' || symbol == 'C' || symbol == 'G' || symbol == 'T' || symbol == 'N');
            th_dna_occ[symbol]++;

            read_buffer[seq_start+n-translation] = symbol;
            read_buffer[qua_start+n-translation] = read_buffer[seq_start+n];
          }
          
          for ( ; n < qua_len; ++n)
          {
            read_buffer[qua_start+n-translation] = read_buffer[seq_start+n];
          }

          read_buffer[seq_end-translation] = 0;
          read_buffer[qua_end-translation] = 0;
          seq_end -= translation;
          qua_len -= translation;
          qua_end -= translation;
        }

        // Check whether make transfer
        bool make_transfer = false;
        bool possible_transfer = true;
        for (j = seq_start, q = qua_start; j < seq_end; ++j, ++q)
        {
          if (trans_amb_codes[read_buffer[j]] == 1)
            continue;

          if (trans_amb_codes[read_buffer[j]] == 0)
          {
            possible_transfer = false;
            break;
          }
          else
          {
            if (read_buffer[q] < 33 || read_buffer[q] > 40)
            {
              possible_transfer = false;
              break;
            }
            make_transfer = true;
          }
        }

        if (make_transfer && possible_transfer)
        {
          for (j = k = seq_start; j < seq_end; ++j)
          {
            if ((tmp = trans_amb_codes[read_buffer[j]]) > 1)
            {
              read_buffer[j + dna_seq_length + 2] = (uchar)(128 + (tmp << 3) - 16 + (read_buffer[j + dna_seq_length + 2] - 33));
            }
            else
            {
              read_buffer[k++] = read_buffer[j];
            }
          }
          read_buffer[k] = '\0';
          seq_end = k;
        }

        seq_len = seq_end - seq_start;

        if (!is_delta)
        {
          for (k = seq_start; k < seq_end; ++k)
          {
            ++th_dna_occ[read_buffer[k]];
          }
        }

        if (qua_len > local_max_qua_len)
          local_max_qua_len = qua_len;

        if (qua_len < local_min_qua_len)
          local_min_qua_len = qua_len;

        if (seq_len > local_max_seq_len)
          local_max_seq_len = seq_len;

        if (seq_len < local_min_seq_len)
          local_min_seq_len = seq_len;

        // Identify quality alphabet
        for (j = qua_start; j < qua_end; ++j)
          th_alphabet_qua[read_buffer[j]] = true;
      }

      // Reduce each thread record analisis result
      #pragma omp critical
      { 
        for (uint32 k = 0; k < 256; ++k)
        {
          dna_occ[k] += th_dna_occ[k];
          alphabet_qua[k] |= th_alphabet_qua[k];
        }

        if (local_max_qua_len > max_quality_length)
          max_quality_length = local_max_qua_len;

        if (local_min_qua_len > min_quality_length)
          min_quality_length = local_min_qua_len;

        if (local_max_seq_len > max_sequence_length)
          max_sequence_length = local_max_seq_len;

        if (local_min_seq_len > min_sequence_length)
          min_sequence_length = local_min_seq_len;

        if (max_sequence_length > global_max_sequence_length)
          global_max_sequence_length = max_sequence_length;

        if (is_delta)
          fastq_flags |= FLAG_USE_DELTA;

        if (is_delta_constant)
          fastq_flags |= FLAG_DELTA_CONSTANT;

        if (has_no_begin_nuc)
          fastq_flags |= FLAG_DELTA_NO_BEGIN_NUC;
      }

      // Only one thread prepares all variables and flags to start compression
      #pragma omp single 
      {
        bool is_length_variable = min_quality_length != max_quality_length;
        is_length_variable &= min_sequence_length != max_sequence_length;

        if (is_length_variable)
        {
          fastq_flags |= FLAG_VARIABLE_LENGTH;
        }

        symbols.clear();
        qualities.clear();

        int32 sym_idx = 0;
        int32 qua_idx = 0;
        for (uint32 i = 0; i < 256; ++i)
        {
          if (dna_occ[i])
          {
            symbols.push_back((uchar) i);
            sym_code[i] = (uchar) sym_idx++;
          }
          if (alphabet_qua[i])
          {
            qualities.push_back((uchar) i);
            qua_code[i] = (uchar) qua_idx++;
          }
        }        
      }

      // Analize, compress and store record information
      #pragma omp sections
      {

        // Store title header information
        #pragma omp section
        {
          h_title_bit_stream.Create(1);
          StoreTitleHeader(h_title_bit_stream, fields);
        }

        // Store title information
        #pragma omp section
        {
          b_title_bit_stream.Create(1);
          StoreTitleBody(b_title_bit_stream, fields, read_buffer, records, no_records, rec_start_pos);
        }

        // Analize DNA data and store information
        #pragma omp section
        {
          AnalyzeDNA(fastq_flags, dna_occ, symbols, sym_stats);
          dna_bit_stream.Create(1);
          StoreDNA(dna_bit_stream, read_buffer, records, no_records, fastq_flags, symbols, sym_stats, sym_code);
          dna_occ.clear();
        }

        // Analyze quality data and store information
        #pragma omp section
        {
          AnalyzeQuality(read_buffer, records, no_records, max_quality_length, qualities, qua_code, quality_stats, raw_quality_stats);
          quality_bit_stream.Create(1);
          StoreQuality(quality_bit_stream, read_buffer, records, no_records, max_quality_length, qua_code, qualities, quality_stats, qua_huf_codes, raw_qua_huf_codes);
        }

        // Store meta information about the current data to be written to output file
        #pragma omp section
        {
          p_write_buff_bit_stream.Create(1);
          p_write_buff_bit_stream.PutWord(no_records);
          p_write_buff_bit_stream.PutWord(max_quality_length);
          p_write_buff_bit_stream.PutWord(global_max_sequence_length);
          p_write_buff_bit_stream.PutByte((uchar) symbols.size());
          p_write_buff_bit_stream.PutByte((uchar) QUALITY_PLAIN);
          p_write_buff_bit_stream.PutByte((uchar) qualities.size());
        }
      }
    }

    p_write_buff_bit_stream.PutWord(fastq_flags);
    p_write_buff_bit_stream.FlushPartialWordBuffer();

    uint32 quality_len_bits = BitStream::BitLength(max_quality_length);
    if ((fastq_flags & FLAG_VARIABLE_LENGTH) != 0)
    {
      for (uint32 i = 0; i < no_records; ++i)
      {
        uint32 qua_len = records[i].seq_end_pos - records[i].title_end_pos;
        p_write_buff_bit_stream.PutBits(qua_len, quality_len_bits);
      }
      p_write_buff_bit_stream.FlushPartialWordBuffer();
    }

    // Prepare for next read from the FASTQ working region
    p_bytes_read += (records[no_records - 1].seq_end_pos * 2) - records[no_records - 1].title_end_pos + 3;
    r_buffer_curr_pos = p_bytes_read + p_wr_start;

    if ((r_buffer_curr_pos + r_buffer_size) > p_wr_end)
    {
      if(p_rank == g_size - 1)
        overlap = 0;
      
      // Chunk to read will be g_size = rest of bytes in working region
      r_buffer_size = p_wr_end - r_buffer_curr_pos;
    }

    // Release read_buffer allocated memory
    free(read_buffer);

    // Release records allocated memory
    delete[] records;
    records = NULL;
    no_records = 0;
    rec_start_pos = 0;

    // Release fields allocated memory
    fields.clear();
    n_fields = 0;

    // Release DNA and quality stats allocated memory
    if (sym_stats)
    {
      delete[] sym_stats;
      sym_stats = NULL;
    }
    if (quality_stats)
    {
      delete[] quality_stats;
      delete[] raw_quality_stats;
      quality_stats = NULL;
      raw_quality_stats = NULL;
    }

    // Release quality Huffman codes allocated memory
    if (qua_huf_codes)
    {
      delete[] qua_huf_codes;
      delete[] raw_qua_huf_codes;
      qua_huf_codes = NULL;
      raw_qua_huf_codes = NULL;
    }

    uint32 w_title_h_len = h_title_bit_stream.GetIO_Buffer_Pos();
    uint32 w_title_b_len = b_title_bit_stream.GetIO_Buffer_Pos();
    uint32 w_dna_len = dna_bit_stream.GetIO_Buffer_Pos();
    uint32 w_quality_len = quality_bit_stream.GetIO_Buffer_Pos();
    uint32 w_info_buff_len = p_write_buff_bit_stream.GetIO_Buffer_Pos();
    
    // Allocate an specific ammount of memory for the copy buffer
    p_bytes_to_copy = w_title_h_len + w_title_b_len + w_dna_len + w_quality_len + w_info_buff_len;
    copy_buffer = (uchar*) malloc((p_bytes_to_copy) * sizeof(uchar));

    // Copy SubBlock's data into copy buffer
    // --------------------------------------------------------------------------------------------
    #pragma omp parallel sections shared(copy_buffer) num_threads(no_threads)
    {
      // Copy general information about the compressed chunk to copy buffer
      #pragma omp section
      {
        uchar *tem_buffer = p_write_buff_bit_stream.GetIO_Buffer();
        std::copy(tem_buffer, tem_buffer+w_info_buff_len, copy_buffer);
        p_write_buff_bit_stream.Close();
      }

      // Copy title header information to copy buffer
      #pragma omp section
      {
        uchar *tem_buffer = h_title_bit_stream.GetIO_Buffer();
        uint32 bytes_in_w_buff = w_info_buff_len;
        std::copy(tem_buffer, tem_buffer+w_title_h_len, copy_buffer+bytes_in_w_buff);
        h_title_bit_stream.Close();
      }

      // Copy title information to copy buffer
      #pragma omp section
      {
        uchar *tem_buffer = b_title_bit_stream.GetIO_Buffer();
        uint32 bytes_in_w_buff = w_info_buff_len + w_title_h_len;
        std::copy(tem_buffer, tem_buffer+w_title_b_len, copy_buffer+bytes_in_w_buff);
        b_title_bit_stream.Close();
      }

      // Copy DNA information to copy buffer
      #pragma omp section
      {
        uchar *tem_buffer = dna_bit_stream.GetIO_Buffer();
        uint32 bytes_in_w_buff = w_info_buff_len + w_title_h_len + w_title_b_len;
        std::copy(tem_buffer, tem_buffer+w_dna_len, copy_buffer+bytes_in_w_buff);
        dna_bit_stream.Close();
      }

      // Copy quality information
      #pragma omp section
      {
        uchar *tem_buffer = quality_bit_stream.GetIO_Buffer();
        uint32 bytes_in_w_buff = w_info_buff_len + w_title_h_len + w_title_b_len + w_dna_len;
        std::copy(tem_buffer, tem_buffer+w_quality_len, copy_buffer+bytes_in_w_buff);
        quality_bit_stream.Close();
      }
    }

    // Add information to the header for the current block
    p_block_header.sb_offset_list.push_back(p_bytes_to_copy);
    p_block_header.BESO = floor(log2(*std::max_element(p_block_header.sb_offset_list.begin(), p_block_header.sb_offset_list.end()))) + 1;
    
    int32 h_size = p_block_header.HeaderSize();

    // If the write buffer will be full with current subblock, 
    // then build block with splitted subblock and wirte it to NGSC file
    // else copy subblock to write buffer
    // --------------------------------------------------------------------------------------------
    if ((p_bytes_written+p_bytes_to_copy+h_size) > w_buffer_size)
    {
      BitStream header_bit_stream;
      p_block_header.split_sb |= LSBS;

      // Get the number of bytes needed to fill the write_buffer.
      uint32 bytes_fill_buffer = w_buffer_size - (p_bytes_written + h_size);
      // Copy bytes_fill_buffer bytes to write_buffer, to fill it up.
      std::copy(copy_buffer, copy_buffer+bytes_fill_buffer, write_buffer+p_bytes_written);
      // Make space to insert header
      std::copy_backward(write_buffer, write_buffer+(p_bytes_written+bytes_fill_buffer), write_buffer+w_buffer_size);

      // Build header
      MakeHeader(header_bit_stream, p_block_header);
      uchar *header_buffer = header_bit_stream.GetIO_Buffer();
      // Add header to Block
      std::copy(header_buffer, header_buffer+h_size, write_buffer);

      header_bit_stream.Close();

      // Write complete block of compress data to NGSC
      MPI_File_write_shared(output_NGSC, write_buffer, w_buffer_size, MPI_CHAR, &status);
      // Take the time of completion of the writing operation
      timestamps.push_back(MPI_Wtime());
      
      free(write_buffer);
      write_buffer = (uchar*) malloc((w_buffer_size) * sizeof(uchar));

      if(!write_buffer)
      {
        printf("[E] ERROR: p_Rank %d can't allocate memory for the write_buffer.\n", p_rank);
        printf("           Calling MPI_Finalize() and exit(3). 3 = mem alloc returned NULL)\n");
        MPI_Finalize();
        exit(3);
      }

      // Copy rest of the bytes in copy_buffer to write_buffer.
      std::copy(copy_buffer+bytes_fill_buffer, copy_buffer+p_bytes_to_copy, write_buffer);
      p_bytes_written = p_bytes_to_copy - bytes_fill_buffer;
      p_block_header.split_sb |= FSBS;
      p_block_header.split_sb &= ~LSBS;
      p_block_header.sb_offset_list.clear();
    }
    else
    {
      // Copy subblock in copy_buffer to write_buffer
      std::copy(copy_buffer, copy_buffer+p_bytes_to_copy, write_buffer+p_bytes_written);
      p_bytes_written += p_bytes_to_copy;

    }

    free(copy_buffer);
  } 

  // If the write_buffer is not empty
  // --------------------------------------------------------------------------------------------
  if(p_bytes_written > 0)
  {
    BitStream header_bit_stream;
    MakeHeader(header_bit_stream, p_block_header);
    uchar *header_buffer = header_bit_stream.GetIO_Buffer();
    int32 header_size = header_bit_stream.GetIO_Buffer_Pos();

    std::copy_backward(write_buffer, write_buffer+p_bytes_written, write_buffer+(p_bytes_written+header_size));
    std::copy(header_buffer, header_buffer+header_size, write_buffer);

    p_bytes_written += header_size;

    header_bit_stream.Close();

    MPI_File_write_shared(output_NGSC, write_buffer, p_bytes_written, MPI_CHAR, &status);
    timestamps.push_back(MPI_Wtime());
    free(write_buffer);
  }

  // Gather information to build footer
  // --------------------------------------------------------------------------------------------
  if (p_rank == 0) 
  {
    all_info = (ProcessCompressionInfo*) malloc(g_size * sizeof(ProcessCompressionInfo));
  } 
  p_compress_info.n_blocks = timestamps.size();
  p_compress_info.n_subblocks = p_subblock_count;
  p_compress_info.last_block_size = p_bytes_written;
  p_compress_info.wr_overlap = wr_ov_used;

  // create an MPI type for struct ProcessCompressionInfo
  const int32 nitems = 4;
  int32 blocklengths[4] = {1,1,1,1};
  MPI_Datatype types[4] = {MPI_INT, MPI_INT, MPI_UNSIGNED, MPI_INT};
  MPI_Datatype mpi_p_info;
  MPI_Aint offsets[4];
  MPI_Aint addr[5];

  MPI_Get_address(&p_compress_info, &addr[0]);
  MPI_Get_address(&p_compress_info.n_blocks, &addr[1]);
  MPI_Get_address(&p_compress_info.n_subblocks, &addr[2]);
  MPI_Get_address(&p_compress_info.last_block_size, &addr[3]);
  MPI_Get_address(&p_compress_info.wr_overlap, &addr[4]);

  offsets[0] = addr[1] - addr[0];
  offsets[1] = addr[2] - addr[0];
  offsets[2] = addr[3] - addr[0];
  offsets[3] = addr[4] - addr[0];

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_p_info);
  MPI_Type_commit(&mpi_p_info);

  // Gather information from all processes
  MPI_Gather(&p_compress_info, 1, mpi_p_info, &(all_info[0]), 1, mpi_p_info, 0, MPI_COMM_WORLD);

  int32  n_total_blocks     = 0;    // Total number of Blocks in NGSC.
  int32  n_total_subblocks  = 0;    // Total number of SubBlocks in NGSC.
  int32  *displs            = NULL; // Displacements.
  int32  *timestamps_counts = NULL; // Number of timestamps in the toc list of each p_rank.
  double *TOCW_buffer       = NULL; // Buffer to save all timestamps gathered from each p_rank.
  uint32 *lb_sizes          = NULL; // Size of each p_rank's last Block written.
  int32  *all_wr_overlaps   = NULL; // Overlap used at the beginin of each p_rank's WR.

  if (p_rank == 0)
  {
    displs            = (int32*) malloc(g_size * sizeof(int32));
    timestamps_counts = (int32*) malloc(g_size * sizeof(int32));
    all_wr_overlaps   = (int32*) malloc(g_size * sizeof(int32));
    lb_sizes          = (uint32*) malloc(g_size * sizeof(uint32));

    for (int32 i = 0; i < g_size; ++i)
    {
      displs[i] = n_total_blocks;
      n_total_blocks += all_info[i].n_blocks;
      n_total_subblocks += all_info[i].n_subblocks;
      timestamps_counts[i] = all_info[i].n_blocks;
      all_wr_overlaps[i] = all_info[i].wr_overlap;
      lb_sizes[i]	= all_info[i].last_block_size;
    }

    TOCW_buffer = (double*) malloc(n_total_blocks * sizeof(double));
  }
  // Gather time of completion lists from all processes
  MPI_Gatherv(timestamps.data(), timestamps.size(), MPI_DOUBLE, TOCW_buffer, timestamps_counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Sort list of time of completion to stablish order of writing
  if (p_rank == 0)
  {
    std::map<double,int32> blocks_order;
    for (int32 i = 0; i < g_size; ++i)
    {
      for (int32 j = 0; j < timestamps_counts[i]; ++j)
      {
        int32 ts_pos = j + displs[i];
        blocks_order[TOCW_buffer[ts_pos]] = i;
      }
    }
    
    BitStream footer_bit_stream;

    // Build footer
    MakeFooter(footer_bit_stream, g_size, FASTQ_size, n_total_blocks, n_total_subblocks, 
                all_wr_overlaps, blocks_order, lb_sizes);

    uchar *footer_buffer = footer_bit_stream.GetIO_Buffer();
    int32 footer_size = footer_bit_stream.GetIO_Buffer_Pos();
    
    // Write footer to NGSC file
    MPI_File_write_shared(output_NGSC, footer_buffer, footer_size, MPI_CHAR, &status);
    // Close stream
    footer_bit_stream.Close();

    // Release memory
    free(all_info);
    free(displs);
    free(timestamps_counts);
    free(TOCW_buffer);
    free(lb_sizes);
    free(all_wr_overlaps);  
  }
 
  // Stop timer
  p_timer_end = MPI_Wtime();

  if (p_rank == 0)
    printf("\nRANK\tCOMP_TIME\tN_BLOCK\tN_SUBBLOCKS\n----------------------------------------------\n");

  MPI_Barrier(MPI_COMM_WORLD);
  printf("%d\t%f\t%ld\t%d\n", p_rank, p_timer_end-p_timer_start, timestamps.size(), p_subblock_count);

  // Wait for all ranks, before closing the I/O files.
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_File_close(&input_FASTQ);
  MPI_File_close(&output_NGSC);

  MPI_Finalize();
  return 0;
}