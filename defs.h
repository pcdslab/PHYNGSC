/*
  This file is part of phyNGSC (Hybrid MPI-OpenMP Strategy for Compression).  
  phyNGSC uses methods developed for DSRC version 1.00 (distributed under GNU GPL 2 licence)
  to underline de compression portion of the strategy.
  
  phyNGSC Authors: Sandino Vargas-Perez and Fahad Saeed
  DSRC Authors: Sebastian Deorowicz and Szymon Grabowski
*/

#ifndef _DEFS_H
#define _DEFS_H

#define BIT(x)							(1 << (x))
#define MIN(x,y)						((x) <= (y) ? (x) : (y))
#define MAX(x,y)						((x) >= (y) ? (x) : (y))
#define ABS(x)							((x) >=  0  ? (x) : -(x))
#define SIGN(x)							((x) >=  0  ?  1  : -1)
#define INVALID_BYTE					0xAA
#define INVALID_DWORD					0xAAAAAAAA
#define READ_BUFFER_SIZE				BIT(23)
#define WRITE_BUFFER_SIZE				BIT(23)

#if defined (_WIN32)
#	pragma warning(disable : 4996) // D_SCL_SECURE
#	pragma warning(disable : 4244) // conversion uint64 to uint32
#endif

typedef unsigned char uchar;
typedef int int32;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;


#if defined(DEBUG) || defined(_DEBUG)
#	include <assert.h>
#	define my_assert(x) assert(x)
#	define D_RESERVE_BYTES_PER_BLOCK		1
#	define D_RESERVE_BYTES_PER_SUPERBLOCK	1
#	define D_COMPUTE_RECORDS_CRC_PER_BLOCK	1
#else
#	define my_assert(x)
#	define D_RESERVE_BYTES_PER_BLOCK		0
#	define D_RESERVE_BYTES_PER_SUPERBLOCK	0
#	define D_COMPUTE_RECORDS_CRC_PER_BLOCK	0
#endif

#define COMPILE_TIME_ASSERT(COND, MSG) typedef char static_assertion_##MSG[(!!(COND))*2-1]
#define COMPILE_TIME_ASSERT1(X, L) COMPILE_TIME_ASSERT(X,static_assertion_at_line_##L)
#define COMPILE_TIME_ASSERT2(X, L) COMPILE_TIME_ASSERT1(X,L)
#define STATIC_ASSERT(X)    COMPILE_TIME_ASSERT2(X,__LINE__)


// --------------------------------------------------------------------------------------------
enum MemoryModeEnum
{
	MEM_MODE_NONE = 0, 
	MEM_MODE_READ, 
	MEM_MODE_WRITE, 
}; 

enum QualityEnum
{
	QUALITY_PLAIN = 0, 
	QUALITY_RLE, 
	QUALITY_PLAIN_TRUNC
};

enum FastqFlags
{
	FLAG_TRY_LZ					= BIT(0),
	FLAG_DNA_PLAIN				= BIT(1),
	FLAG_CONST_NUM_FIELDS		= BIT(2),
	FLAG_PLUS_ONLY				= BIT(3),
	FLAG_USE_DELTA				= BIT(4),
	FLAG_DELTA_CONSTANT			= BIT(5),
	FLAG_DELTA_NO_BEGIN_NUC		= BIT(6),
	FLAG_VARIABLE_LENGTH		= BIT(7),
	FLAG_LINE_BREAKS			= BIT(8),
};

enum SplittedSubBlock
{
	LSBS 						= BIT(0), // Last SubBlock splitted?
	FSBS 						= BIT(1)  // First SubBlock splitted?
};

// --------------------------------------------------------------------------------------------
struct Field;
struct Record;
struct ProcessCompressionInfo;
struct BlockHeader;

class BitMemory;
class BitStream;
class HuffmanEncoder;

#endif

