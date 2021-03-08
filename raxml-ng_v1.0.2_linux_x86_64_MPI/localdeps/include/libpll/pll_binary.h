/*
 Copyright (C) 2016 Diego Darriba, Pierre Barbera

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Contact: Diego Darriba <Diego.Darriba@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */
#ifndef PLLMOD_BIN_H_
#define PLLMOD_BIN_H_

#ifndef PLL_H_
#define PLL_H_
#include "pll.h"
#endif

#define PLLMOD_BIN_BLOCK_PARTITION  0
#define PLLMOD_BIN_BLOCK_CLV        1
#define PLLMOD_BIN_BLOCK_TREE       2
#define PLLMOD_BIN_BLOCK_CUSTOM     3
#define PLLMOD_BIN_BLOCK_REPEATS    4

#define PLLMOD_BIN_ACCESS_SEQUENTIAL  0
#define PLLMOD_BIN_ACCESS_RANDOM      1
#define PLLMOD_BIN_ACCESS_SEEK       -1

#define PLLMOD_BIN_INVALID_OFFSET    -1

#define PLLMOD_BIN_ATTRIB_UPDATE_MAP              (1<<0)
#define PLLMOD_BIN_ATTRIB_PARTITION_DUMP_CLV      (1<<1)
#define PLLMOD_BIN_ATTRIB_PARTITION_DUMP_WGT      (1<<2)
#define PLLMOD_BIN_ATTRIB_ALIGNED                 (1<<3)
#define PLLMOD_BIN_ATTRIB_PARTITION_LOAD_SKELETON (1<<4)

#define PLLMOD_BIN_ERROR_BLOCK_MISMATCH         4001
#define PLLMOD_BIN_ERROR_BLOCK_LENGTH           4002
#define PLLMOD_BIN_ERROR_BINARY_IO              4003
#define PLLMOD_BIN_ERROR_INVALID_INDEX          4010
#define PLLMOD_BIN_ERROR_INVALID_SIZE           4011
#define PLLMOD_BIN_ERROR_LOADSTORE              4012
#define PLLMOD_BIN_ERROR_MISSING_BLOCK          4020

/*
 * This is the main header of the binary stream.
 * Access type can be sequential or random. If sequential, right after the
 * header comes the first block, and we can jump from one block to the next one
 * using `block_len`. If access is random, after the header comes a table with
 * the offsets of the different blocks, such can we access directly.
 */
typedef struct
{
  unsigned int n_blocks;      //! number of blocks in the file
  unsigned int max_blocks;    //! maximum number of blocks (size of block map)
  unsigned int access_type;   //! PLLMOD_BIN_ACCESS_{SEQUENTIAL|RANDOM}
  char pad[4];                //! padding
  long map_offset;            //! offset of the block map
} pll_binary_header_t;

/* block map for random access */
typedef struct
{
  long block_id;             //! user-defined block id
  long block_offset;         //! offset in the file
} pll_block_map_t;

/*
 * Header stored before each block
 * If the binary file was created for random access, it may be important that
 * attributes contain PLLMOD_BIN_ATTRIB_UPDATE_MAP such that the block map
 * is updated. Otherwise the block will be only accessible sequentially after
 * reading the previous block.
 */
typedef struct
{
  long block_id;             //! user-defined block id
  unsigned int type;         //! block type PLLMOD_BIN_BLOCK_...
  unsigned int attributes;   //! custom block attributes
  unsigned int alignment;    //! if memory should be aligned
  char pad[4];               //! padding
  size_t block_len;          //! block length
} pll_block_header_t;

PLL_EXPORT FILE * pllmod_binary_create(const char * filename,
                                       pll_binary_header_t * header,
                                       unsigned int access_type,
                                       unsigned int n_blocks);

PLL_EXPORT FILE * pllmod_binary_open(const char * filename,
                                     pll_binary_header_t * header);

PLL_EXPORT FILE * pllmod_binary_append_open(const char * filename,
                                            pll_binary_header_t * header);

PLL_EXPORT int pllmod_binary_close(FILE * bin_file);

PLL_EXPORT pll_block_map_t * pllmod_binary_get_map(FILE * bin_file,
                                                   unsigned int * n_blocks);

PLL_EXPORT int pllmod_binary_partition_dump(FILE * bin_file,
                                         int block_id,
                                         pll_partition_t * partition,
                                         unsigned int attributes);

PLL_EXPORT pll_partition_t * pllmod_binary_partition_load(FILE * bin_file,
                                                    int block_id,
                                                    pll_partition_t * partition,
                                                    unsigned int * attributes,
                                                    long int offset);

PLL_EXPORT int pllmod_binary_repeats_dump(FILE * bin_file,
                                      int block_id,
                                      pll_partition_t * partition,
                                      unsigned int attributes);

PLL_EXPORT int pllmod_binary_repeats_load(FILE * bin_file,
                                      int block_id,
                                      pll_partition_t * partition,
                                      unsigned int * attributes,
                                      long int offset);

PLL_EXPORT int pllmod_binary_pernoderepeats_dump(FILE * bin_file,
                                      int block_id,
                                      pll_partition_t * partition,
                                      unsigned int clv_index,
                                      unsigned int attributes);

PLL_EXPORT int pllmod_binary_pernoderepeats_load(FILE * bin_file,
                                      int block_id,
                                      pll_partition_t * partition,
                                      unsigned int clv_index,
                                      unsigned int * attributes,
                                      long int offset);

PLL_EXPORT int pllmod_binary_clv_dump(FILE * bin_file,
                                      int block_id,
                                      pll_partition_t * partition,
                                      unsigned int clv_index,
                                      unsigned int attributes);

PLL_EXPORT int pllmod_binary_clv_load(FILE * bin_file,
                                      int block_id,
                                      pll_partition_t * partition,
                                      unsigned int clv_index,
                                      unsigned int * attributes,
                                      long int offset);

PLL_EXPORT int pllmod_binary_utree_dump(FILE * bin_file,
                                        int block_id,
                                        pll_unode_t * tree,
                                        unsigned int tip_count,
                                        unsigned int attributes);

PLL_EXPORT pll_unode_t * pllmod_binary_utree_load(FILE * bin_file,
                                                  int block_id,
                                                  unsigned int * attributes,
                                                  long int offset);

PLL_EXPORT int pllmod_binary_custom_dump(FILE * bin_file,
                                        int block_id,
                                        void * data,
                                        size_t size,
                                        unsigned int attributes);

PLL_EXPORT void * pllmod_binary_custom_load(FILE * bin_file,
                                           int block_id,
                                           size_t * size,
                                           unsigned int * type,
                                           unsigned int * attributes,
                                           long int offset);

#endif /* PLLMOD_BIN_H_ */
