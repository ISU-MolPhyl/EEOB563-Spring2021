/*
 Copyright (C) 2016 Diego Darriba, Alexey Kozlov

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
#ifndef PLL_MSA_H_
#define PLL_MSA_H_

#ifndef PLL_H_
#define PLL_H_
#include "pll.h"
#endif

#define PLLMOD_MSA_STATS_NONE        (0ul)
#define PLLMOD_MSA_STATS_DUP_TAXA    (1ul<<0)
#define PLLMOD_MSA_STATS_DUP_SEQS    (1ul<<1)
#define PLLMOD_MSA_STATS_GAP_PROP    (1ul<<2)
#define PLLMOD_MSA_STATS_GAP_SEQS    (1ul<<3)
#define PLLMOD_MSA_STATS_GAP_COLS    (1ul<<4)
#define PLLMOD_MSA_STATS_INV_PROP    (1ul<<5)
#define PLLMOD_MSA_STATS_INV_COLS    (1ul<<6)
#define PLLMOD_MSA_STATS_FREQS       (1ul<<7)
#define PLLMOD_MSA_STATS_SUBST_RATES (1ul<<8)
#define PLLMOD_MSA_STATS_ALL         (~0ul)

#define PLLMOD_MSA_MAX_ERRORS        100


typedef struct msa_stats
{
  unsigned int states;

  unsigned long dup_taxa_pairs_count;
  unsigned long * dup_taxa_pairs;

  unsigned long dup_seqs_pairs_count;
  unsigned long * dup_seqs_pairs;

  double gap_prop;
  unsigned long gap_seqs_count;
  unsigned long * gap_seqs;
  unsigned long gap_cols_count;
  unsigned long * gap_cols;

  double inv_prop;
  unsigned long inv_cols_count;
  unsigned long * inv_cols;

  double * freqs;
  double * subst_rates;
} pllmod_msa_stats_t;

typedef struct msa_errors
{
  unsigned long invalid_char_count;
  char * invalid_chars;
  unsigned long * invalid_char_seq;
  unsigned long * invalid_char_pos;
  int status;
} pllmod_msa_errors_t;

PLL_EXPORT double * pllmod_msa_empirical_frequencies(pll_partition_t * partition);
PLL_EXPORT double * pllmod_msa_empirical_subst_rates(pll_partition_t * partition);
PLL_EXPORT double pllmod_msa_empirical_invariant_sites(pll_partition_t *partition);

PLL_EXPORT pllmod_msa_errors_t * pllmod_msa_check(const pll_msa_t * msa,
                                                  const pll_state_t * tipmap);

PLL_EXPORT void pllmod_msa_destroy_errors(pllmod_msa_errors_t * errs);

PLL_EXPORT pllmod_msa_stats_t * pllmod_msa_compute_stats(const pll_msa_t * msa,
                                                         unsigned int states,
                                                         const pll_state_t * tipmap,
                                                         const unsigned int * weights,
                                                         unsigned long stats_mask);

PLL_EXPORT void pllmod_msa_destroy_stats(pllmod_msa_stats_t * stats);

PLL_EXPORT pll_msa_t * pllmod_msa_filter(pll_msa_t * msa,
                                         unsigned long * remove_seqs,
                                         unsigned long remove_seqs_count,
                                         unsigned long * remove_cols,
                                         unsigned long remove_cols_count,
                                         unsigned int inplace);

PLL_EXPORT pll_msa_t ** pllmod_msa_split(const pll_msa_t * msa,
                                         const unsigned int * site_part,
                                         unsigned int part_count);

// This could be moved to a general file I/O module if we decide to have one
PLL_EXPORT int pllmod_msa_save_phylip(const pll_msa_t * msa,
                                      const char * out_fname);

#endif /* PLL_MSA_H_ */
