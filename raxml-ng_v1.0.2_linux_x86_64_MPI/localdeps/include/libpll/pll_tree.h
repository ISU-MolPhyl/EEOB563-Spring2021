/*
 Copyright (C) 2016 Diego Darriba

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
#ifndef PLL_TREE_H_
#define PLL_TREE_H_

#ifndef PLL_H_
#define PLL_H_
#include "pll.h"
#endif

#include "stdbool.h"

/**
 * PLL Tree utils module
 * Prefix: pll_tree_, pll_utree_, pll_rtree_
 */
#define PLLMOD_TREE_DEFAULT_BRANCH_LENGTH 0.1

/* error codes (for this module, 3000-4000) ; B = 2^10+2^11*/
/* TBR errors (B + {2^2,2^1,2^0}) */
#define PLLMOD_TREE_ERROR_TBR_LEAF_BISECTION   3073 // B + {001}
#define PLLMOD_TREE_ERROR_TBR_OVERLAPPED_NODES 3074 // B + {010}
#define PLLMOD_TREE_ERROR_TBR_SAME_SUBTREE     3075 // B + {011}
#define PLLMOD_TREE_ERROR_TBR_MASK             3079 // B + {111}

/* NNI errors (B + {2^4,2^3}) */
#define PLLMOD_TREE_ERROR_NNI_INVALID_MOVE     3080 // B + {01...}
#define PLLMOD_TREE_ERROR_NNI_MASK             3096 // B + {11...}

/* SPR errors (B + {2^6,2^5}) */
#define PLLMOD_TREE_ERROR_SPR_INVALID_NODE     3104 // B + {01...}
#define PLLMOD_TREE_ERROR_SPR_MASK             3168 // B + {11...}

/* general errors (B + {2^8,2^7}) */
#define PLLMOD_TREE_ERROR_INTERCHANGE_LEAF     3200 // B + {01...}
#define PLLMOD_TREE_ERROR_INVALID_REARRAGE     3328 // B + {10...}
#define PLLMOD_TREE_ERROR_INVALID_TREE_SIZE    3456 // B + {10...}
#define PLLMOD_TREE_ERROR_INVALID_TREE         3584 // B + {10...}
#define PLLMOD_TREE_ERROR_INVALID_SPLIT        3712 // B + {10...}
#define PLLMOD_TREE_ERROR_EMPTY_SPLIT          3840 // B + {10...}
#define PLLMOD_TREE_ERROR_INVALID_THRESHOLD    3968 // B + {10...}
#define PLLMOD_TREE_ERROR_POLYPHYL_OUTGROUP    3970 // B + {10...}

#define PLLMOD_TREE_REARRANGE_SPR  0
#define PLLMOD_TREE_REARRANGE_NNI  1
#define PLLMOD_TREE_REARRANGE_TBR  2

#define PLLMOD_TREEINFO_PARTITION_ALL -1

#define HASH_KEY_UNDEF ((unsigned int) -1)

typedef unsigned int pll_split_base_t;
typedef pll_split_base_t * pll_split_t;

typedef struct split_system_t
{
  pll_split_t * splits;
  double * support;
  unsigned int split_count;
  double max_support;
} pll_split_system_t;

typedef unsigned int hash_key_t;

typedef struct bitv_hash_entry
{
  hash_key_t key;
  pll_split_t bit_vector;
  unsigned int *tree_vector;
  unsigned int tip_count;
  double support;
  unsigned int bip_number;

  struct bitv_hash_entry *next;
} bitv_hash_entry_t;

typedef struct
{
  unsigned int table_size;
  bitv_hash_entry_t **table;
  unsigned int entry_count;
  unsigned int bit_count;     /* number of bits per entry */
  unsigned int bitv_len;      /* bitv length */
} bitv_hashtable_t;

typedef struct consensus_data_t
{
  pll_split_t split;
  unsigned int bit_count;
  double support;
} pll_consensus_data_t;


typedef struct consensus_utree_t
{
  pll_unode_t * tree;
  pll_consensus_data_t * branch_data;
  unsigned int tip_count;
  unsigned int branch_count;
} pll_consensus_utree_t;

typedef struct string_hash_entry
{
  hash_key_t key;
  int node_number;
  char * word;
  struct string_hash_entry *next;
} string_hash_entry_t;

typedef struct
{
  char **labels;
  unsigned int table_size;
  string_hash_entry_t **table;
  unsigned int entry_count;
} string_hashtable_t;

typedef struct pll_tree_edge
{
  union
  {
    struct
    {
      pll_unode_t * parent;
      pll_unode_t * child;
    } utree;
    struct
    {
      pll_rnode_t * parent;
      pll_rnode_t * child;
    } rtree;
  } edge;
    double length;
} pll_tree_edge_t;

typedef struct
{
  int rearrange_type;
  int rooted;
  double  likelihood;

  union {
    struct {
      void * prune_edge;
      void * regraft_edge;
      double prune_bl;        //! length of the pruned branch
      double prune_left_bl;   //! length of the removed branch when pruning
      double prune_right_bl;  //! length of the removed branch when pruning
      double regraft_bl;      //! length of the splitted branch when regrafting
    } SPR;
    struct {
      void * edge;
      double left_left_bl;
      double left_right_bl;
      double right_left_bl;
      double right_right_bl;
      double edge_bl;
      int type;
    } NNI;
    struct {
      void * bisect_edge;
      pll_tree_edge_t reconn_edge;
      double bisect_left_bl;
      double bisect_right_bl;
      double reconn_parent_left_bl;
      double reconn_parent_right_bl;
      double reconn_child_left_bl;
      double reconn_child_right_bl;
    } TBR;
  };
} pll_tree_rollback_t;

typedef struct treeinfo_edge
{
  unsigned int left_index;
  unsigned int right_index;
  unsigned int pmatrix_index;
  double brlen;
} pllmod_treeinfo_edge_t;

typedef struct treeinfo_topology
{
  unsigned int edge_count;
  unsigned int brlen_set_count;
  unsigned int root_index;
  pllmod_treeinfo_edge_t * edges;
  double ** branch_lengths;
} pllmod_treeinfo_topology_t;

typedef struct treeinfo
{
  // dimensions
  unsigned int tip_count;
  unsigned int partition_count;

  /* 0 = linked/shared, 1 = linked with scaler, 2 = unlinked */
  int brlen_linkage;
  double * linked_branch_lengths;

  pll_unode_t * root;
  pll_utree_t * tree;

  unsigned int subnode_count;
  pll_unode_t ** subnodes;

  // partitions & partition-specific stuff
  pll_partition_t ** partitions;
  double * alphas;
  int * gamma_mode; /* discrete GAMMA rates computation mode (mean, median) */
  unsigned int ** param_indices;
  int ** subst_matrix_symmetries;
  double ** branch_lengths;
  double * brlen_scalers;
  double * partition_loglh;
  int * params_to_optimize;

  // partition that have been initialized (useful for parallelization)
  unsigned int init_partition_count;
  unsigned int * init_partition_idx;
  pll_partition_t ** init_partitions;

  /* tree topology constraint */
  unsigned int * constraint;

  /* precomputation buffers for derivatives (aka "sumtable") */
  double ** deriv_precomp;

  // invalidation flags
  char ** clv_valid;
  char ** pmatrix_valid;

  // buffers
  pll_unode_t ** travbuffer;
  unsigned int * matrix_indices;
  pll_operation_t * operations;

  // partition on which all operations should be performed
  int active_partition;

  // general-purpose counter
  unsigned int counter;

  // parallelization stuff
  void * parallel_context;
  void (*parallel_reduce_cb)(void *, double *, size_t, int);
} pllmod_treeinfo_t;

typedef struct
{
  unsigned int node_count;
  pll_unode_t ** nodes;

  unsigned int partition_count;
  unsigned int * partition_indices;

  pll_utree_t * tree;
  double ** probs;
} pllmod_ancestral_t;

/* Topological rearrangements */
/* functions at pll_tree.c */

PLL_EXPORT int pllmod_utree_tbr(pll_unode_t * b_edge,
                                pll_tree_edge_t * r_edge,
                                pll_tree_rollback_t * rollback_info);

PLL_EXPORT int pllmod_utree_spr(pll_unode_t * p_edge,
                                pll_unode_t * r_edge,
                                pll_tree_rollback_t * rollback_info);

/* type = {PLL_NNI_NEXT, PLL_NNI_NEXTNEXT} */
PLL_EXPORT int pllmod_utree_nni(pll_unode_t * edge,
                                int type,
                                pll_tree_rollback_t * rollback_info);

PLL_EXPORT int pllmod_tree_rollback(pll_tree_rollback_t * rollback_info);

PLL_EXPORT pll_unode_t * pllmod_utree_serialize(pll_unode_t * tree,
                                                unsigned int tip_count);

PLL_EXPORT pll_utree_t * pllmod_utree_expand(pll_unode_t * serialized_tree,
                                             unsigned int tip_count);

PLL_EXPORT int pllmod_utree_draw_support(pll_utree_t * ref_tree,
                                         const double * support,
                                         pll_unode_t ** node_map,
                                         char * (*cb_serialize)(double));


/* Topological operations */

/* functions at rtree_operations.c */

PLL_EXPORT int pllmod_rtree_spr(pll_rnode_t * p_node,
                                pll_rnode_t * r_tree,
                                pll_rnode_t ** root,
                                pll_tree_rollback_t * rollback_info);

PLL_EXPORT int pllmod_rtree_get_sibling_pointers(pll_rnode_t * node,
                                                 pll_rnode_t ***self,
                                                 pll_rnode_t ***sister);

PLL_EXPORT pll_rnode_t * pllmod_rtree_prune(pll_rnode_t * node);

PLL_EXPORT int pllmod_rtree_regraft(pll_rnode_t * node,
                                    pll_rnode_t * tree);

/* functions at utree_operations.c */

PLL_EXPORT int pllmod_utree_bisect(pll_unode_t * edge,
                                   pll_unode_t ** parent_subtree,
                                   pll_unode_t ** child_subtree);

PLL_EXPORT pll_tree_edge_t pllmod_utree_reconnect(pll_tree_edge_t * edge,
                                                  pll_unode_t * pruned_edge);

PLL_EXPORT pll_unode_t * pllmod_utree_prune(pll_unode_t * edge);

PLL_EXPORT int pllmod_utree_regraft(pll_unode_t * edge,
                                    pll_unode_t * tree);

PLL_EXPORT int pllmod_utree_interchange(pll_unode_t * edge1,
                                        pll_unode_t * edge2);

PLL_EXPORT pll_unode_t * pllmod_utree_create_node(unsigned int clv_index,
                                                  int scaler_index,
                                                  char * label,
                                                  void * data);

PLL_EXPORT int pllmod_utree_connect_nodes(pll_unode_t * parent,
                                          pll_unode_t * child,
                                           double length);

/* Topological search */

/* functions at rtree_operations.c */

PLL_EXPORT int pllmod_rtree_nodes_at_node_dist(pll_rnode_t * root,
                                               pll_rnode_t ** outbuffer,
                                               unsigned int * node_count,
                                               int min_distance,
                                               int max_distance);

/* functions at utree_operations.c */

PLL_EXPORT int pllmod_utree_nodes_at_node_dist(pll_unode_t * node,
                                            pll_unode_t ** outbuffer,
                                            unsigned int * node_count,
                                            unsigned int min_distance,
                                            unsigned int max_distance);

PLL_EXPORT int pllmod_utree_nodes_at_edge_dist(pll_unode_t * edge,
                                               pll_unode_t ** outbuffer,
                                               unsigned int * node_count,
                                               unsigned int min_distance,
                                               unsigned int max_distance);

/* Tree construction */
/* functions at pll_tree.c */

PLL_EXPORT pll_utree_t * pllmod_utree_create_random(unsigned int taxa_count,
                                                    const char * const* names,
                                                    unsigned int random_seed);

PLL_EXPORT int pllmod_utree_extend_random(pll_utree_t * tree,
                                          unsigned int ext_taxa_count,
                                          const char * const* ext_names,
                                          unsigned int random_seed);

PLL_EXPORT
pll_utree_t * pllmod_utree_create_parsimony(unsigned int taxon_count,
                                            unsigned int seq_length,
                                            char * const * names,
                                            char * const * sequences,
                                            const unsigned int * site_weights,
                                            const pll_state_t * map,
                                            unsigned int states,
                                            unsigned int attributes,
                                            unsigned int random_seed,
                                            unsigned int * score);

pll_utree_t * pllmod_utree_create_parsimony_multipart(unsigned int taxon_count,
                                                      char * const * taxon_names,
                                                      unsigned int partition_count,
                                                      pll_partition_t * const * partitions,
                                                      unsigned int random_seed,
                                                      unsigned int * score);


/* Discrete operations */
/* functions at utree_distances.c */

PLL_EXPORT unsigned int pllmod_utree_rf_distance(pll_unode_t * t1,
                                                 pll_unode_t * t2,
                                                 unsigned int tip_count);

/* check that node ids and tip labels agree in both trees */
PLL_EXPORT int pllmod_utree_consistency_check(pll_utree_t * t1,
                                              pll_utree_t * t2);

/* if 2 different trees are parsed from newick node ids might have been set
   in a different order, so this function sets node ids in t2 such that
   node ids and tip labels agree in both trees */
PLL_EXPORT int pllmod_utree_consistency_set(pll_utree_t * t1,
                                            pll_utree_t * t2);

PLL_EXPORT unsigned int pllmod_utree_split_rf_distance(pll_split_t * s1,
                                                       pll_split_t * s2,
                                                       unsigned int tip_count);

PLL_EXPORT pll_split_t * pllmod_utree_split_create(const pll_unode_t * tree,
                                                   unsigned int tip_count,
                                                   pll_unode_t ** split_to_node_map);

PLL_EXPORT pll_split_t pllmod_utree_split_from_tips(unsigned int * subtree_tip_ids,
                                                    unsigned int subtree_size,
                                                    unsigned int tip_count);

PLL_EXPORT void pllmod_utree_split_normalize_and_sort(pll_split_t * s,
                                                      unsigned int tip_count,
                                                      unsigned int n_splits,
                                                      int keep_first);

PLL_EXPORT void pllmod_utree_split_show(pll_split_t split,
                                        unsigned int tip_count);

PLL_EXPORT void pllmod_utree_split_destroy(pll_split_t * split_list);

PLL_EXPORT unsigned int pllmod_utree_split_lightside(pll_split_t split,
                                                     unsigned int tip_count);

PLL_EXPORT unsigned int pllmod_utree_split_hamming_distance(pll_split_t s1,
                                                            pll_split_t s2,
                                                            unsigned int tip_count);

PLL_EXPORT pll_split_t * pll_utree_split_newick_string(char * s,
                                                       unsigned int tip_count,
                                                       string_hashtable_t * names_hash);

PLL_EXPORT
bitv_hashtable_t * pllmod_utree_split_hashtable_create(unsigned int tip_count,
                                                       unsigned int slot_count);

PLL_EXPORT bitv_hash_entry_t *
pllmod_utree_split_hashtable_insert_single(bitv_hashtable_t * splits_hash,
                                           pll_split_t split,
                                           double support);

PLL_EXPORT bitv_hashtable_t *
pllmod_utree_split_hashtable_insert(bitv_hashtable_t * splits_hash,
                                    pll_split_t * splits,
                                    unsigned int tip_count,
                                    unsigned int split_count,
                                    const double * support,
                                    int update_only);

PLL_EXPORT bitv_hash_entry_t *
pllmod_utree_split_hashtable_lookup(bitv_hashtable_t * splits_hash,
                                    pll_split_t split,
                                    unsigned int tip_count);

PLL_EXPORT
void pllmod_utree_split_hashtable_destroy(bitv_hashtable_t * hash);


/* functions in consensus.c */

PLL_EXPORT int pllmod_utree_compatible_splits(const pll_split_t s1,
                                              const pll_split_t s2,
                                              unsigned int split_len,
                                              unsigned int tip_count);

PLL_EXPORT pll_consensus_utree_t * pllmod_utree_from_splits(
                                        const pll_split_system_t * split_system,
                                        unsigned int tip_count,
                                        char * const * tip_labels);

PLL_EXPORT pll_split_system_t * pllmod_utree_split_consensus(
                                                bitv_hashtable_t * splits_hash,
                                                unsigned int tip_count,
                                                double threshold);

PLL_EXPORT pll_consensus_utree_t * pllmod_utree_weight_consensus(
                                                    pll_utree_t * const * trees,
                                                    const double * weights,
                                                    double threshold,
                                                    unsigned int tree_count);

PLL_EXPORT pll_consensus_utree_t * pllmod_utree_consensus(
                                                    const char * trees_filename,
                                                    double threshold,
                                                    unsigned int * tree_count);

PLL_EXPORT void pllmod_utree_split_system_destroy(pll_split_system_t * split_system);

PLL_EXPORT void pllmod_utree_consensus_destroy(pll_consensus_utree_t * tree);


/* Additional utilities */
/* functions at pll_tree.c */

PLL_EXPORT int pllmod_utree_set_clv_minimal(pll_unode_t * root,
                                         unsigned int tip_count);

PLL_EXPORT int pllmod_utree_traverse_apply(pll_unode_t * root,
                                        int (*cb_pre_trav)(pll_unode_t *,
                                                           void *),
                                        int (*cb_in_trav)(pll_unode_t *,
                                                          void *),
                                        int (*cb_post_trav)(pll_unode_t *,
                                                            void *),
                                        void *data);

PLL_EXPORT int pllmod_utree_is_tip(const pll_unode_t * node);

PLL_EXPORT void pllmod_utree_set_length(pll_unode_t * edge,
                                     double length);

PLL_EXPORT void pllmod_utree_set_length_recursive(pll_utree_t * tree,
                                                  double length,
                                                  int missing_only);


PLL_EXPORT void pllmod_utree_scale_branches(pll_utree_t * tree,
                                            double branch_length_scaler);

PLL_EXPORT void pllmod_utree_scale_branches_all(pll_unode_t * root,
                                                double branch_length_scaler);

PLL_EXPORT void pllmod_utree_scale_subtree_branches(pll_unode_t * root,
                                                    double branch_length_scaler);

PLL_EXPORT pll_utree_t * pllmod_utree_resolve_multi(const pll_utree_t * multi_tree,
                                                    unsigned int random_seed,
                                                    int * clv_index_map);

PLL_EXPORT int pllmod_utree_root_inplace(pll_utree_t * tree);

PLL_EXPORT int pllmod_utree_outgroup_root(pll_utree_t * tree,
                                          unsigned int * outgroup_tip_ids,
                                          unsigned int outgroup_size,
                                          int add_root_node);

PLL_EXPORT int pllmod_utree_collapse_branches(pll_utree_t * tree,
                                              double min_brlen);

PLL_EXPORT double pllmod_utree_compute_lk(pll_partition_t * partition,
                                       pll_unode_t * tree,
                                       const unsigned int * params_indices,
                                       int update_pmatrices,
                                       int update_partials);

PLL_EXPORT int pllmod_rtree_traverse_apply(pll_rnode_t * root,
                                           int (*cb_pre_trav)(pll_rnode_t *,
                                                              void *),
                                           int (*cb_in_trav)(pll_rnode_t *,
                                                             void *),
                                           int (*cb_post_trav)(pll_rnode_t *,
                                                               void *),
                                           void *data);

/* treeinfo */

PLL_EXPORT pllmod_treeinfo_t * pllmod_treeinfo_create(pll_unode_t * root,
                                                     unsigned int tips,
                                                     unsigned int partitions,
                                                     int brlen_linkage);

PLL_EXPORT
int pllmod_treeinfo_set_parallel_context(pllmod_treeinfo_t * treeinfo,
                                         void * parallel_context,
                                         void (*parallel_reduce_cb)(void *,
                                                                    double *,
                                                                    size_t,
                                                                    int op));

PLL_EXPORT int pllmod_treeinfo_init_partition(pllmod_treeinfo_t * treeinfo,
                                           unsigned int partition_index,
                                           pll_partition_t * partition,
                                           int params_to_optimize,
                                           int gamma_mode,
                                           double alpha,
                                           const unsigned int * param_indices,
                                           const int * subst_matrix_symmetries);

PLL_EXPORT int pllmod_treeinfo_set_active_partition(pllmod_treeinfo_t * treeinfo,
                                                    int partition_index);

PLL_EXPORT int pllmod_treeinfo_set_root(pllmod_treeinfo_t * treeinfo,
                                        pll_unode_t * root);

PLL_EXPORT
int pllmod_treeinfo_get_branch_length_all(const pllmod_treeinfo_t * treeinfo,
                                          const pll_unode_t * edge,
                                          double * lengths);

PLL_EXPORT int pllmod_treeinfo_set_branch_length(pllmod_treeinfo_t * treeinfo,
                                                 pll_unode_t * edge,
                                                 double length);

PLL_EXPORT
int pllmod_treeinfo_set_branch_length_all(pllmod_treeinfo_t * treeinfo,
                                          pll_unode_t * edge,
                                          const double * lengths);

PLL_EXPORT
int pllmod_treeinfo_set_branch_length_partition(pllmod_treeinfo_t * treeinfo,
                                                pll_unode_t * edge,
                                                int partition_index,
                                                double length);

PLL_EXPORT
pll_utree_t * pllmod_treeinfo_get_partition_tree(const pllmod_treeinfo_t * treeinfo,
                                                 int partition_index);

PLL_EXPORT
pllmod_treeinfo_topology_t * pllmod_treeinfo_get_topology(const pllmod_treeinfo_t * treeinfo,
                                                          pllmod_treeinfo_topology_t * topol);

PLL_EXPORT
int pllmod_treeinfo_set_topology(pllmod_treeinfo_t * treeinfo,
                                 const pllmod_treeinfo_topology_t * topol);

PLL_EXPORT
int pllmod_treeinfo_destroy_topology(pllmod_treeinfo_topology_t * topol);

PLL_EXPORT int pllmod_treeinfo_destroy_partition(pllmod_treeinfo_t * treeinfo,
                                                 unsigned int partition_index);

PLL_EXPORT void pllmod_treeinfo_destroy(pllmod_treeinfo_t * treeinfo);

PLL_EXPORT int pllmod_treeinfo_update_prob_matrices(pllmod_treeinfo_t * treeinfo,
                                                    int update_all);

PLL_EXPORT void pllmod_treeinfo_invalidate_all(pllmod_treeinfo_t * treeinfo);

PLL_EXPORT int pllmod_treeinfo_validate_clvs(pllmod_treeinfo_t * treeinfo,
                                             pll_unode_t ** travbuffer,
                                             unsigned int travbuffer_size);

PLL_EXPORT void pllmod_treeinfo_invalidate_pmatrix(pllmod_treeinfo_t * treeinfo,
                                                   const pll_unode_t * edge);

PLL_EXPORT void pllmod_treeinfo_invalidate_clv(pllmod_treeinfo_t * treeinfo,
                                               const pll_unode_t * edge);

PLL_EXPORT double pllmod_treeinfo_compute_loglh(pllmod_treeinfo_t * treeinfo,
                                                int incremental);

PLL_EXPORT double pllmod_treeinfo_compute_loglh_flex(pllmod_treeinfo_t * treeinfo,
                                                     int incremental,
                                                     int update_pmatrices);

PLL_EXPORT double pllmod_treeinfo_compute_loglh_persite(pllmod_treeinfo_t * treeinfo,
                                                        int incremental,
                                                        double ** persite_lnl);

PLL_EXPORT
int pllmod_treeinfo_scale_branches_all(pllmod_treeinfo_t * treeinfo, double scaler);

PLL_EXPORT
int pllmod_treeinfo_scale_branches_partition(pllmod_treeinfo_t * treeinfo,
                                             unsigned int partition_idx,
                                             double scaler);

PLL_EXPORT
int pllmod_treeinfo_normalize_brlen_scalers(pllmod_treeinfo_t * treeinfo);

PLL_EXPORT int pllmod_treeinfo_set_tree(pllmod_treeinfo_t * treeinfo,
                                        pll_utree_t * tree);

PLL_EXPORT int pllmod_treeinfo_set_constraint_clvmap(pllmod_treeinfo_t * treeinfo,
                                                     const int * clv_index_map);

PLL_EXPORT int pllmod_treeinfo_set_constraint_tree(pllmod_treeinfo_t * treeinfo,
                                                   const pll_utree_t * cons_tree);

PLL_EXPORT int pllmod_treeinfo_check_constraint(pllmod_treeinfo_t * treeinfo,
                                                pll_unode_t * subtree,
                                                pll_unode_t * regraft_edge);

PLL_EXPORT pllmod_ancestral_t * pllmod_treeinfo_compute_ancestral(pllmod_treeinfo_t * treeinfo);

PLL_EXPORT void pllmod_treeinfo_destroy_ancestral(pllmod_ancestral_t * ancestral);

/* tbe_functions.c */

typedef struct refsplit_info
{
  unsigned int p;
  bool subtree_res;
  unsigned int left_leaf_idx;
  unsigned int right_leaf_idx;
} pllmod_tbe_split_info_t;

PLL_EXPORT
pllmod_tbe_split_info_t * pllmod_utree_tbe_nature_init(pll_unode_t * ref_root,
                                                       unsigned int tip_count,
                                                       const pll_unode_t ** split_to_node_map);

/* Compute Transfer Support (Lemoine et al., Nature 2018) for every split in ref_splits.
 * Sarahs implementation of the algorithm from the Nature paper. */
PLL_EXPORT int pllmod_utree_tbe_nature(pll_split_t * ref_splits,
                                       pll_split_t * bs_splits,
                                       pll_unode_t* bs_root,
                                       unsigned int tip_count,
                                       double * support,
                                       pllmod_tbe_split_info_t* split_info);

/* This is an old, naive and rather inefficient TBE computation method by Alexey.
 * Keep it here just in case */
PLL_EXPORT int pllmod_utree_tbe_naive(pll_split_t * ref_splits,
                                      pll_split_t * bs_splits,
                                      unsigned int tip_count,
                                      double * support);

#endif /* PLL_TREE_H_ */
