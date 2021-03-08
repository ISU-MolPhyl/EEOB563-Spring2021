/*
    Copyright (C) 2015 Diego Darriba

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
#ifndef PLL_OPTIMIZE_H_
#define PLL_OPTIMIZE_H_

#ifndef PLL_H_
#define PLL_H_
#include "pll.h"
#endif

/* Parameters mask */
#define PLLMOD_OPT_PARAM_ALL                 (~0)
#define PLLMOD_OPT_PARAM_SUBST_RATES         (1<<0)
#define PLLMOD_OPT_PARAM_ALPHA               (1<<1)
#define PLLMOD_OPT_PARAM_PINV                (1<<2)
#define PLLMOD_OPT_PARAM_FREQUENCIES         (1<<3)
#define PLLMOD_OPT_PARAM_BRANCHES_SINGLE     (1<<4)
#define PLLMOD_OPT_PARAM_BRANCHES_ALL        (1<<5)
#define PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE  (1<<6)
#define PLLMOD_OPT_PARAM_TOPOLOGY            (1<<7)
#define PLLMOD_OPT_PARAM_FREE_RATES          (1<<8)
#define PLLMOD_OPT_PARAM_RATE_WEIGHTS        (1<<9)
#define PLLMOD_OPT_PARAM_BRANCH_LEN_SCALER   (1<<10)
/* !!! NOTE: all params in user code must be defined as
 *           PLLMOD_OPT_PARAM_USER<<0, PLLMOD_OPT_PARAM_USER<<1 etc. !!! */
#define PLLMOD_OPT_PARAM_USER                (1<<16)

/* L-BFGS-B bound type */
#define PLLMOD_OPT_LBFGSB_BOUND_NONE  0
#define PLLMOD_OPT_LBFGSB_BOUND_LOWER 1
#define PLLMOD_OPT_LBFGSB_BOUND_BOTH  2
#define PLLMOD_OPT_LBFGSB_BOUND_UPPER 3

/* Parameter defaults */
#define PLLMOD_OPT_DEFAULT_RATE_RATIO        1
#define PLLMOD_OPT_DEFAULT_FREQ_RATIO        1
#define PLLMOD_OPT_DEFAULT_PINV           0.01
#define PLLMOD_OPT_DEFAULT_ALPHA           0.5
#define PLLMOD_OPT_DEFAULT_BRANCH_LEN      0.1

/* Default parameter limits */
#define PLLMOD_OPT_MIN_BRANCH_LEN       1.0e-4
#define PLLMOD_OPT_MAX_BRANCH_LEN          100.
#define PLLMOD_OPT_TOL_BRANCH_LEN       1.0e-4
#define PLLMOD_OPT_MIN_SUBST_RATE       1.0e-3
#define PLLMOD_OPT_MAX_SUBST_RATE        1000.
#define PLLMOD_OPT_MIN_FREQ             1.0e-3
#define PLLMOD_OPT_MAX_FREQ               100.
#define PLLMOD_OPT_MIN_ALPHA            0.0201 //+ PLL_LBFGSB_ERROR
#define PLLMOD_OPT_MAX_ALPHA              100.
#define PLLMOD_OPT_MIN_PINV                  0
#define PLLMOD_OPT_MAX_PINV               0.99
#define PLLMOD_OPT_LNL_UNLIKELY         -1e+80

/* mixture models limits */
#define PLLMOD_OPT_MIN_RATE               0.02
#define PLLMOD_OPT_MAX_RATE               100.
#define PLLMOD_OPT_MIN_RATE_WEIGHT      1.0e-3
#define PLLMOD_OPT_MAX_RATE_WEIGHT        100.

/* branch length optimization methods */
#define PLLMOD_OPT_BLO_NEWTON_FAST        0 /* standard Newton-Raphson (NR) */
#define PLLMOD_OPT_BLO_NEWTON_SAFE        1 /* NR with per-branch LH check */
#define PLLMOD_OPT_BLO_NEWTON_FALLBACK    2 /* NR-FAST with fallback to NR-SAFE */
#define PLLMOD_OPT_BLO_NEWTON_GLOBAL      3 /* NR variant which looks for local optima */

#define PLLMOD_OPT_BLO_NEWTON_OLDFAST     10 /* old Newton-Raphson (a la IQTree) */
#define PLLMOD_OPT_BLO_NEWTON_OLDSAFE     11 /* old NR with per-branch LH check */

/* error codes (for this module, 2000-3000) */
#define PLLMOD_OPT_ERROR_PARAMETER           2000
#define PLLMOD_OPT_ERROR_TAXA_MISMATCH       2010
#define PLLMOD_OPT_ERROR_SEQLEN_MISMATCH     2020
#define PLLMOD_OPT_ERROR_ALIGN_UNREADABLE    2030
#define PLLMOD_OPT_ERROR_LBFGSB_UNKNOWN      2100
#define PLLMOD_OPT_ERROR_NEWTON_DERIV        2210
#define PLLMOD_OPT_ERROR_NEWTON_LIMIT        2220
#define PLLMOD_OPT_ERROR_NEWTON_UNKNOWN      2230
#define PLLMOD_OPT_ERROR_NEWTON_WORSE_LK     2240
#define PLLMOD_OPT_ERROR_NEWTON_BAD_RADIUS   2250
#define PLLMOD_OPT_ERROR_BRENT_INIT          2310

/* special options */
#define PLLMOD_OPT_BRLEN_OPTIMIZE_ALL  -1

/* Structure with information necessary for evaluating the likelihood */

/* Custom parameters structures provided by PLL for the
 * high level optimization functions (L-BFGS-B + Brent). */
typedef struct
{
  pll_partition_t * partition;
  pll_operation_t * operations;
  double * branch_lengths;
  unsigned int * matrix_indices;
  int rooted;
  const unsigned int * params_indices;
  union {
      struct {
        unsigned int root_clv_index;
        int scaler_index;
      } rooted_t;
      struct {
        unsigned int parent_clv_index;
        int parent_scaler_index;
        unsigned int child_clv_index;
        int child_scaler_index;
        unsigned int edge_pmatrix_index;
      } unrooted_t;
    } where;

  char __padding__[4];
  double alpha_value;
} pll_likelihood_info_t;

typedef struct
{
  pll_likelihood_info_t lk_params;
  unsigned int highest_freq_state;
  unsigned int highest_weight_state;
  //const unsigned int * params_indices;     /* indices according to rate cats */
  unsigned int params_index;         /* individual index to optimize */
  unsigned int which_parameters;
  int * subst_params_symmetries;
  double factr;
  double pgtol;

  double * sumtable;
} pll_optimize_options_t;

/* Custom parameters structure provided by PLL for the
 * high level optimization functions (Newton-Raphson). */
 typedef struct
{
  pll_partition_t * partition;
  pll_unode_t * tree;
  const unsigned int * params_indices;
  double * sumtable;
  double branch_length_min;
  double branch_length_max;
  double tolerance;
  int max_newton_iters;
  int opt_method;          /* see PLLMOD_OPT_BLO_* constants above */
} pll_newton_tree_params_t;

typedef struct
{
  pll_unode_t * tree;
  unsigned int partition_count;
  pll_partition_t ** partitions;
  unsigned int ** params_indices;
  double ** precomp_buffers;
  double ** brlen_buffers;
  double * brlen_orig;
  double * brlen_guess;
  int * converged;
  double * brlen_scalers;
  double branch_length_min;
  double branch_length_max;
  double tolerance;
  int max_newton_iters;
  int opt_method;          /* see PLLMOD_OPT_BLO_* constants above */
  int brlen_linkage;
  void * parallel_context;
  void (*parallel_reduce_cb)(void *,
                             double *,
                             size_t,
                             int);
} pll_newton_tree_params_multi_t;

/******************************************************************************/

/* functions in opt_algorithms.c */

/* core Newton-Raphson optimization function (multiple variables) */
PLL_EXPORT int pllmod_opt_minimize_newton_multi(unsigned int xnum,
                                                double xmin,
                                                double * xguess,
                                                double xmax,
                                                double tolerance,
                                                unsigned int max_iters,
                                                int * converged,
                                                void * params,
                                                void (deriv_func)(void *,
                                                          double *,
                                                          double *, double *));

/* core Newton-Raphson optimization function */
PLL_EXPORT double pllmod_opt_minimize_newton(double xmin,
                                             double xguess,
                                             double xmax,
                                             double tolerance,
                                             unsigned int max_iters,
                                             void * params,
                                             void (deriv_func)(void *,
                                                               double,
                                                               double *, double *));

/* core Newton-Raphson optimization function */
PLL_EXPORT double pllmod_opt_minimize_newton_old(double x1,
                                                 double xguess,
                                                 double x2,
                                                 double tolerance,
                                                 unsigned int max_iters,
                                                 void *params,
                                                 void (deriv_func)(void *,
                                                                   double,
                                                                   double *, double *));

/* core L-BFGS-B optimization function */
PLL_EXPORT double pllmod_opt_minimize_lbfgsb(double *x,
                                            double *xmin,
                                            double *xmax,
                                            int *bound,
                                            unsigned int n,
                                            double factr,
                                            double pgtol,
                                            void *params,
                                            double (*target_funk)(
                                                    void *,
                                                    double *));

/* core Brent optimization function */
PLL_EXPORT double pllmod_opt_minimize_brent(double xmin,
                                           double xguess,
                                           double xmax,
                                           double xtol,
                                           double *fx,
                                           double *f2x,
                                           void *params,
                                           double (*target_funk)(
                                               void *,
                                               double));

/* core Expectation-Maximization (EM) function */
PLL_EXPORT void pllmod_opt_minimize_em(double *w,
                                       unsigned int w_count,
                                       double *sitecat_lh,
                                       unsigned int *site_w,
                                       unsigned int l,
                                       void * params,
                                       double (*update_sitecatlk_funk)(
                                           void *,
                                           double *));

/******************************************************************************/

/* functions in pll_optimize.c */
PLL_EXPORT void pllmod_opt_derivative_func(void * parameters,
                                          double proposal,
                                          double *df, double *ddf);

/* high level optimization functions */
PLL_EXPORT double pllmod_opt_optimize_onedim(pll_optimize_options_t * p,
                                             double min,
                                             double max);

PLL_EXPORT double pllmod_opt_optimize_multidim(pll_optimize_options_t * p,
                                               double *umin,
                                               double *umax);

PLL_EXPORT double pllmod_opt_optimize_branch_lengths_iterative (
                                            pll_partition_t * partition,
                                            pll_unode_t * tree,
                                            const unsigned int * params_indices,
                                            double branch_length_min,
                                            double branch_length_max,
                                            double tolerance,
                                            int smoothings,
                                            int keep_update);

PLL_EXPORT double pllmod_opt_optimize_branch_lengths_local (
                                            pll_partition_t * partition,
                                            pll_unode_t * tree,
                                            const unsigned int * params_indices,
                                            double branch_length_min,
                                            double branch_length_max,
                                            double tolerance,
                                            int smoothings,
                                            int radius,
                                            int keep_update);


PLL_EXPORT double pllmod_opt_optimize_branch_lengths_local_multi (
                                              pll_partition_t ** partitions,
                                              size_t partition_count,
                                              pll_unode_t * tree,
                                              unsigned int ** params_indices,
                                              double ** sumtable_buffers,
                                              double ** brlen_buffers,
                                              double * brlen_scalers,
                                              double branch_length_min,
                                              double branch_length_max,
                                              double lh_epsilon,
                                              int max_iters,
                                              int radius,
                                              int keep_update,
                                              int opt_method,
                                              int brlen_linkage,
                                              void * parallel_context,
                                              void (*parallel_reduce_cb)(void *,
                                                                         double *,
                                                                         size_t,
                                                                         int));

PLL_EXPORT int pllmod_opt_minimize_brent_multi(unsigned int xnum,
                                               int * opt_mask,
                                               double * xmin,
                                               double * xguess,
                                               double * xmax,
                                               double xtol,
                                               double * xopt,
                                               double * fx,
                                               double * f2x,
                                               void * params,
                                               double (*target_funk)(
                                                       void *,
                                                      double *,
                                                      double *,
                                                      int *),
                                               int global_range);

PLL_EXPORT double pllmod_opt_minimize_lbfgsb_multi(unsigned int xnum,
                                                   double ** x,
                                                   double ** xmin,
                                                   double ** xmax,
                                                   int ** bound,
                                                   unsigned int * n,
                                                   unsigned int nmax,
                                                   double factr,
                                                   double pgtol,
                                                   void * params,
                                                   double (*target_funk)(void *,
                                                                         double **,
                                                                         double *,
                                                                         int *));



#endif /* PLL_OPTIMIZE_H_ */
