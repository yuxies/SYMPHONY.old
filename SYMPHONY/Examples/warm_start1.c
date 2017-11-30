/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* The author of this file is Menal Guzelsoy                                 */
/*                                                                           */
/* (c) Copyright 2005-2015 Lehigh University. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifdef USE_OSI_INTERFACE

#include "OsiSymSolverInterface.hpp"

int main(int argc, char **argv)
{

  OsiSymSolverInterface si;
  si.parseCommandLine(argc, argv);
  si.loadProblem();
  si.setSymParam(OsiSymKeepWarmStart, true);
  si.setSymParam(OsiSymFindFirstFeasible, true);
   /* set node selection rule to DEPTH_FIRST_SEARCH */
  si.setSymParam(OsiSymSearchStrategy, 3);

  si.initialSolve();

  si.setSymParam(OsiSymFindFirstFeasible, false);
   /* set node selection rule to BEST_FIRST_SEARCH */
  si.setSymParam(OsiSymSearchStrategy, 4);

  si.resolve();

  return(0);

}

#else

#include "symphony.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>

int main(int argc, char **argv)
{
   sym_environment *env = sym_open_environment();
   sym_parse_command_line(env, argc, argv);
   sym_load_problem(env);

   sym_set_int_param(env, "keep_warm_start", TRUE);
   sym_set_int_param(env, "warm_start_type", 1);
   sym_set_int_param(env, "sensitivity_analysis", TRUE);
   sym_set_int_param(env, "sensitivity_rhs", TRUE);
   sym_set_int_param(env, "sensitivity_bounds", TRUE);
   sym_set_int_param(env, "do_reduced_cost_fixing", 0);
   sym_set_int_param(env, "do_primal_heuristic", 0);
   sym_set_int_param(env, "generate_cgl_cuts", 0);
   sym_set_int_param(env, "should_use_rel_br", 0);
   sym_set_int_param(env, "prep_level", 0);
   sym_set_int_param(env, "node_limit", 500);

   /* Write out the problem at hand */
   char *file_name = (char *) malloc(CSIZE * 80);
   sym_get_str_param(env, "problem_name", &file_name);

   char *file_name1 = (char *) malloc(CSIZE * 80);
   sprintf(file_name1, "%s_base", file_name);
   sym_write_mps(env, file_name1);

   sym_solve(env);

   /* Read instance again and do further processing for WS from scratch */
   sym_environment *env2 = sym_open_environment();
   sym_parse_command_line(env2, argc, argv);
   sym_load_problem(env2);

   /* Perturbation BEGIN */
   /* Read parameters specific to WS: percent_decrease, percent_increase, srand_seed */
   double percent_increase = 0.0;
   double percent_decrease = 0.0;
   unsigned int srand_seed = 1;
   int i;
   for (i = 0; i < argc; i++) {
      if (strcmp("--up_change", argv[i]) == 0)
         sscanf(argv[i+1], "%lf", &percent_increase);
      else if (strcmp("--down_change", argv[i]) == 0)
         sscanf(argv[i+1], "%lf", &percent_decrease);
      else if (strcmp("--srand_seed", argv[i]) == 0)
         sscanf(argv[i+1], "%u", &srand_seed);
   }

   // Assumption: no new rows or columns added into new instance
   int numcols, numrows;
   sym_get_num_rows(env2, &numrows);
   sym_get_num_cols(env2, &numcols);

   /* Generate a sequence of random binary digits to decide which rhs to
    * increase (1) and which demand to decrese (0) */
   int *rand_seq = (int *) malloc(ISIZE * numrows);
   srand(srand_seed);
   for (i = 0; i < numrows; i++) {
      rand_seq[i] = rand()%2;
   }

   int num_rhs_to_change = numrows, counter = 0;
   double temp;

   char *newsense = (char *) malloc(CSIZE * numrows);
   sym_get_row_sense(env2, newsense);
   double *newrange = (double *) malloc(DSIZE * numrows);
   sym_get_row_range(env2, newrange);
   double *rhs = (double *) malloc(DSIZE * numrows);
   sym_get_rhs(env2, rhs);
   //FIXME: get this tolerance form env2
   double zerotol = 1e-7;

   while (counter < num_rhs_to_change) {
      if (fabs(rhs[counter]) > zerotol) {
         temp = ((rand_seq[counter]) ? (1 + double(percent_increase/100)) : (1 - double(percent_decrease/100)));
         if (newsense[counter] == 'L' && !newrange[counter]) {
            sym_set_row_upper(env2, counter, temp * rhs[counter]);
         } else if (newsense[counter] == 'G' && !newrange[counter]) {
            sym_set_row_lower(env2, counter, temp * rhs[counter]);
         } else if (newsense[counter] == 'E') {
            sym_set_row_lower(env2, counter, temp * rhs[counter]);
            sym_set_row_upper(env2, counter, temp * rhs[counter]);
         } else {
            // Ranged constraint!
            printf("main(): Unsupported row sense during perturbation!\n");
            exit(1);
         }
      }
      counter++;
   }
   /* Perturbation END */

   /* Write out the problem at hand */
   char *file_name3 = (char *) malloc(CSIZE * 80);
   sprintf(file_name3, "%s_modPerturbed", file_name);
   sym_write_mps(env2, file_name3);

   /* Finding cut coefficients BEGIN */
   // Getting some data from previous tree
   int num_leaf_nodes;
   sym_get_num_leaves(env, &num_leaf_nodes);
   int *lb_cnt = (int *) malloc(ISIZE * num_leaf_nodes);
   int *ub_cnt = (int *) malloc(ISIZE * num_leaf_nodes);
   int **lb_ind = (int **) malloc(sizeof(int *) * num_leaf_nodes);
   int **ub_ind = (int **) malloc(sizeof(int *) * num_leaf_nodes);
   double **lb_val = (double **) malloc(sizeof(double *) * num_leaf_nodes);
   double **ub_val = (double **) malloc(sizeof(double *) * num_leaf_nodes);
   sym_get_branchdesc_bounds(env, lb_cnt, lb_ind, lb_val, ub_cnt, ub_ind, ub_val);
   CoinPackedMatrix duals_by_row = CoinPackedMatrix();
   sym_get_leaf_duals_by_row(env, duals_by_row);
   CoinPackedMatrix pos_djs_by_row = CoinPackedMatrix();
   sym_get_leaf_pos_djs_by_row(env, pos_djs_by_row);
   CoinPackedMatrix neg_djs_by_row = CoinPackedMatrix();
   sym_get_leaf_neg_djs_by_row(env, neg_djs_by_row);

   struct timeval begin_tval, end_tval;
   long elapsed_tval;
   // Finding cut LHS
   // Getting constraint matrix from base instance
   //    NOTE: all constraints are in <= form due to SYMPHONY's internal change
   gettimeofday(&begin_tval, NULL);
   printf("About to find cut LHS\n");
   int nz;
   sym_get_num_elements(env, &nz);
   int *matbeg = (int *) malloc(ISIZE * (numcols+1));
   int *matind = (int *) malloc(ISIZE * nz);
   double *matval = (double *) malloc(DSIZE * nz);
   sym_get_matrix(env, &nz, matbeg, matind, matval);

   double *cutlhs = (double *) calloc(DSIZE, numcols);
   int *numelems = (int *) malloc(ISIZE * numcols);
   int **indices = (int **) malloc(sizeof(int *) * numcols);
   double **values = (double **) malloc(sizeof(double *) * numcols);

   for (i = 0; i < numcols; i++) {
      numelems[i] = matbeg[i+1] - matbeg[i];
      indices[i] = (int *) malloc(ISIZE * numelems[i]);
      values[i] = (double *) malloc(DSIZE * numelems[i]);
      memcpy(indices[i], &matind[matbeg[i]], ISIZE * numelems[i]);
      memcpy(values[i], &matval[matbeg[i]], DSIZE * numelems[i]);
   }
   int colind;
   double col_val = 1;

   for (i = 0; i < numcols; i++) {
      colind = i;
      sym_get_coeff_for_new_rhs(env, numelems[i], indices[i], values[i],
        1, &colind, &col_val, 1, &colind, &col_val, &cutlhs[i]);
   }

   int newrownz = 0;
   for (i = 0; i < numcols; i++) {
      if (fabs(cutlhs[i]) > zerotol)
         newrownz++;
   }

   counter = 0;
   int *newrowind = (int *) malloc(ISIZE * newrownz);
   double *newrowval = (double *) malloc(DSIZE * newrownz);
   for (i = 0; i < numcols; i++) {
      if (fabs(cutlhs[i]) > zerotol) {
         newrowind[counter] = i;
         newrowval[counter] = cutlhs[i];
         counter++;
      }
   }
   assert(counter == newrownz);
   printf("Found cut LHS\n\n");
   gettimeofday(&end_tval, NULL);
   elapsed_tval = ((end_tval.tv_sec - begin_tval.tv_sec)*1000000L + 
            (end_tval.tv_usec - begin_tval.tv_usec));
   printf("Microsec time for finding cut LHS = %ld\n\n", elapsed_tval);

   // Finding cut RHS
   printf("About to find cut RHS\n");
   gettimeofday(&begin_tval, NULL);

   sym_get_rhs(env2, rhs);
   int *rhs_ind = (int *) malloc(ISIZE * numrows);
   int *col_ind = (int *) malloc(ISIZE * numcols);
   for (i = 0; i < numrows; i++) {
      rhs_ind[i] = i;
   }
   for (i = 0; i < numcols; i++) {
      col_ind[i] = i;
   }
   // Changing rhs such that it represents all "<=" type cons
   // coz SYMPHONY changed so in base instance solving!
   for (i = 0; i < numrows; i++) {
      if (newsense[i] == 'G' && !newrange[i]) {
         rhs[i] *= -1;
      } else if ((newsense[i] == 'R') || (newsense[i] == 'G' && newrange[i])
            || (newsense[i] == 'L' && newrange[i])) {
         // Ranged constraint!
         printf("main(): Unsupported row sense!\n");
         exit(1);
      }
   }

   double *newlb_val = (double *) malloc(DSIZE * numcols);
   double *newub_val = (double *) malloc(DSIZE * numcols);
   sym_get_col_lower(env, newlb_val);
   sym_get_col_upper(env, newub_val);

   double cutrhs = SYM_INFINITY;
   double *templb_val = (double *) malloc(DSIZE * numcols);
   double *tempub_val = (double *) malloc(DSIZE * numcols);
   double temprhs;
   int k;

   for (i = 0; i < num_leaf_nodes; i++) {
      memcpy(templb_val, newlb_val, DSIZE * numcols);
      memcpy(tempub_val, newub_val, DSIZE * numcols);
      for (k = 0; k < lb_cnt[i]; k++) {
         if (lb_val[i][k] >= templb_val[lb_ind[i][k]] + zerotol) {
            templb_val[lb_ind[i][k]] = lb_val[i][k];
         }
      }
      for (k = 0; k < ub_cnt[i]; k++) {
         if (ub_val[i][k] <= tempub_val[ub_ind[i][k]] - zerotol) {
            tempub_val[ub_ind[i][k]] = ub_val[i][k];
         }
      }

      sym_get_coeff_for_new_rhs(env, numrows, rhs_ind, rhs,
            numcols, col_ind, templb_val, numcols, col_ind, tempub_val, &temprhs);

      if (temprhs < cutrhs)
         cutrhs = temprhs;
   }
   printf("Found cut RHS = %f\n\n", cutrhs);
   gettimeofday(&end_tval, NULL);
   elapsed_tval = ((end_tval.tv_sec - begin_tval.tv_sec)*1000000L + 
            (end_tval.tv_usec - begin_tval.tv_usec));
   printf("Microsec time for finding cut RHS = %ld\n\n", elapsed_tval);
  /* Finding cut coefficients END */

   sym_add_row(env2, newrownz, newrowind, newrowval, 'G', cutrhs, 0);

   /* Write out the problem at hand */
   char *file_name4 = (char *) malloc(CSIZE * 80);
   sprintf(file_name4, "%s_modAdded", file_name);
   sym_write_mps(env2, file_name4);
   sym_write_lp(env2, file_name4);

   sym_set_int_param(env2, "keep_warm_start", TRUE);
   sym_set_int_param(env2, "warm_start_type", 1);
   sym_set_int_param(env2, "sensitivity_analysis", TRUE);
   sym_set_int_param(env2, "sensitivity_rhs", TRUE);
   sym_set_int_param(env2, "sensitivity_bounds", TRUE);
   sym_set_int_param(env2, "do_reduced_cost_fixing", 0);
   sym_set_int_param(env2, "do_primal_heuristic", 0);
   sym_set_int_param(env2, "generate_cgl_cuts", 0);
   sym_set_int_param(env2, "should_use_rel_br", 0);
   sym_set_int_param(env2, "prep_level", 0);
   sym_set_int_param(env2, "node_limit", 500);

   sym_solve(env2);

   // Freeing memory
   free(file_name4);
   free(tempub_val);
   free(templb_val);
   free(newub_val);
   free(newlb_val);
   free(col_ind);
   free(rhs_ind);
   free(newrowval);
   free(newrowind);
   for (i = 0; i < numcols; i++) {
      free(values[i]);
      free(indices[i]);
   }
   free(values);
   free(indices);
   free(numelems);
   free(cutlhs);
   free(matval);
   free(matind);
   free(matbeg);
   for (i = 0; i < num_leaf_nodes; i++) {
      free(ub_val[i]);
      free(lb_val[i]);
      free(ub_ind[i]);
      free(lb_ind[i]);
   }
   free(ub_val);
   free(lb_val);
   free(ub_ind);
   free(lb_ind);
   free(ub_cnt);
   free(lb_cnt);
   free(file_name3);
   free(rhs);
   free(newrange);
   free(newsense);
   free(rand_seq);
   free(file_name1);

   sym_close_environment(env2);
   sym_close_environment(env);

   return(0);
}

#endif
