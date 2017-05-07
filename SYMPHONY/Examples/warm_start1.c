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
   sym_set_int_param(env, "verbosity", 5);
   sym_set_int_param(env, "node_limit", 20);

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

   /* Write out the problem at hand */
   char *file_name2 = (char *) malloc(CSIZE * 80);
   sprintf(file_name2, "%s_modImmediate", file_name);
   sym_write_mps(env2, file_name2);

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
   int *leaf_depth = (int *) malloc(ISIZE * num_leaf_nodes);
   sym_get_leaf_depths(env, leaf_depth);
   int *lb_cnt = (int *) malloc(ISIZE * num_leaf_nodes);
   int *ub_cnt = (int *) malloc(ISIZE * num_leaf_nodes);
   int **lb_ind = (int **) malloc(sizeof(int *) * num_leaf_nodes);
   int **ub_ind = (int **) malloc(sizeof(int *) * num_leaf_nodes);
   double **lb_val = (double **) malloc(sizeof(double *) * num_leaf_nodes);
   double **ub_val = (double **) malloc(sizeof(double *) * num_leaf_nodes);
   sym_get_branchdesc_bounds(env, lb_cnt, lb_ind, lb_val, ub_cnt, ub_ind, ub_val);

   // Finding cut RHS
   sym_get_rhs(env2, rhs);
   int *rhs_ind = (int *) malloc(ISIZE * numrows);
   for (i = 0; i < numrows; i++) {
      rhs_ind[i] = i;
   }

   double *newlb_val = (double *) malloc(DSIZE * numcols);
   double *newub_val = (double *) malloc(DSIZE * numcols);
   sym_get_col_lower(env, newlb_val);
   sym_get_col_upper(env, newub_val);

   int nzlb_num = 0, nzub_num = 0;

   printf("About to find cut RHS\n");
   double cutrhs = SYM_INFINITY;
   int j;
   for (i = 0; i < num_leaf_nodes; i++) {
      for (j = 0; j < lb_cnt[i]; j++) {
         if (lb_val[i][j] >= newlb_val[lb_ind[i][j]]) {
            newlb_val[lb_ind[i][j]] = lb_val[i][j];
         }
      }
      for (j = 0; j < ub_cnt[i]; j++) {
         if (ub_val[i][j] <= newub_val[ub_ind[i][j]]) {
            newub_val[ub_ind[i][j]] = ub_val[i][j];
         }
      }
      nzlb_num = 0;
      nzub_num = 0;
      // Compressing lb and ub to capture only nz bounds
      //    NOTE: TODO: SYM_INFINITY bounds still exist?
      for (j = 0; j < numcols; j++) {
         if (fabs(newlb_val[j]) > zerotol) {
            nzlb_num++;
         }
         if (fabs(newub_val[j] > zerotol)) {
            nzub_num++;
         }
      }
      int *nzlb_ind = (int *) malloc(ISIZE * nzlb_num);
      int *nzub_ind = (int *) malloc(ISIZE * nzub_num);
      double *nzlb_val = (double *) malloc(DSIZE * nzlb_num);
      double *nzub_val = (double *) malloc(DSIZE * nzub_num);
      int counter_lb = 0, counter_ub = 0;
      for (j = 0; j < numcols; j++) {
         assert(counter_lb <= nzlb_num);
         assert(counter_ub <= nzub_num);
         if (fabs(newlb_val[j]) > zerotol) {
            nzlb_ind[counter_lb] = j;
            nzlb_val[counter_lb] = newlb_val[j];
            counter_lb++;
         }
         if (fabs(newub_val[j] > zerotol)) {
            nzub_ind[counter_ub] = j;
            nzub_val[counter_ub] = newub_val[j];
            counter_ub++;
         }
      }

      // Changing rhs such that it represents all "<=" type cons
      // coz SYMPHONY changed such in base instance solving and the dual
      // signs are accordingly!
      for (j = 0; j < numrows; j++) {
         if (newsense[j] == 'G' && !newrange[j]) {
            rhs[j] *= -1;
         } else if ((newsense[j] == 'R') || (newsense[j] == 'G' && newrange[j])
                        || (newsense[j] == 'L' && newrange[j])) {
            // Ranged constraint!
            printf("main(): Unsupported row sense!\n");
            exit(1);
         }
      }
      sym_get_coeff_for_new_rhs(env, numrows, rhs_ind, rhs,
            nzlb_num, nzlb_ind, nzlb_val, nzub_num, nzub_ind, nzub_val, &temp);
      if (temp < cutrhs)
         cutrhs = temp;
      free(nzub_val);
      free(nzlb_val);
      free(nzub_ind);
      free(nzlb_ind);
   }
   printf("Found cut RHS = %f\n", cutrhs);

   // Finding cut LHS
   // Getting constraint matrix from base instance
   //    NOTE: all constraints are in <= form due to SYMPHONY's internal change
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
   int lbcnt = 1, ubcnt = 1;
   int new_lb_ind, new_ub_ind;
   double new_lb_val = 1, new_ub_val = 1;
   int newrownz = 0;
#pragma omp parallel for
   for (i = 0; i < numcols; i++) {
      printf("About to find var %d coeff\n", i);
      new_lb_ind = i;
      new_ub_ind = i;
      sym_get_coeff_for_new_rhs(env, numelems[i], indices[i], values[i],
         lbcnt, &new_lb_ind, &new_lb_val, ubcnt, &new_ub_ind, &new_ub_val, &cutlhs[i]);
      printf("var %d, coeff = %f\n", i, cutlhs[i]);
   }
   for (i = 0; i < numcols; i++) {
      if (fabs(cutlhs[i]) > zerotol)
         newrownz++;
   }
   printf("\n");
   counter = 0;
   int *newrowind = (int *) malloc(ISIZE * newrownz);
   double *newrowval = (double *) malloc(DSIZE * newrownz);
   for (i = 0; i < numcols; i++) {
      assert(counter <= newrownz);
      if (fabs(cutlhs[i]) > zerotol) {
         newrowind[counter] = i;
         newrowval[counter] = cutlhs[i];
         counter++;
      }
   }
   /* Finding cut coefficients END */

   sym_add_row(env2, newrownz, newrowind, newrowval, 'G', cutrhs, 0);

   /* Write out the problem at hand */
   char *file_name4 = (char *) malloc(CSIZE * 80);
   sprintf(file_name4, "%s_modAdded", file_name);
   sym_write_mps(env2, file_name4);

   sym_set_int_param(env2, "node_limit", 20);

   sym_solve(env2);

   sym_close_environment(env2);
   sym_close_environment(env);

   return(0);
}

#endif
