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
#include <time.h>

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

   time_t begin_time, end_time;
   long elapsed_time;
   clock_t begin_tick, end_tick;
   long elapsed_tick;
   begin_time = time(NULL);
   begin_tick = clock();
   // Finding cut LHS
   // Getting constraint matrix from base instance
   //    NOTE: all constraints are in <= form due to SYMPHONY's internal change
   printf("About to find cut LHS\n");
   int nz;
   sym_get_num_elements(env, &nz);
   int *matbeg = (int *) malloc(ISIZE * (numcols+1));
   int *matind = (int *) malloc(ISIZE * nz);
   double *matval = (double *) malloc(DSIZE * nz);
   sym_get_matrix(env, &nz, matbeg, matind, matval);

   int *eye_matbeg = (int *) malloc(ISIZE * (numcols+1));
   int *eye_matind = (int *) malloc(ISIZE * numcols);
   double *eye_matval = (double *) malloc(DSIZE * numcols);
   eye_matbeg[0] = 0;
   for (i = 0; i < numcols; i++) {
      eye_matbeg[i+1] = eye_matbeg[i] + 1;
      eye_matind[i] = i;
      eye_matval[i] = 1;
   }
   double *cutlhs = (double *) calloc(DSIZE, numcols);

   //last argument: 1 ==> cut LHS
   //               2 ==> cut RHS
   //TODO: work around the last argument later! It is here right now because
   //      in case of RHS, only "b" vector is passed on instead of a
   //      [b..b 'num_leaf_nodes times'] matrix to avoid redundant multiplication 
   //      operations with the duals matrix
   sym_get_coeff_for_new_rhs(env, matbeg, matind, matval,
      eye_matbeg, eye_matind, eye_matval, eye_matbeg, eye_matind, eye_matval,
      cutlhs, numcols, 1);

   int newrownz = 0;
   for (i = 0; i < numcols; i++) {
      if (fabs(cutlhs[i]) > zerotol)
         newrownz++;
   }

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
   printf("Found cut LHS\n\n");

   // Finding cut RHS
   sym_get_rhs(env2, rhs);
   int *rhs_ind = (int *) malloc(ISIZE * numrows);
   for (i = 0; i < numrows; i++) {
      rhs_ind[i] = i;
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

   printf("About to find cut RHS\n");
   int k;
   //matrices for LBs and UBs stored in arrays representing column ordered format
   double *templb_mat = (double *) malloc(DSIZE * numcols*num_leaf_nodes);
   double *tempub_mat = (double *) malloc(DSIZE * numcols*num_leaf_nodes);

   int *lb_matbeg = (int *) malloc(ISIZE * (num_leaf_nodes+1));
   int *ub_matbeg = (int *) malloc(ISIZE * (num_leaf_nodes+1));
   lb_matbeg[0] = ub_matbeg[0] = 0;
   int lb_nz, ub_nz;

   for (i = 0; i < num_leaf_nodes; i++) {
      memcpy(&templb_mat[i*numcols], &newlb_val[0], DSIZE*numcols);
      memcpy(&tempub_mat[i*numcols], &newub_val[0], DSIZE*numcols);
      for (k = 0; k < lb_cnt[i]; k++) {
         if (lb_val[i][k] >= templb_mat[i*numcols + lb_ind[i][k]]) {
            templb_mat[i*numcols + lb_ind[i][k]] = lb_val[i][k];
         }
      }
      for (k = 0; k < ub_cnt[i]; k++) {
         if (ub_val[i][k] <= tempub_mat[i*numcols + ub_ind[i][k]]) {
            tempub_mat[i*numcols + ub_ind[i][k]] = ub_val[i][k];
         }
      }
      lb_nz = ub_nz = 0;
      for (k = 0; k < numcols; k++) {
         if (fabs(templb_mat[i*numcols + k]) > zerotol) {
            lb_nz++;
         }
         if (fabs(tempub_mat[i*numcols + k]) > zerotol) {
            ub_nz++;
         }
      }
      lb_matbeg[i+1] = lb_matbeg[i] + lb_nz;
      ub_matbeg[i+1] = ub_matbeg[i] + ub_nz;
   }
   
   int *lb_matind = (int *) malloc(ISIZE * lb_matbeg[num_leaf_nodes]);
   int *ub_matind = (int *) malloc(ISIZE * ub_matbeg[num_leaf_nodes]);
   double *lb_matval = (double *) malloc(DSIZE * lb_matbeg[num_leaf_nodes]);
   double *ub_matval = (double *) malloc(DSIZE * ub_matbeg[num_leaf_nodes]);
   lb_nz = ub_nz = 0;
   for (i = 0; i < num_leaf_nodes; i++) {
      for (k = 0; k < numcols; k++) {
         if (fabs(templb_mat[i*numcols + k]) > zerotol) {
            lb_matind[lb_nz] = k;
            lb_matval[lb_nz] = templb_mat[i*numcols + k];
            lb_nz++;
         }
         if (fabs(tempub_mat[i*numcols + k]) > zerotol) {
            ub_matind[ub_nz] = k;
            ub_matval[ub_nz] = tempub_mat[i*numcols + k];
            ub_nz++;
         }
      }
   }

   double *cutrhs_array = (double *) calloc(DSIZE, num_leaf_nodes);
   sym_get_coeff_for_new_rhs(env, NULL, rhs_ind, rhs,
         lb_matbeg, lb_matind, lb_matval, ub_matbeg, ub_matind, ub_matval,
         cutrhs_array, num_leaf_nodes, 2);
   double cutrhs = SYM_INFINITY;
   for (i = 0; i < num_leaf_nodes; i++) {
      if (cutrhs_array[i] < cutrhs) {
         cutrhs = cutrhs_array[i];
      }
   }
   printf("Found cut RHS = %f\n\n", cutrhs);
   /* Finding cut coefficients END */

   sym_add_row(env2, newrownz, newrowind, newrowval, 'G', cutrhs, 0);
   end_tick = clock();
   end_time = time(NULL);
   elapsed_time = (long) (end_time - begin_time);
   elapsed_tick = (long) (end_tick - begin_tick) / CLOCKS_PER_SEC;
   printf("Time for finding cut = %ld\n\n", elapsed_time);
   printf("CPU time for finding cut = %ld\n\n", elapsed_tick);

   /* Write out the problem at hand */
   char *file_name4 = (char *) malloc(CSIZE * 80);
   sprintf(file_name4, "%s_modAdded", file_name);
   sym_write_mps(env2, file_name4);

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
   free(cutrhs_array);
   free(ub_matval);
   free(lb_matval);
   free(ub_matind);
   free(lb_matind);
   free(ub_matbeg);
   free(lb_matbeg);
   free(tempub_mat);
   free(templb_mat);
   free(newub_val);
   free(newlb_val);
   free(rhs_ind);
   free(newrowval);
   free(newrowind);
   free(cutlhs);
   free(eye_matval);
   free(eye_matind);
   free(eye_matbeg);  
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
