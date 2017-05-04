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

   /* Read parameters specific to UC WS: percent_decrease, percent_increase, srand_seed */
   // default values
   double percent_increase = 10.0;
   double percent_decrease = 10.0;
   unsigned int srand_seed = 1;
   int i, j;
   for (i = 0; i < argc; i++) {
      if (strcmp("--up_change", argv[i]) == 0)
         sscanf(argv[i+1], "%lf", &percent_increase);
      else if (strcmp("--down_change", argv[i]) == 0)
         sscanf(argv[i+1], "%lf", &percent_decrease);
      else if (strcmp("--srand_seed", argv[i]) == 0)
         sscanf(argv[i+1], "%u", &srand_seed);
   }

   sym_load_problem(env);

   int numcols, var_ind;
   double var_bound, var_coeff;
   sym_get_num_cols(env, &numcols);

   // Add variable bounds as constraints
   for (i = 0; i < numcols; i++) {
      var_ind = i;
      var_bound = 0;
      var_coeff = 1;
      sym_add_row(env, 1, &var_ind, &var_coeff, 'G', var_bound, 0);
   }
   for (i = 0; i < numcols; i++) {
      var_ind = i;
      var_bound = 1;
      var_coeff = 1;
      sym_add_row(env, 1, &var_ind, &var_coeff, 'L', var_bound, 0);
   }
   
   sym_set_int_param(env, "keep_warm_start", TRUE);
   sym_set_int_param(env, "warm_start_type", 1);
   sym_set_int_param(env, "sensitivity_analysis", TRUE);
   sym_set_int_param(env, "sensitivity_rhs", TRUE);
   sym_set_int_param(env, "sensitivity_bounds", TRUE);
//   sym_set_int_param(env, "do_reduced_cost_fixing", 0);
   sym_set_int_param(env, "do_primal_heuristic", 0);
   sym_set_int_param(env, "generate_cgl_cuts", 0);
   sym_set_int_param(env, "should_use_rel_br", 0);
   sym_set_int_param(env, "prep_level", 0);
//   sym_set_int_param(env, "verbosity", 5);
//   sym_set_dbl_param(env, "gap_limit", 0.25);
   sym_set_int_param(env, "node_limit", 20);

   /* Write out the problem at hand */
   char *file_name = (char *) malloc(CSIZE * 80);
   sym_get_str_param(env, "problem_name", &file_name);

   sym_solve(env);

   sym_environment *env2 = sym_open_environment();   
   sym_parse_command_line(env2, argc, argv);   
   sym_load_problem(env2);

   /* Perturbation */
   int numrows = 0;
   int num_rhs_to_change, counter = 0;
   double temp;

   sym_get_num_rows(env, &numrows);
   int *rand_seq = (int *) malloc(ISIZE * numrows);
   double *rhs = (double *) malloc(DSIZE * numrows);
   sym_get_rhs(env, rhs);

   /* Generate a sequence of random binary digits to decide which rhs to 
    * increase (1) and which demand to decrese (0) */
   srand(srand_seed);
   for (i = 0; i < numrows; i++) {
      rand_seq[i] = rand()%2;
   }

   num_rhs_to_change = numrows;

   while (counter < num_rhs_to_change) {
      if (rhs[counter]) {
         temp = ((rand_seq[counter]) ? (1 + double(percent_increase/100)) : (1 - double(percent_decrease/100)));
         sym_set_row_upper(env2, counter, temp * rhs[counter]);
      }
      counter++;
   }

   int nz;
   sym_get_num_elements(env, &nz);
   int *matbeg = (int *) malloc(ISIZE * (numcols+1));
   int *matind = (int *) malloc(ISIZE * nz);
   double *matval = (double *) malloc(DSIZE * nz);
   sym_get_matrix(env, &nz, matbeg, matind, matval);
   double *newrhsval = (double *) malloc(DSIZE * numrows);
   sym_get_rhs(env2, newrhsval);
   int *newrhsind = (int *) malloc(ISIZE * numrows);
   int newrownz = 0;
   double zerotol = 1e-7;
   for (i = 0; i < numrows; i++) {
      newrhsind[i] = i;
   }


   double cutrhs = SYM_INFINITY;
   int *numelems = (int *) malloc(ISIZE * numcols);
   int **indices = (int **) malloc(sizeof(int *) * numcols);
   double **values = (double **) malloc(sizeof(double *) * numcols);
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
   double *temp_val = (double *) malloc(DSIZE * numcols);
   sym_get_branchdesc_bounds(env, lb_cnt, lb_ind, lb_val, ub_cnt, ub_ind, ub_val);

   printf("About to find cut RHS\n");
   for (i = 0; i < num_leaf_nodes; i++) {
      if (ub_cnt[i]) {
         memcpy(temp_val, ub_val[i], DSIZE*numcols);
         for (j = 0; j < numcols; j++) {
            temp_val[j] *= -1;
         }
      } else {
         temp_val = NULL;
      }
      sym_get_coeff_for_new_rhs(env, numrows, newrhsind, newrhsval, 
         lb_cnt[i], lb_ind[i], lb_val[i], ub_cnt[i], ub_ind[i], temp_val, &temp);
      if (temp < cutrhs)
         cutrhs = temp;
   }
   printf("Found cut RHS = %f", cutrhs);
   double *cutlhs = (double *) calloc(DSIZE, numcols);
   printf("\n");
   for (i = 0; i < numcols; i++) {
      numelems[i] = matbeg[i+1] - matbeg[i];
      indices[i] = (int *) malloc(ISIZE * numelems[i]);
      values[i] = (double *) malloc(DSIZE * numelems[i]);
      memcpy(indices[i], &matind[matbeg[i]], ISIZE * numelems[i]);
      memcpy(values[i], &matval[matbeg[i]], DSIZE * numelems[i]);
   }
   int lbcnt = 1, ubcnt = 1;
   int new_lb_ind, new_ub_ind;
   double new_lb_val = 1, new_ub_val = -1;
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
         newrowval[counter] = -cutlhs[i];
         counter++;
      }
   }

   /* Write out the problem at hand */
   char *file_name2 = (char *) malloc(CSIZE * 80);
   sprintf(file_name2, "%s_mod1", file_name);
   sym_write_mps(env2, file_name2);

   sym_add_row(env2, newrownz, newrowind, newrowval, 'L', -cutrhs, 0);

   /* Write out the problem at hand */
   char *file_name3 = (char *) malloc(CSIZE * 80);
   sprintf(file_name3, "%s_mod2", file_name);
   sym_write_mps(env2, file_name3);

   sym_solve(env2);

   sym_close_environment(env2);
   sym_close_environment(env);
  
   return(0);
}

#endif
