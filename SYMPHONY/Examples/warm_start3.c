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

#include <cstdio>

#ifdef USE_OSI_INTERFACE

#include "OsiSymSolverInterface.hpp"
#include <iostream>
int main(int argc, char **argv)
{
   int termcode;
   OsiSymSolverInterface si;

   si.parseCommandLine(argc, argv);
   si.loadProblem();

   si.setSymParam(OsiSymTimeLimit, 10);
   si.setSymParam(OsiSymKeepWarmStart, 1);
   si.initialSolve();
   termcode = si.isProvenOptimal();

   while (!termcode){
      printf("Starting problem again from warm start...\n");
      si.resolve();
      termcode = si.isProvenOptimal();
   }

   return(0);
}

#else

#include "symphony.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
  
int main(int argc, char **argv)
{
   sym_environment *env = sym_open_environment();   
   sym_parse_command_line(env, argc, argv);   
   sym_load_problem(env);

   sym_set_int_param(env, "keep_warm_start", TRUE);
   sym_set_int_param(env, "do_reduced_cost_fixing", 0);
   sym_set_int_param(env, "do_primal_heuristic", 0);
   sym_set_int_param(env, "generate_cgl_cuts", 0);
   sym_set_int_param(env, "prep_level", 0);
   sym_set_dbl_param(env, "gap_limit", 0.0);
   sym_set_int_param(env, "should_use_rel_br", 0);
   sym_set_int_param(env, "sensitivity_analysis", TRUE);
   sym_set_int_param(env, "sensitivity_rhs", TRUE);
   sym_solve(env);

   // ***** NOTE: parameters to change: card_G, card_T, card_N, 
          //  percent_increase, percent_decrease, up_reserve, down_reserve ***** //
   /* modify RHS values of sup_demand, reserve_up, reserve_down, flow_bal cons randomly */
   int card_G = 5, card_T = 24, card_N = 5;
   double percent_increase = 1.0, percent_decrease = 1.0, up_reserve = 0.10, down_reserve = 0.05;

   int num_rows = 0, i, j, total_demands = card_N * card_T;
   int *rand_seq = (int *) malloc(ISIZE * total_demands);
   int r1, num_rhs_to_change, counter = 0;
   unsigned int srand_seed = 10;

   sym_get_num_rows(env, &num_rows);
   double *rhs = (double *) malloc(DSIZE * num_rows);
   sym_get_rhs(env, rhs);

   /* Total demands: rhs[6*card_G*card_T + 3*card_T] to rhs[6*card_G*cardT + 3*card_T + total_demands].
   * Generate a sequence of random binary digits to decide which demand to 
   * increase (1) and which demand to decrese (0) */
   srand(srand_seed);
   for (i = 0; i < total_demands; i++) {
      rand_seq[i] = rand()%2;
   }

   // Generate new_rhs_card_T based on rand_seq
   double *new_rhs_card_T = (double *) calloc(card_T, DSIZE);
   double temp2;
   int temp = 6*card_G*card_T + 3*card_T;
   for (i = 0; i < card_T; i++) {
      for (j = 0; j < card_N; j++) {
         new_rhs_card_T[i] += rhs[temp + j*card_T + i] * 
            ((rand_seq[j*card_T + i]) ? (1 + double(percent_increase/100)) : (1 - double(percent_decrease/100)));
      }
   }

   sym_environment *env2 = sym_open_environment();   
   sym_parse_command_line(env2, argc, argv);   
   sym_load_problem(env2);

   num_rhs_to_change = (3)*card_T + card_N*card_T;

   while (counter < num_rhs_to_change) {
      // index of con for rhs change
      r1 = 6*card_G*card_T + counter;
      if (counter < card_T) {
         sym_set_row_upper(env2, r1, new_rhs_card_T[counter]);
         sym_set_row_lower(env2, r1, new_rhs_card_T[counter]);
      } else if (counter < 2*card_T) {
         sym_set_row_upper(env2, r1, -new_rhs_card_T[counter - card_T]*(1 + up_reserve));
      } else if (counter < 3*card_T) {
         sym_set_row_upper(env2, r1, new_rhs_card_T[counter - 2*card_T]*(1 - down_reserve));
      } else {
         temp2 = ((rand_seq[counter - 3*card_T]) ? (1 + double(percent_increase/100)) : (1 - double(percent_decrease/100)));
         sym_set_row_upper(env2, r1, temp2 * rhs[r1]);
         sym_set_row_lower(env2, r1, temp2 * rhs[r1]);
      }
      counter++;
   }

   int nz, numcols;
   sym_get_num_elements(env, &nz);
   sym_get_num_cols(env, &numcols);
   int *matbeg = (int *) malloc(ISIZE * (numcols+1));
   int *matind = (int *) malloc(ISIZE * nz);
   double *matval = (double *) malloc(DSIZE * nz);
   sym_get_matrix(env, &nz, matbeg, matind, matval);
   double *newrhs = (double *) malloc(DSIZE * num_rows);
   sym_get_rhs(env2, newrhs);
   int newrhsnz = 0;
   double zerotol = 1e-7;
   for (i = 0; i < num_rows; i++) {
      if (fabs(newrhs[i]) > zerotol)
         newrhsnz++;
   }
   int *newrhsind = (int *) malloc(ISIZE * newrhsnz);
   double *newrhsval = (double *) malloc(DSIZE * newrhsnz);
   counter = 0;
   for (i = 0; i < num_rows; i++) {
      assert(counter <= newrhsnz);
      if (fabs(newrhs[i]) > zerotol) {
         newrhsind[counter] = i;
         newrhsval[counter] = newrhs[i];
         counter++;
      }
   }

   printf("\n About to start finding cut LHS and RHS \n");
   printf("About to find cut RHS\n");
   double cutrhs;
   sym_get_coeff_for_new_rhs(env, newrhsnz, newrhsind, newrhsval, 0, NULL, NULL, 0, NULL, NULL, &cutrhs);
   printf("Found cut RHS = %f\n", cutrhs);

   int newrownz = 0;
   int *numelems = (int *) malloc(ISIZE * numcols);
   int **indices = (int **) malloc(sizeof(int *) * numcols);
   double **values = (double **) malloc(sizeof(double *) * numcols);
   double *cutlhs = (double *) calloc(DSIZE, numcols);
   for (i = 0; i < numcols; i++) {
      numelems[i] = matbeg[i+1] - matbeg[i];
      indices[i] = (int *) malloc(ISIZE * numelems[i]);
      values[i] = (double *) malloc(DSIZE * numelems[i]);
      memcpy(indices[i], &matind[matbeg[i]], ISIZE * numelems[i]);
      memcpy(values[i], &matval[matbeg[i]], DSIZE * numelems[i]);
   }
#pragma omp parallel for
   for (i = 0; i < numcols; i++) {
      printf("About to find var %d coeff\n", i);
      sym_get_coeff_for_new_rhs(env, numelems[i], indices[i], values[i], 0, NULL, NULL, 0, NULL, NULL, &cutlhs[i]);
      printf("var %d, coeff = %f\n", i, cutlhs[i]);
   }
   printf("\n Done finding cut LHS and RHS \n");
   for (i = 0; i < numcols; i++) {
      if (fabs(cutlhs[i]) > zerotol)
         newrownz++;
   }
   int *newrowind = (int *) malloc(ISIZE * newrownz);
   double *newrowval = (double *) malloc(DSIZE * newrownz);
   counter = 0;
   for (i = 0; i < numcols; i++) {
      assert(counter <= newrownz);
      if (fabs(cutlhs[i]) > zerotol) {
         newrowind[counter] = i;
         newrowval[counter] = cutlhs[i];
         counter++;
      }
   }
   sym_add_row(env2, newrownz, newrowind, newrowval, 'L', cutrhs, 0);
   printf("\n Done adding the new cut \n");

   /* solve the modified problem with warm starting */
//   sym_set_int_param(env, "prep_level", 0);
//   sym_set_int_param(env, "do_reduced_cost_fixing", 0);
//   sym_set_int_param(env, "do_primal_heuristic", 0);
//   sym_set_int_param(env, "generate_cgl_cuts", 0);
//   sym_set_dbl_param(env, "gap_limit", 0.0);

   /* Write out the problem at hand */
   char *file_name = (char *) malloc(CSIZE * 80);
   sym_get_str_param(env2, "problem_name", &file_name);
   sprintf(file_name, "%s_%u_%.1f_%.1f", file_name, srand_seed, percent_increase, percent_decrease);
   sym_write_mps(env2, file_name);

   sym_solve(env2);

   sym_close_environment(env2);
   sym_close_environment(env);

   return(0);
}

#endif
