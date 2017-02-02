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
   
int main(int argc, char **argv)
{    
   sym_environment *env = sym_open_environment();   
   sym_parse_command_line(env, argc, argv);   
   sym_load_problem(env);
   
   sym_set_int_param(env, "keep_warm_start", TRUE);
//   sym_set_int_param(env, "do_reduced_cost_fixing", 1);
//   sym_set_int_param(env, "do_primal_heuristic", 1);
//   sym_set_int_param(env, "generate_cgl_cuts", 1);
//   sym_set_int_param(env, "generate_cgl_gomory_cuts", -1);
//   sym_set_int_param(env, "generate_cgl_knapsack_cuts", -1);
//   sym_set_int_param(env, "generate_cgl_twomir_cuts", -1);
   sym_set_int_param(env, "prep_level", 0);
//   sym_set_int_param(env, "verbosity", 5);
   sym_set_dbl_param(env, "gap_limit", 0);

   /* Write out the problem at hand */
   char *file_name = (char *) malloc(CSIZE * 80);
   sym_get_str_param(env, "problem_name", &file_name);

   /*
   char *file_name2 = (char *) malloc(CSIZE * 80);
   sprintf(file_name2, "%s_basic", file_name);
   sym_write_mps(env, file_name2);
   sym_write_lp(env, file_name2);
   */

   sym_solve(env);

   sym_set_row_upper(env, 0, 6000000);

   sym_set_int_param(env, "prep_level", 0);
//   sym_set_int_param(env, "do_primal_heuristic", 0);
//   sym_set_int_param(env, "generate_cgl_cuts", 1);
//   sym_set_int_param(env, "generate_cgl_gomory_cuts", 3);
//   sym_set_int_param(env, "generate_cgl_knapsack_cuts", 3);
//   sym_set_int_param(env, "generate_cgl_probing_cuts", 3);
//   sym_set_int_param(env, "generate_cgl_flowcover_cuts", 3);
//   sym_set_int_param(env, "generate_cgl_twomir_cuts", 3);
//   sym_set_int_param(env, "generate_cgl_clique_cuts", 3);
//   sym_set_int_param(env, "warm_start_node_level", 1);
//   sym_set_dbl_param(env, "warm_start_node_level_ratio", 1.0);
//   sym_set_int_param(env, "debug_lp", 1);
//   sym_set_int_param(env, "verbosity", 5);
   sym_set_dbl_param(env, "gap_limit", 0.0);

   /* Write out the problem at hand */
   char *file_name3 = (char *) malloc(CSIZE * 80);
   sprintf(file_name3, "%s_mod", file_name);
//   sym_write_mps(env, file_name3);
//   sym_write_lp(env, file_name3);

   sym_warm_solve(env);

   sym_close_environment(env);
  
   return(0);
}  

#endif
