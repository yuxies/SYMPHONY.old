/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* the Vehicle Routing Problem and the Traveling Salesman Problem.           */
/*                                                                           */
/* (c) Copyright 2000-2003 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This application was developed by Ted Ralphs (tkralphs@lehigh.edu)        */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#define COMPILING_FOR_MASTER

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the main() for the master process.
\*===========================================================================*/

#if 1

#include "OsiSymSolverInterface.hpp"

int main(int argc, char **argv)
{
   OsiSymSolverInterface si;

   si.setSymParam(OsiSymVerbosity, 3);
   si.setSymParam(OsiSymGranularity, 0.9999); 
   //si.setSymParam(OsiSymGenerateMipCuts, false);
   
   /* Parse the command line */
   si.parseCommandLine(argc, argv);
   
   /* Read in the problem */
   si.loadProblem();

   /* Find a priori problem bounds */
   si.findInitialBounds();

   /* Solve the problem */
   si.branchAndBound();

   return(0);
}

#else

#include "master.h"

int main(int argc, char **argv)
{
   problem *p = sym_open_environment();

   sym_parse_command_line(p, argc, argv);

   sym_load_problem(p);

   sym_find_initial_bounds(p);

   sym_solve(p);

   sym_close_environment(p);

   return(0);
}

#endif