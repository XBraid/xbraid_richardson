
/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. Email xbraid-support@llnl.gov for support.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free Software
 * Foundation) version 2.1 dated February 1999.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 ***********************************************************************EHEADER*/

/**
 * Example:       ex-01-load_balance.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-01
 *
 * Help with:     this is a simple example showing usage of refinement
 *
 * Sample run:    mpirun -np 2 ex-01
 *
 * Description:   solve the scalar ODE 
 *                   u' = lambda u, 
 *                   with lambda=-1 and y(0) = 1
 *                in a very simplified XBraid setting.
 *                
 *                When run with the default 10 time steps and no refinement, the solution is:
 *                $ ./ex-01
 *                $ cat ex-01.out.00*
 *                  1.00000000000000e+00
 *                  6.66666666666667e-01
 *                  4.44444444444444e-01
 *                  2.96296296296296e-01
 *                  1.97530864197531e-01
 *                  1.31687242798354e-01
 *                  8.77914951989026e-02
 *                  5.85276634659351e-02
 *                  3.90184423106234e-02
 *                  2.60122948737489e-02
 *                  1.73415299158326e-02
 *
 *                The load balance option enables temporal load balancing. The load balance
 *                is defined by the solution value. A generic sleep function is used to simulate
 *                the load imbalance. 
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "braid.h"

/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

/* App structure can contain anything, and be named anything as well */
typedef struct _braid_App_struct
{
   int       rank;
   int       lbalance;
   double    tol;
   double    tstart;
} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
   double value;
} my_Vector;

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   double dt;                 /* step size */

   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   dt = tstop-tstart;

   /* Use backward Euler to propagate solution */
   (u->value) = 1./(1. + dt)*(u->value);
 
   int microseconds = (int) abs((1000.0*u->value));
   usleep(microseconds);  

   /* XBraid only accepts weights on level 0.*/
   int level;
   braid_StepStatusGetLevel(status,&level);
   if (level == 0)   
   {
      braid_StepStatusSetWFactor(status,1000.0*u->value);
   }

   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   if (t == 0.0) /* Initial condition */
   {
      (u->value) = 1.0;
   }
   else /* All other time points set to arbitrary value */
   {
      (u->value) = 0.456;
   }
   *u_ptr = u;

   return 0;
}

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v;

   v = (my_Vector *) malloc(sizeof(my_Vector));
   (v->value) = (u->value);
   *v_ptr = v;

   return 0;
}

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u);
   return 0;
}

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{
   (y->value) = alpha*(x->value) + beta*(y->value);
   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   double dot;

   dot = (u->value)*(u->value);
   *norm_ptr = sqrt(dot);
   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int        index;
   char       filename[255];
   FILE      *file;
   double     t, wfactor;
   
   braid_AccessStatusGetT(astatus, &t);
   braid_AccessStatusGetTIndex(astatus, &index);
   braid_AccessStatusGetWFactor(astatus, &wfactor);
   sprintf(filename, "%s.%03d", "ex-01-loadbalance.out", app->rank);
   file = fopen(filename, "a");
   if ( t == app->tstart )
      fprintf(file, "##### (Rank, time, value) ##### \n " ) ;
   fprintf(file, "%d %.14e %.14e  \n", app->rank, t, (u->value));
   fflush(file);
   fclose(file);

   return 0;
}

int
my_BufSize(braid_App          app,
           int                *size_ptr,
           braid_BufferStatus bstatus)
{
   *size_ptr = sizeof(double);
   return 0;
}

int
my_BufPack(braid_App          app,
           braid_Vector       u,
           void               *buffer,
           braid_BufferStatus bstatus)
{
   double *dbuffer = buffer;

   dbuffer[0] = (u->value);
   braid_BufferStatusSetSize( bstatus, sizeof(double) );

   return 0;
}

int
my_BufUnpack(braid_App          app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus bstatus)
{
   double    *dbuffer = buffer;
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   (u->value) = dbuffer[0];
   *u_ptr = u;

   return 0;
}

int 
my_GlobalAccess(braid_App      app,
                MPI_Comm       comm,
                braid_Vector  *u,
                braid_Int     *umask,
                braid_Int      npoints,
                braid_GlobalAccessStatus gstatus)
{
   
   int ilower,iupper,gupper,myid;
   MPI_Comm_rank(comm,&myid);
   braid_GlobalAccessStatusGetInterval(gstatus, &ilower, &iupper, &gupper);
   printf(" %d : %d %d %d \n " , myid, ilower, iupper, gupper );
   return 0;
}


/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   braid_Core    core;
   my_App       *app;
   double        tstart, tstop, tol;
   int           ntime, rank, arg_index, print_usage;
   int           lbalance, comm_min, output;

   /* Define time domain: ntime intervals */
   ntime  = 1000;
   tstart = 0.0;
   tstop  = 5.0;
   print_usage = 0;
   tol = 1.0e-6;
   output = 1;
   lbalance = 0;
   comm_min = 0;
   
   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   /* Parse command line */
   arg_index = 0;
   while( arg_index < argc ){
      if( strcmp(argv[arg_index], "-nt") == 0 ){
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-lbalance") == 0 ){
         arg_index++;
         lbalance = atoi(argv[arg_index++]);
      }      
      else if( strcmp(argv[arg_index], "-comm_min") == 0 ){
         arg_index++;
         comm_min = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-tol") == 0 ){
         arg_index++;
         tol = atof(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-no_output") == 0 ){
         arg_index++;
         output = 0;
      }
      else if( strcmp(argv[arg_index], "-help") == 0 ){
         print_usage = 1;
         break;
      }
      else{
         if(arg_index > 1){
            printf("UNUSED command line paramter %s\n", argv[arg_index]);
         }
         arg_index++;
      }
   }

   if((print_usage) && (rank == 0)){
      printf("\n");
      printf("Usage: %s [<options>]\n", argv[0]);
      printf("\n");
      printf(" General XBraid configuration parameters\n");
      printf(" ---------------------------------------\n");
      printf("  -nt  <n>                           : number of time steps (default: 100)\n"); 
      printf("  -tol <tol>                         : set the stopping tolerance (default: 1e-6)\n");
      printf("  -lbalance <n>                      : set the type of load balance (default: 0)\n");
      printf("                                     : 0 - block data distribution \n");
      printf("                                     : 1 - weighted load balancing \n");
      printf("  -comm_min                          : use communication minimization \n");
      printf("  -no_output                         : do not save the solution in output files\n");
      printf("  -help                              : print this help and exit\n");
      printf("\n");
   }

   if( print_usage ){
      MPI_Finalize();
      return (0);
   }


   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->rank)   = rank;
   (app->lbalance)   = lbalance;
   (app->tstart)  = tstart;
   /* initialize XBraid and set options */
   braid_Init(MPI_COMM_WORLD, MPI_COMM_WORLD, tstart, tstop, ntime, app,
             my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
   
   /* Set some typical Braid parameters */
   braid_SetPrintLevel( core, 1);
   braid_SetMaxLevels(core, 15);
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, 2);
   braid_SetLoadBalance(core, lbalance);
   braid_SetCommMin(core, comm_min);
   braid_SetStorage(core, 0);
   braid_SetGlobalAccess(core, my_GlobalAccess);

   if (!output)
      braid_SetAccessLevel(core, 0);
   
   /* Run simulation, and then clean up */
   braid_Drive(core);

   braid_Destroy(core);
   free(app);
   MPI_Finalize();

   return (0);
}
