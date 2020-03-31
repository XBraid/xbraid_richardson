
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


/*
   Example 01_a

   Compile with: make ex-03

   Sample run:   mpirun -np 2 ex-03

   Description:

   Simple linear example. Solution is u = 3.0*t

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "braid.h"

/*--------------------------------------------------------------------------
 * My integration routines
 *--------------------------------------------------------------------------*/

/* can put anything in my app and name it anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm  comm;
   double    tstart;
   double    tstop;
   int       ntime;
   int       factor;
   double    a; // van der pol constant 
} my_App;

/* can put anything in my vector and name it anything as well */
typedef struct _braid_Vector_struct
{
   double x, dx ;

} my_Vector;

/* Implicit solve K = F( X + dt K ) 
 * Did it this way to use same code for SDIRK 1 -- 3 ( same as mfem ) 
 */
int 
my_implicit_solve(double x,
                  double y,
                  double *k0,
                  double *k1,
                  double dt,
                  double a)
{

   double J00,J01,J10,J11,f0,f1,det;            

   int iter = 0;
   int max_iter = 1000; //app->newton_max;
   double tol = 1e-13 ; //app->newton_tol; 
   
     
   f0 = *k0 - ( y + dt*(*k1) );
   f1 = *k1 - ( a*( y+dt*(*k1))*(1-( x+dt*(*k0))*( x+dt*(*k0))) - ( x+dt*(*k0)) );
   double resid = sqrt( f0*f0 + f1*f1 ); 

   while ( resid > tol && iter < max_iter )
   { 
      J00 = 1;
      J01 = -dt;
      J10 = -dt*( -2*a*( y + dt*(*k1) )*( x + dt*(*k0) ) - 1 );
      J11 = 1 - dt*( a*(1 - ( x + dt*(*k0) )*( x + dt*(*k0)) ) ) ;
      det = 1.0/( J00*J11 - J10*J01 );
           
      *k0 = *k0 -  det*( J11*f0 - J01*f1 );
      *k1 = *k1 -  det*( J00*f1 - J10*f0 );

      f0 = *k0 - ( y + dt*(*k1) );
      f1 = *k1 - ( a*( y+dt*(*k1))*(1-( x+dt*(*k0))*( x+dt*(*k0))) - ( x+dt*(*k0)) );
      
      resid = sqrt( f0*f0 + f1*f1 ); 
      iter++;
   }

   return 0;
}

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   double dt = tstop - tstart ; 
   double k0 = 0;//( ustop->x - u->x ) / dt  ; //Initial guess
   double k1 = 0;//( ustop->dx - u->dx ) / dt ; //Initial guess
   
   if ( 1 ) //Backward Euler
   {
      my_implicit_solve( u->x, u->dx , &k0, &k1, dt, app->a );
      u->x = u->x + dt*k0;
      u->dx = u->dx + dt*k1;
   }
   else if ( 0 )//SDIRK2 
   {
      double gamma = (3. + sqrt(3.))/6.; //A stable 3rd order 
      double gamma = (2. - sqrt(2.))/2.; ( L stable 2nd order )
      my_implicit_solve(u->x, u->dx, &k0, &k1, gamma*dt, app->a );
      
      double temp_x = u->x + (1.0-2.0*gamma)*dt*k0;
      double temp_dx = u->dx + (1.0-2.0*gamma)*dt*k1;

      u->x = u->x + k0*dt/2.0;
      u->dx = u->dx + k1*dt/2.0;

      my_implicit_solve(temp_x, temp_dx, &k0, &k1, gamma*dt, app->a);
      u->x = u->x + k0*(dt/2.0);
      u->dx = u->dx + k1*(dt/2.0);
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
   if (t == 0.0)
   {
      /* Initial condition */
      (u->x) =  2.0;
      (u->dx) =  0.0;
   }
   else
   {
      /* Initialize all other time points */
      (u->x) = rand()/RAND_MAX;
      (u->dx) = rand()/RAND_MAX;
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
   (v->x) = (u->x);
   (v->dx) = (u->dx);
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
   (y->x) = alpha*(x->x) + beta*(y->x);
   (y->dx) = alpha*(x->dx) + beta*(y->dx);

   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   double dot;

   dot = (u->x)*(u->x) + (u->dx)*(u->dx);
   *norm_ptr = sqrt(dot);

   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int        iter, nrefine,level,caller,done;
   double     t;
   braid_AccessStatusGetLevel(astatus, &level);
   braid_AccessStatusGetCallingFunction(astatus, &caller);
   braid_AccessStatusGetT(astatus, &t);
   braid_AccessStatusGetIter(astatus, &iter);
   braid_AccessStatusGetNRefine(astatus, &nrefine );
   braid_AccessStatusGetDone(astatus, &done);   
   
   if (done)
   {
      FILE *fp;
      char filename[20];
      sprintf(filename, "temp2.csv");
      fp=fopen(filename, "a");
      if(fp == NULL)
      {
         exit(-1);
      }
      fprintf(fp, " %.15f , %.15f , %.15f \n " , t, u->x, u->dx);
      fclose(fp);
   }
   
   

   return 0;
}

int
my_BufSize(braid_App          app,
           int                *size_ptr,
           braid_BufferStatus bstatus)
{
   *size_ptr = 2*sizeof(double);
   return 0;
}

int
my_BufPack(braid_App          app,
           braid_Vector       u,
           void               *buffer,
           braid_BufferStatus bstatus)
{
   double *dbuffer = buffer;

   dbuffer[0] = (u->x);
   dbuffer[1] = (u->dx);
   braid_BufferStatusSetSize( bstatus, 2*sizeof(double) );

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
   (u->x) = dbuffer[0];
   (u->dx) = dbuffer[1];
   *u_ptr = u;

   return 0;
}

int
my_GlobalAccess(braid_App    app,
                MPI_Comm     comm,
                braid_Vector *u,
                braid_Int    *umask,
                braid_Int     npoints,
                braid_GlobalAccessStatus gstatus)
{

   
   braid_Int *rfactors,i, npoints1;
   braid_Real *estimate, temp_est;
   braid_GlobalAccessStatusGetRFactors(gstatus, &npoints, &rfactors);
   braid_GlobalAccessStatusGetErrorEst(gstatus, &npoints1, &estimate);
   for ( i = 0 ; i < npoints; i++ )
   {
      temp_est = -1;
      if (i < npoints1 )
         temp_est = estimate[i];
   
      if ( temp_est > 1e-6 )
      {
         rfactors[i] = 16; 
      }
    
      //if ( umask[i] == 0 )
         //printf(" %d %f %d %f\n " , i, U[i]->value , rfactors[i], temp_est);
      //else
         //printf(" %d not stored %d %f \n ", i ,rfactors[i], temp_est);
   }

   
   return _braid_error_flag;
}


/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/
#include <execinfo.h>

int main (int argc, char *argv[])
{
   braid_Core    core;
   my_App       *app;
   MPI_Comm      comm;
   double        tstart, tstop;
   int           ntime;
   int           tpts       = 10000000;
   int           max_levels = 3;
   int           nrelax     = 1;
   int           nrelax0    = -1;
   double        tol        = 1.0e-8;
   int           cfactor    = 4;
   int           refine_time = 0;
   int           arg_index;
   int           min_coarse = 2000;
   int           factor = 2;
   double        a = 1000;
   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   /* ntime time intervals with spacing 1 */
   comm   = MPI_COMM_WORLD;
   ntime  = 10000;
   tstart = 0.0;
   tstop  = 30;
   /* Parse command line */

   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         int  myid;
         MPI_Comm_rank(comm, &myid);
         if ( myid == 0 )
         {
            printf("\n");
            printf("  -ml  <max_levels> : set max levels\n");
            printf("  -nu  <nrelax>     : set num F-C relaxations\n");
            printf("  -nu0 <nrelax>     : set num F-C relaxations on level 0\n");
            printf("  -tol <tol>        : set stopping tolerance\n");
            printf("  -cf  <cfactor>    : set coarsening factor\n");
            printf("  -mi  <max_iter>   : set max iterations\n");
            printf("  -lb               : use load balencing\n");
            printf("  -rt               : use temporal refinment\n");
            printf("  -tpts <tpts>      : cutoff for time refinement\n");
            printf("  -nt <ntime>       : number of time points\n");
           
            printf("\n");
         }
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-ml") == 0 )
      {
         arg_index++;
         max_levels = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu") == 0 )
      {
         arg_index++;
         nrelax = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu0") == 0 )
      {
         arg_index++;
         nrelax0 = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tol") == 0 )
      {
         arg_index++;
         tol = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-cf") == 0 )
      {
         arg_index++;
         cfactor = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tpts") == 0 )
      {
         arg_index++;
         tpts = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nt") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mi") == 0 )
      {
         arg_index++;
         factor = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-lb") == 0 )
      {
         arg_index++;
      }
      else if ( strcmp(argv[arg_index], "-rt") == 0 )
      {
         arg_index++;
         refine_time = 1;
      }
      else
      {
         arg_index++;
         /*break;*/
      }
   }

   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->comm)   = comm;
   (app->tstart) = tstart;
   (app->tstop)  = tstop;
   (app->ntime)  = ntime;
   (app->factor) = factor;
   (app->a)      = a;

   braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
              my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm,
              my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

   braid_SetPrintLevel( core, 1);
   braid_SetMaxLevels(core, max_levels);
   braid_SetNRelax(core, -1, nrelax);
   braid_SetAccessLevel(core,2);
   if (nrelax0 > -1)
   {
      braid_SetNRelax(core,  0, nrelax0);
   }
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   /*braid_SetCFactor(core,  0, 10);*/
   braid_SetMaxIter(core, 1000);
   braid_SetMinCoarse(core, min_coarse);
   braid_SetStorage(core, 0);
   braid_SetRefine( core, 1 );
   braid_SetTPointsCutoff( core, tpts);
   braid_SetGlobalAccess(core, my_GlobalAccess);
   braid_SetErrorEstimation(core,1,1,2); 

   int myid;

   MPI_Comm_rank( comm, &myid );
   braid_Drive(core);

   braid_Destroy(core);
   free( app );
   /* Finalize MPI */
   MPI_Finalize();

   return (0);
}
