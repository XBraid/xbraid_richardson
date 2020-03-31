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
 * Example:       ex-02.c
 *
 * Interface:     C
 *
 * Requires:      only C-language support
 *
 * Compile with:  make ex-02
 *
 * Help with:     ex-02 -help
 *
 * Sample run:    mpirun -np 2 ex-02
 *
 * Description:   Solves the 1D heat equation, with only parallelism in time
 *                Space-time domain:  [0, PI] x [0, 2*PI]
 *                Exact solution:     u(t,x) = sin(x)*cos(t)
 *
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "braid.h"
#include "braid_test.h"
#include "ex-02-lib.c"

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

/* can put anything in my app and name it anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm  comm;
   double    tstart;       /* Define the temporal domain */
   double    tstop;
   int       ntime;
   double    xstart;       /* Define the spatial domain */
   double    xstop;
   int       nspace;
   double    matrix[3];    /* the three point spatial discretization stencil */
   double *  g;            /* temporary vector for inversions and mat-vecs */
   char      filename[40];
   int       richardson;
   int       cfactor;
   int       sequential;
   int       adaptive;
   int       VCycles;
   double    AdaptTol;

} my_App;

/* Can put anything in my vector and name it anything as well */
typedef struct _braid_Vector_struct
{
   int     size;
   double *values;

} my_Vector;

/* create and allocate a vector */
void
create_vector(my_Vector **u,
              int size)
{
   (*u) = (my_Vector *) malloc(sizeof(my_Vector));
   ((*u)->size)   = size;
   ((*u)->values) = (double *) malloc(size*sizeof(double));
}


/*--------------------------------------------------------------------------
 * My integration routines
 *--------------------------------------------------------------------------*/



int my_Step(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     fstop,
            braid_Vector     u,
            braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   int level;
   double deltaX, deltaT;

   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   deltaT = tstop - tstart;
   deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);

   /* Take backward Euler step
    * Note: if an iterative solver were used, ustop->values would
    *       contain the XBraid's best initial guess. */
   take_step(u->values, u->size, tstop, app->xstart, deltaX, deltaT,
             app->matrix, app->g);

   return 0;
}


int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;
   int    i;
   int nspace = (app->nspace);
   double deltaX = (app->xstop - app->xstart) / (nspace - 1.0);

   /* Allocate vector */
   create_vector(&u, nspace);

   /* Initialize vector (with correct boundary conditions) */
   if(t == 0.0)
   {
      /* Get the solution at time t=0 */
      get_solution(u->values, u->size, 0.0, app->xstart, deltaX);
   }
   else
   {
      /* Use random values for u(t>0), this measures asymptotic convergence rate */
      srand(0);
      for(i=0; i < nspace; i++)
      {
         (u->values)[i] = ((double)rand())/RAND_MAX;
      }
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
   int size = (u->size);
   int i;

   create_vector(&v, size);
   for (i = 0; i < size; i++)
   {
      (v->values)[i] = (u->values)[i];
   }
   *v_ptr = v;

   return 0;
}

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u->values);
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
   int i;
   int size = (y->size);

   for (i = 0; i < size; i++)
   {
      (y->values)[i] = alpha*(x->values)[i] + beta*(y->values)[i];
   }

   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   int    i;
   int size   = (u->size);
   double dot = 0.0;
   double deltaX = (app->xstop - app->xstart) / (u->size - 1.0);

   for (i = 0; i < size; i++)
   {
      dot += (u->values)[i]*(u->values)[i];
   }
   *norm_ptr = sqrt(dot*(deltaX));

   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int        level, done,index;
   double     t, error;

   braid_AccessStatusGetT(astatus, &t);
   braid_AccessStatusGetLevel(astatus, &level);
   braid_AccessStatusGetDone(astatus, &done);
   braid_AccessStatusGetTIndex(astatus, &index);
   /* IF on the finest level AND print_level is high enough AND at the final time,
    * THEN print out the discretization error */
   FILE *fp;
   fp = fopen(app->filename,"a");
   if( (level == 0) && done )
   {
      error = compute_error_norm(u->values, app->xstart, app->xstop, u->size, t);
      fprintf(fp, "  %1.10e,   %1.10e\n", t, error );
   }
   fclose(fp);
   return 0;
}

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
   int size = (app->nspace);
   *size_ptr = (size+1)*sizeof(double);
   return 0;
}

int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus)
{
   double *dbuffer = buffer;
   int i, size = (u->size);

   dbuffer[0] = size;
   for (i = 0; i < size; i++)
   {
      dbuffer[i+1] = (u->values)[i];
   }

   braid_BufferStatusSetSize( bstatus,  (size+1)*sizeof(double));

   return 0;
}

int
my_BufUnpack(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus)
{
   my_Vector *u = NULL;
   double    *dbuffer = buffer;
   int        i, size;

   size = dbuffer[0];
   create_vector(&u, size);

   for (i = 0; i < size; i++)
   {
      (u->values)[i] = dbuffer[i+1];
   }
   *u_ptr = u;

   return 0;
}

int my_GlobalAccess(braid_App                app,
                    MPI_Comm                 comm,
                    braid_Vector            *u_ptr,
                    braid_Int               *umask,
                    braid_Int                npoints,
                    braid_GlobalAccessStatus status )
{

   braid_Int *rfactors,i, npoints1, iter, nrefine;
   braid_Real *estimate, *tgrid;
   braid_GlobalAccessStatusGetRFactors(status, &npoints, &rfactors);
   braid_GlobalAccessStatusGetIter(status,&iter);
   braid_GlobalAccessStatusGetNRefine(status, &nrefine);
   braid_GlobalAccessStatusGetErrorEst(status, &npoints1, &estimate);
   braid_GlobalAccessStatusGetTGrid(status, &npoints, &tgrid );
   
   if ( !app->adaptive ) 
   {   
      braid_GlobalAccessStatusSetRefine(status,-1); //No refine, and allow to quit 
   }
   else if ( iter + nrefine == 0  || ( app->VCycles > 0  && (  ( iter + nrefine) % app->VCycles != 0 ) ) )
   {
      braid_GlobalAccessStatusSetRefine(status,0); //No refine + force another iteration 
   }           
   else
   {
      for (i = 0 ; i < npoints; i++ )
      {
         double error = sqrt( estimate[i] / app->AdaptTol);
         rfactors[i] = (int) ceil( error );
      }
      braid_GlobalAccessStatusSetRefine(status,1); //Refine in time (auto forces another iteration)
   }
   return 0;
}


int sequential_solver(braid_App    app)
{
   int    ntime = app->ntime;
   double tstart = app->tstart;
   double tstop = app->tstop;

   braid_Vector u_init;
   my_Init(app, tstart, &u_init);
   braid_Vector u_small, u_large;
   FILE *fp;
   fp = fopen(app->filename, "a"); 


   int i,j;
   double deltaX, deltaT;
   deltaT = (tstop - tstart) / ( (double) ntime );
   deltaX = (app->xstop - app->xstart) / (u_init->size - 1.0);

   if ( !(app->adaptive) && !(app->richardson) )
   {
      for ( i = 0; i < ntime ; i++ )
      {
         tstart += deltaT;
         take_step( u_init->values, u_init->size, tstart, app->xstart, deltaX, deltaT, app->matrix, app->g);
         double error = compute_error_norm(u_init->values, app->xstart, app->xstop, u_init->size, tstart);
         fprintf(fp, "%1.4e,%1.4e \n" , tstart, error);
         
      }
      my_Free(app, u_init );
   }
   else
   {
      double cfactor =  (double) app->cfactor;
      my_Init(app,tstart,&u_small);
      my_Init(app,tstart,&u_large);

      while ( fabs( tstart - tstop ) > 1e-10  && tstart < tstop  )
      {
         take_step( u_large->values, u_large->size, tstart+deltaT, app->xstart, deltaX, deltaT, app->matrix, app->g);
         take_step( u_small->values, u_small->size, tstart+deltaT/2.0, app->xstart, deltaX, deltaT/2.0, app->matrix, app->g); 
         double aerror = compute_error_norm(u_small->values, app->xstart, app->xstop, u_small->size, tstart+deltaT/2.0);
         take_step( u_small->values, u_small->size, tstart+deltaT, app->xstart, deltaX, deltaT/2.0, app->matrix, app->g);
         double berror = compute_error_norm(u_small->values, app->xstart, app->xstop, u_small->size, tstart+deltaT);
         int done = 0;
         if  (app->adaptive)
         {
            double error = 0;
            for ( j = 0; j < u_init->size; j++ )
            {
               error += pow( (u_small->values)[j] - (u_large->values)[j] , 2.0 );
            }
            error = sqrt(error*deltaX);
            double rerror = sqrt( (app->AdaptTol) / error );
            
            if ( rerror <= 1 ) // repeat the time step
            {
               done = 1;
               deltaT *= 0.5;
               for ( j = 0; j < u_init->size; j++ )
               {
                  (u_small->values)[j] = (u_init->values)[j];
                  (u_large->values)[j] = (u_init->values)[j];
               }
            }
            else
            {
               fprintf(fp, "%1.4e,%1.4e \n" , tstart + deltaT/2.0, aerror);
               fprintf(fp, "%1.4e,%1.4e \n" , tstart + deltaT, berror);
               tstart += deltaT;
               double bb = ( 0.2 > rerror ) ? 0.2 : rerror ;
               bb = ( 1.5 < rerror ) ? 1.5 : rerror ;
               deltaT = bb*0.9*(deltaT);
               if ( tstart + deltaT > tstop )
                  deltaT = tstop - tstart ; 
            }
         }

         if ( !done  )
         {
            double order =  1.0;
            double a = pow ( cfactor , order ) / ( pow ( cfactor, order ) - 1 );
            double b = 1.0 / ( pow( cfactor, order ) - 1 );
            for ( j = 0; j < u_init->size  ; j++ )
            {
               if (app->richardson)
               {
                  u_small->values[j] = a*u_small->values[j] - b*u_large->values[j] ;
               }
               u_large->values[j] = u_small->values[j];
               u_init->values[j] = u_small->values[j];
            }
         }
      }

      double error = compute_error_norm(u_small->values, app->xstart, app->xstop, u_small->size, tstop);
      my_Free(app, u_init );
      my_Free(app, u_small);
      my_Free(app, u_large);
     printf("  Discretization error at final time:  %1.4e , %1.4e\n", tstart, error);
   }
   fclose(fp);
   return 0;
}

int
adaptive_sequential_solver( braid_App app)
{
   return 0;
}


/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   braid_Core    core;
   my_App       *app;
   MPI_Comm      comm;
   int           rank, arg_index;

   /* Define space-time domain */
   double    tstart        =  0.0;
   double    tstop         =  1.0;
   int       ntime         =  64;
   double    xstart        =  -1.0;
   double    xstop         =  1.0;
   int       nspace        =  33;

   /* Define XBraid parameters
    * See -help message for descriptions */
   int       max_levels    = 200;
   int       nrelax        = 1;
   int       skip          = 0;
   double    tol           = 1.0e-07;
   int       max_iter      = 30;
   int       min_coarse    = 3;
   int       print_level   = 1;
   int       access_level  = 1;
   int       use_sequential= 0;

   int       sequential = 0;
   int       richardson = 0;
   int       adaptive   = 0;
   int       Vcycles    = 0;
   int       cfactor    = 2;
   double    AdaptTol   = 1e-6;
   char      filename[40] = "Def_output"; 
   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   comm   = MPI_COMM_WORLD;
   MPI_Comm_rank(comm, &rank);

   /* Parse command line */
   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         if ( rank == 0 )
         {
            printf("\n");
            printf(" Solve the 1D heat equation on space-time domain:  [0, PI] x [0, 2*PI]\n");
            printf(" with exact solution u(t,x) = sin(x)*cos(t) \n\n");
            printf("   -ntime <ntime>       : set num points in time\n");
            printf("   -nspace <nspace>     : set num points in space\n\n");
            printf("   -tol  <tol>          : set stopping tolerance (scaled later by sqrt(dt) sqrt(dx))\n");
            printf("   -cf   <cfactor>      : set coarsening factor on all levels \n");
         
	    printf("Eg. ./ex-02 --> Runs standard mgrit\n");
	    printf("Eg. ./ex-02 -adaptive 1 --> Runs mgrit with automatic refinement. \n");


	 
	 }
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-sequential") == 0 )
      {
         arg_index++;
         sequential = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-fname") == 0 )
      {
         arg_index++;
         strcpy(filename, argv[arg_index++]);

      }

      else if ( strcmp(argv[arg_index], "-richardson") == 0 )
      {
         arg_index++;
         richardson = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-adaptive") == 0 )
      {
         arg_index++;
         adaptive = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-adapt_tol") == 0 )
      {
         arg_index++;
         AdaptTol = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-Vcycls") == 0 )
      {
         arg_index++;
         Vcycles = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nspace") == 0 )
      {
         arg_index++;
         nspace = atoi(argv[arg_index++]);
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
      else
      {
         printf("ABORTING: incorrect command line parameter %s\n", argv[arg_index]);
         MPI_Finalize();
         return (0);

      }
   }

   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->g)             = (double*) malloc( nspace*sizeof(double) );
   (app->comm)          = comm;
   (app->tstart)        = tstart;
   (app->tstop)         = tstop;
   (app->ntime)         = ntime;
   (app->xstart)        = xstart;
   (app->xstop)         = xstop;
   (app->nspace)        = nspace;
   (app->adaptive)      = adaptive;
   (app->richardson)    = richardson;
   (app->VCycles)       = Vcycles;
   (app->AdaptTol)      = AdaptTol;
   (app->cfactor)       = cfactor;
   strcpy(app->filename, filename);

   /* Initialize storage for sc_info, for tracking space-time grids visited during the simulation */
   if ( sequential  )
   {
      double localtime = MPI_Wtime();
      sequential_solver(app);

      double finaltime = MPI_Wtime();
      printf(" Wall time = %1.4f seconds \n" , finaltime-localtime );
      return 0;
   }

   /* Initialize Braid */
   braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
              my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm,
              my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

   /* The first step before running simulations, is always to verify the wrapper tests */
   braid_SetPrintLevel( core, print_level);
   braid_SetAccessLevel( core, access_level);
   braid_SetMaxLevels(core, max_levels);
   braid_SetMinCoarse( core, min_coarse );
   braid_SetSkip(core, skip);
   braid_SetNRelax(core, -1, nrelax);
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetMaxIter(core, max_iter);
   braid_SetSeqSoln(core, use_sequential);
   braid_SetStorage(core, 0);
   braid_SetRefine(core,1);
   braid_SetErrorEstimation(core, adaptive, richardson, 2);
   braid_SetGlobalAccess(core, my_GlobalAccess);
   braid_Drive(core);

   /* Clean up */
   braid_Destroy(core);
   free( app->g);
   free( app );

   /* Finalize MPI */
   MPI_Finalize();

   return (0);
}
