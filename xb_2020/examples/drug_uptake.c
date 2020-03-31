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
   double    tolerance;
   int       ntime;
   int       richardson;
   int       BE;
   int       adaptive;
   int       rfactor;
   char      filename[20];
   int       VCycles;
} my_App;

/* can put anything in my vector and name it anything as well */
typedef struct _braid_Vector_struct
{
   double y0,y1;

} my_Vector;

   
   int 
my_implicit_solve(double y0,
                  double y1,
                  double *k0,
                  double *k1,
                  double dt,
                  double tstop)
{
   double a = 6.0; //app->a;
   double b = 0.6; //app->b;
   double uptake = ( fmod(24*tstop,24) <= 1 ) ? 2 : 0;
   
   *k0 = ( -a*y0 + uptake ) / ( 1.0 + a*dt );
   *k1 = (  a*( y0 + dt*(*k0) ) - b*y1 ) / ( 1 + b*dt ) ; 
   return 0;
}

int 
take_step( braid_App app,
           braid_Vector u, 
           braid_Real tstart, 
           braid_Real tstop )
{
   
   double k0, k1; 
   double dt = tstop - tstart;
   double gamma = (2. - sqrt(2.))/2.;// ( L stable 2nd order )
   
   if ( app->BE ) //Backward Euler
   {
      my_implicit_solve( u->y0, u->y1 , &k0, &k1, tstop-tstart, tstop  );
      u->y0 = u->y0 + dt*k0;
      u->y1 = u->y1 + dt*k1;
   }
   else //SDIRK2 
   {
      //double gamma = (3. + sqrt(3.))/6.; //A stable 3rd order 
      my_implicit_solve(u->y0, u->y1, &k0, &k1, gamma*dt, tstart + gamma*dt );
      
      double temp_y0 = u->y0 + (1.0-2.0*gamma)*dt*k0;
      double temp_y1 = u->y1 + (1.0-2.0*gamma)*dt*k1;

      u->y0 = u->y0 + k0*dt/2.0;
      u->y1 = u->y1 + k1*dt/2.0;

      my_implicit_solve(temp_y0, temp_y1, &k0, &k1, gamma*dt, tstart + (1.0 - gamma)*dt);
      u->y0 = u->y0 + k0*(dt/2.0);
      u->y1 = u->y1 + k1*(dt/2.0);
   }

   return 0;

}
int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v;

   v = (my_Vector *) malloc(sizeof(my_Vector));
   (v->y0) = (u->y0);
   (v->y1) = (u->y1);
   *v_ptr = v;

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
   take_step(app, u, tstart, tstop);

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
      (u->y0) =  0.0;
      (u->y1) =  0.0;
   }
   else
   {
      (u->y0) = 0.456;
      (u->y1) = 0.456;
   }
   *u_ptr = u;

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
   (y->y0) = alpha*(x->y0) + beta*(y->y0);
   (y->y1) = alpha*(x->y1) + beta*(y->y1);

   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   double dot;

   dot = (u->y1)*(u->y1);
   *norm_ptr = sqrt(dot);

   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int        level,caller, iter, nrefine;
   double     t;
  
   braid_AccessStatusGetLevel(astatus, &level);
   braid_AccessStatusGetCallingFunction(astatus, &caller);
   braid_AccessStatusGetIter(astatus, &iter);
   braid_AccessStatusGetNRefine(astatus, &nrefine );
   braid_AccessStatusGetT(astatus, &t);
   
   if (level > 0 || caller != 3 )
   {
      return 0;
   }
    
   FILE *fp;   
   int file_num = iter+nrefine;
   char filename[20];
   sprintf(filename, "%.2d_", file_num);   
   strcat(filename, app->filename);
   fp = fopen(filename, "a"); 
   
   fprintf(fp," %.13g , %.13g , %.13g \n " , t, u->y0, u->y1);
   fclose(fp); 
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
   dbuffer[0] = (u->y0);
   dbuffer[1] = (u->y1);
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
   (u->y0) = dbuffer[0];
   (u->y1) = dbuffer[1];
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
   
   if ( ( iter + nrefine == 0 )  )
   {
      braid_GlobalAccessStatusSetRefine(status,0);
      return 0;
   }
  
   int num_refine = 0;
   if (estimate) {
   for ( i = 0 ; i < npoints; i++ )
   {
   
      if ( estimate[i] > app->tolerance)
      {
              double error = sqrt( estimate[i] / app->tolerance);            
              rfactors[i] = (int) ceil( error );
              if (rfactors[i] > 1 )
                 num_refine++;
      }
   }
   }
   if ( ( ( iter + nrefine ) % app->VCycles ) != 0 ) 
   {
      int glob_refine = 0;
      MPI_Allreduce(&num_refine, &glob_refine, 1, MPI_INT, MPI_SUM, comm );
      if ( num_refine == 0 )
         braid_GlobalAccessStatusSetRefine(status,-1); //Ask xbraid to try finish
      else
         braid_GlobalAccessStatusSetRefine(status,0); //Force XBraid to keep iterating
   }
   else
      braid_GlobalAccessStatusSetRefine(status, 1); //Ask XBraid to refine

   return 0 ;
}

int 
take_adaptive_step( braid_App app, 
                    braid_Vector u , 
                    braid_Real tstart, 
                    braid_Real tstop,
                    braid_Real *new_dt)
{
   double tol = app->tolerance;
   //Solutions from the large step with size dt;
   braid_Vector large_u, small_u;  
   my_Clone(app, u, &small_u);
   my_Clone(app, u, &large_u); 
   take_step(app, large_u, tstart, tstop );
   take_step(app, small_u, tstart, tstart + (tstop-tstart)/2.0);
   double a0 = small_u->y0;
   double a1 = small_u->y1;
   take_step(app, small_u, tstart + (tstop-tstart)/2.0, tstop); 
   double error = sqrt( tol / fabs(large_u->y1 - small_u->y1) );
   
   if ( error <= 1 )
   {
      *new_dt = -1;
   }
   else
   {     
      double bb = ( 0.2 > error ) ? 0.2 : error ;
      bb = ( 1.5 < error ) ? 1.5 : error ;
      bb *= 0.9*(tstop-tstart);
      bb = ( bb > 1.0/24.0 ) ? 1.0/24.0 : bb ;
      *new_dt = bb;   
      if ( app->richardson )
      {
         u->y0 = 2.0*small_u->y0 - large_u->y0;
         u->y1 = 2.0*small_u->y1 - large_u->y1;
      }
      else
      {
         u->y0 = small_u->y0;
         u->y1 = small_u->y1;
      }   
      
      FILE *fp;
      fp = fopen(app->filename, "a"); 
      fprintf(fp, " %.13g, %.13g , %.13g \n", tstart + (tstop-tstart)/2.0 , a0, a1 );
      fprintf(fp, " %.13g, %.13g , %.13g \n", tstop, u->y0, u->y1 );
      fclose(fp);
   }

   my_Free( app, small_u );
   my_Free( app, large_u );
   return 0;
}

int
adaptive_solver(braid_App app)
{
    //FILE *fp;
    //fp = fopen(app->filename, "a"); 
    int ntime = app->ntime;
    double tstart = app->tstart;
    double tstop = app->tstop;
    double dt = (tstop - tstart) / ( (double) ntime ); 
    double new_dt = 0; 
    double new_tstop;
    braid_Vector u;
    my_Init(app, tstart, &u);
    //fprintf(fp, " %.13g, %.13g , %.13g \n", tstart, u->y0, u->y1 );
    //fclose(fp);
    do 
    {
         new_tstop = ( tstart + dt < tstop ) ? tstart + dt :tstop ;
         take_adaptive_step( app, u, tstart, new_tstop , &new_dt );
         if ( new_dt < 0 )
         {
            dt /= 2.0;
         }
         else
         {
            tstart += dt;

            dt = new_dt;
         }   
    } while ( tstop - tstart > 1e-12 );

    my_Free( app, u ); 
    return 0;
}

int sequential_solver(braid_App app)
{
    //FILE *fp;
    //fp=fopen( app->filename , "a"); 
    int ntime = app->ntime;
    double tstart = app->tstart;
    double tstop = app->tstop;
    double dt = (tstop - tstart) / ( (double) ntime ); 
    braid_Vector u;
    my_Init(app, tstart, &u);
    //fprintf(fp, " %.10g, %.10g , %.10g \n", tstart, u->y0, u->y1 );
    do
    {
         take_step( app, u, tstart, tstart + dt );
         tstart += dt;
         //fprintf(fp," %.10g, %.10g , %.10g \n", tstart, u->y0, u->y1 );
    } while (tstop-tstart > 1e-9 );
    //fclose(fp);
    return 0;
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
   
   char          filename[20] = "default_out";

   double        tstart = 0.0;
   double        tstop = 10.0;
   double        tol = 1e-10;
   int           ntime = 5000; 
   int           max_levels = 100;
   int           nrelax     = 1;
   int           min_coarse = 300;
   
   double        adapt_tol     = 1.0e-10;
   int           cfactor       = 4;
   int           richardson = 0;
   int           arg_index;
   int           adaptive = 0;
   int           sequential = 0;
   int           rfactor = 2;
   int           solver = 1;
   int           Vcy = 1;
   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   comm   = MPI_COMM_WORLD;


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
            printf("  -atol <tol>       : set adaptive tolerance\n");
            printf("  -cf  <cfactor>    : set coarsening factor\n");
            printf("  -richardson       : set richardson extrap on\n");
            printf("  -adapt            : set time adaptation\n");
            printf("  -sequential       : set sequential \n");
            printf("  -rfactor          : refinement factor \n");
            printf("\n");
         }
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-filename") == 0 )
      {
         arg_index++;
         strcpy(filename, argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-Vcycs") == 0 )
      {
         arg_index++;
         Vcy = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-atol") == 0 )
      {
         arg_index++;
         adapt_tol = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-cf") == 0 )
      {
         arg_index++;
         cfactor = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-solver") == 0 )
      {
         arg_index++;
         solver = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-richardson") == 0 )
      {
         arg_index++;
         richardson = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-adapt") == 0 )
      {
         arg_index++;
         adaptive = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nt") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-seq") == 0 )
      {
         arg_index++;
         sequential = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-rf") == 0 )
      {
         arg_index++;
         rfactor = atoi(argv[arg_index++]);
      }
      else
      {
         printf( "bad parameter %s\n", argv[arg_index] );
         arg_index++;
      }
   }

   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->comm)   = comm;
   (app->tstart) = tstart;
   (app->tstop)  = tstop;
   (app->ntime)  = ntime;
   (app->adaptive) = adaptive;
   (app->richardson) = richardson;
   (app->rfactor) = rfactor; 
   (app->tolerance) = adapt_tol;
   (app->BE) = solver; 
   (app->VCycles) = Vcy;
   strcpy(app->filename, filename);



   if (sequential)
   {
      if (adaptive)
         adaptive_solver(app);
      else
         sequential_solver(app);
   }
   else
   {
      braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
               my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm,
               my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

      braid_SetPrintLevel( core, 1);
      braid_SetMaxLevels(core, max_levels);
      braid_SetNRelax(core, -1, nrelax);
      braid_SetAccessLevel(core,2);
      braid_SetAbsTol(core, tol);
      braid_SetCFactor(core, -1, cfactor);
      braid_SetMaxIter(core, 1000);
      braid_SetMinCoarse(core, min_coarse);
      braid_SetStorage(core, 0);
      braid_SetRefine( core, adaptive );
      braid_SetErrorEstimation(core, adaptive, richardson, 2);
      braid_SetGlobalAccess(core, my_GlobalAccess);
      braid_Drive(core);
      braid_Destroy(core);
   }

   free( app );
   MPI_Finalize();

   return (0);
}
