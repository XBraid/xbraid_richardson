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

/** \file braid_dist.c
 * \brief Source code for data distribution routines.  See braid_dist.h for more information.
 *
 */
#include "_braid_dist.h"
#include "_braid.h"
#include "braid_defs.h"
/*--------------------------------------------------------------------------------
 * Build the Balance Structure. 
 * *-------------------------------------------------------------------------------*/

braid_Int
_braid_BuildBalanceStructure(braid_Core            core,
                             braid_Int            *done,
                             _braid_BalanceStruct *bstruct)
{
   _braid_Grid **grids         = _braid_CoreElt(core, grids);
   braid_Int   lbalance        = _braid_CoreElt(core, lbalance );
   braid_Int   refine          = _braid_CoreElt( core, refine );
   braid_Int   gupper          = _braid_GridElt( grids[0], gupper);
   braid_Int   ilower          = _braid_GridElt( grids[0], ilower);
   braid_Int   iupper          = _braid_GridElt( grids[0], iupper);
 
   braid_Int finished = 0; /* Done is set if no load bal and no refine */ 
                
   /* Build Phase 1: Initialization */
   if ( !finished )
      _braid_InitBalanceStruct(bstruct, refine, lbalance, ilower, iupper, gupper);
   
   /* Build Phase 2: Get the refined interval */ 
   if ( !finished )
      _braid_GetRefinedInterval(core, &finished, bstruct);

   /* Build Phase 3: Get the new interval */
   if ( !finished )
      _braid_GetFineInterval(core,&finished,bstruct);

   /* Build Phase 4: Build the communication map */   
   if ( !finished )
      _braid_BuildCommunicationMap(core, &finished, bstruct);
  
   /* Build Phase 5: Get the mapping between the coarse grid indicies 
    * and the new refined grid */
   if ( !finished )
      _braid_GetRefinedGridMap(core, &finished, bstruct);
   
   /* Build Phase 6: Get the mapping between the coarse grid and the 
    * new balanced fine grid */
   if ( !finished )
      _braid_GetFineGridMap(core, &finished, bstruct);

   *done = finished;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------------
 * Build Phase 1: Initialization of the balance structure 
 *-------------------------------------------------------------------------------*/

braid_Int
_braid_InitBalanceStruct(_braid_BalanceStruct *bstruct,
                         braid_Int             refine,
                         braid_Int             lbalance,
                         braid_Int             coarse_ilower,
                         braid_Int             coarse_iupper,
                         braid_Int             coarse_gupper)

{
   _braid_BalanceElt( bstruct , refine )          = refine;
   _braid_BalanceElt( bstruct , lbalance )        = lbalance ;
   _braid_BalanceElt( bstruct , coarse_ilower )   = coarse_ilower;
   _braid_BalanceElt( bstruct , coarse_iupper )   = coarse_iupper;
   _braid_BalanceElt( bstruct , coarse_gupper )   = coarse_gupper;
   
   _braid_BalanceElt( bstruct , refined_ilower )  = -1;
   _braid_BalanceElt( bstruct , refined_iupper )  = -1;
   _braid_BalanceElt( bstruct , refined_gupper )  = -1;
   _braid_BalanceElt( bstruct , refined_npoints)  = -1;
   _braid_BalanceElt( bstruct , refined_ca)       = NULL;
   _braid_BalanceElt( bstruct , refined_fa)       = NULL;
   _braid_BalanceElt( bstruct , refined_ta)       = NULL;
   _braid_BalanceElt( bstruct , refined_taalloc) = NULL;

   _braid_BalanceElt( bstruct , fine_ilower )     = -1;
   _braid_BalanceElt( bstruct , fine_iupper )     = -1;
   _braid_BalanceElt( bstruct , fine_gupper )     = -1;
   _braid_BalanceElt( bstruct , fine_first  )     = -1;
   _braid_BalanceElt( bstruct , fine_next   )     = -1;
   _braid_BalanceElt( bstruct , fine_ca)          = NULL;
   _braid_BalanceElt( bstruct , fine_ta)          = NULL;
   _braid_BalanceElt( bstruct , fine_ta_alloc)    = NULL;

   _braid_BalanceElt( bstruct , right_procs )     = NULL;
   _braid_BalanceElt( bstruct , left_procs  )     = NULL;
   _braid_BalanceElt( bstruct , send_map_alloc )  = NULL;
   _braid_BalanceElt( bstruct , send_map )        = NULL;
   _braid_BalanceElt( bstruct , recv_map )        = NULL;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------------
 * Build Phase 2: Get The refined interval  
 *-------------------------------------------------------------------------------*/

braid_Int 
_braid_GetRefinedInterval( braid_Core             core,
                            braid_Int             *done,
                            _braid_BalanceStruct  *bstruct )
{

     
    MPI_Comm  comm              = _braid_CoreElt( core, comm );    
    braid_Real *wfactors        = _braid_CoreElt(core, wfactors);    
    braid_Int  *rfactors        = _braid_CoreElt(core, rfactors);
    braid_Int   nrefine         = _braid_CoreElt(core, nrefine);
    braid_Int   max_refinements = _braid_CoreElt(core, max_refinements);
    braid_Int   tpoints_cutoff  = _braid_CoreElt(core, tpoints_cutoff);
    braid_Int   iter            = _braid_CoreElt(core, niter);

    braid_Int   c_ilower      = _braid_BalanceElt(bstruct, coarse_ilower);    
    braid_Int   c_iupper      = _braid_BalanceElt(bstruct, coarse_iupper);    
    braid_Int   c_gupper      = _braid_BalanceElt(bstruct, coarse_gupper);
    braid_Int   c_npoints     = c_iupper-c_ilower + 1;
    
    braid_Real *new_wfactors;
    braid_Int myid, i, ii,j, cfactor;
    braid_Int r_gupper, r_npoints, r_iupper;
    MPI_Comm_rank(comm, &myid);
   
    /* Defualt to no refinement */
   _braid_BalanceElt(bstruct, refined_npoints) =  c_npoints;    
   _braid_BalanceElt( bstruct, refined_gupper ) = c_gupper;
   _braid_BalanceElt( bstruct, refined_ilower ) = c_ilower;
   _braid_BalanceElt( bstruct, refined_iupper ) = c_iupper; 
   _braid_GetCFactor(core,0,&cfactor);
   
   /* If we are refining in time */
   if ( _braid_BalanceElt(bstruct, refine ) )
   {
      if(   !((nrefine < max_refinements) && (c_gupper < tpoints_cutoff)) )
      {
         _braid_BalanceElt( bstruct, refine ) = 0;
         _braid_CoreElt(core, refine ) = 0;
         _braid_CoreElt(core, rstopped) = iter;
         return _braid_error_flag;
      }
      
      /* Process the refinement factors */
      r_npoints = 0;
      for (i = c_ilower; i <= c_iupper; i++)
      {
         ii = i - c_ilower;
         if (rfactors[ii] < 1)
         {
            _braid_Error(braid_ERROR_GENERIC, "Refinement factor smaller than one");
            rfactors[ii] = 1;
         }
         r_npoints += rfactors[ii];
      }
      
      /* Get the number of time steps on refined interval */
      MPI_Allreduce(&r_npoints, &r_gupper , 1, braid_MPI_INT, MPI_SUM, comm);
      (r_gupper)--;

      /* If nothing was refined, then the defaults are correct */
      if ( r_gupper == c_gupper  )
      {
         _braid_BalanceElt( bstruct, refine) = 0 ;
         return _braid_error_flag;
      }
      
      /* Get the new ilower and iupper on the refined grid */
      MPI_Scan(&r_npoints, &r_iupper, 1, braid_MPI_INT, MPI_SUM, comm);
      _braid_BalanceElt(bstruct, refined_npoints) =  r_npoints;    
      _braid_BalanceElt( bstruct, refined_gupper ) = r_gupper;
      _braid_BalanceElt( bstruct, refined_ilower ) = r_iupper - r_npoints;
      _braid_BalanceElt( bstruct, refined_iupper ) = r_iupper - 1;
      
      /* If we are refining in time, then we need to get weights for the
       * new time steps. For now, this assumes the weight of any new time
       * step is the same as its parent. No communication is needed to do
       * this
      */      
      if ( _braid_BalanceElt(bstruct, lbalance) ) 
      {
         j = 0;
         new_wfactors = _braid_CTAlloc(braid_Real, r_npoints );
         for ( i = 0; i < c_npoints; i++ )
         {
            for ( ii = 0; ii < rfactors[i]; ii++ )
            {
               new_wfactors[j++] = wfactors[i];
            }
         }
  
         _braid_TFree(wfactors);
         _braid_CoreElt(core, wfactors) = new_wfactors;
         wfactors = _braid_CoreElt(core, wfactors);
         if ( _braid_BalanceElt(bstruct, refined_ilower) == 0 )
         {
            wfactors[0] = 0; 
         }
       }
   }
   
   return _braid_error_flag;
}    

/*--------------------------------------------------------------------------------
 * Build Phase 2b: Get the mapping between the coarse grid and the refined grid 
 *-------------------------------------------------------------------------------*/
braid_Int 
_braid_GetRefinedGridMap(braid_Core            core,
                         braid_Int            *done,
                         _braid_BalanceStruct *bstruct)
{
   
   MPI_Comm     comm     = _braid_CoreElt(core, comm);
  _braid_Grid **grids    = _braid_CoreElt(core, grids);
   braid_Int *rfactors   = _braid_CoreElt(core, rfactors);
   braid_Real *ta        = _braid_GridElt(grids[0], ta);
   
   braid_Int refine      = _braid_BalanceElt(bstruct, refine );
   braid_Int ilower      = _braid_BalanceElt(bstruct, coarse_ilower);
   braid_Int iupper      = _braid_BalanceElt(bstruct, coarse_iupper);
   braid_Int gupper      = _braid_BalanceElt(bstruct, coarse_gupper);
   braid_Int r_ilower    = _braid_BalanceElt(bstruct, refined_ilower);
   braid_Int r_iupper    = _braid_BalanceElt(bstruct, refined_iupper);
   braid_Int r_gupper    = _braid_BalanceElt(bstruct, refined_gupper);
   braid_Int r_npoints   = _braid_BalanceElt(bstruct, refined_npoints); 
   
   braid_Real  *r_ta_alloc, *r_ta; 
   braid_Int *r_ca, *r_fa;
   braid_Int r_ii, ii, i, j, rfactor, ncomms, npoints, cfactor, prevproc;
   MPI_Request *requests;
   MPI_Status  *statuses;
  
   npoints   = iupper - ilower + 1;
   r_ta_alloc = _braid_CTAlloc(braid_Real, r_npoints+2);
   r_ca      = _braid_CTAlloc(braid_Int,  r_npoints);
   r_fa      = _braid_CTAlloc(braid_Int,  npoints+1);
   r_ta      =  &r_ta_alloc[1];

   if (refine)
   {
   
      r_ta[-1]=ta[-1];
      r_ii = 0;
      
      for (i = (ilower-1); i < iupper; i++)
      {
         ii = i-ilower;
         rfactor = rfactors[ii+1];

         for (j = 1; j <= rfactor; j++)
         {
            if (j < rfactor)
            {
               r_ca[r_ii] = -1;
               /* This works because we have ta[-1] */
               r_ta[r_ii] = ta[ii] + (((braid_Real)j)/rfactor)*(ta[ii+1]-ta[ii]);
            }
            else
            {
               r_ca[r_ii] = i+1;
               r_ta[r_ii] = ta[ii+1];
               r_fa[ii+1] = r_ilower + r_ii;
            }

            r_ii++;
         }
      }
   }
   else
   {
      for ( i = r_ilower; i <= r_iupper; i++ )
      {
         r_ca[ i - r_ilower ]  = i;
         r_ta[i-r_ilower] = ta[i-r_ilower];
         r_fa[ i - r_ilower ] = i;
      }
   }

   /* Get the next r_fa value to my right */
   ncomms = 2; /* Upper bound */
   requests = _braid_CTAlloc(MPI_Request, ncomms);
   statuses = _braid_CTAlloc(MPI_Status,  ncomms);
   ncomms = 0;
   if (npoints > 0)
   {
      braid_Real send_buf[2], recv_buf[2];
      send_buf[0]=r_fa[0];
      send_buf[1]=r_ta[0];     
      /* Post r_fa receive */
      r_fa[npoints] = r_gupper+1;
      if (iupper < gupper)
      {
         MPI_Irecv(recv_buf, 2, braid_MPI_REAL, MPI_ANY_SOURCE, 2, comm,
                   &requests[ncomms++]);
      }

      /* Post r_fa send */
      if (ilower > 0)
      {
         _braid_GetProcLeftOrRight( core, 0, -1 , &prevproc);
         MPI_Isend(send_buf, 2, braid_MPI_REAL, prevproc, 2, comm,
                   &requests[ncomms++]);
      }
      MPI_Waitall(ncomms, requests, statuses);
      if ( iupper < gupper )
      {
         r_fa[npoints]=recv_buf[0];
         r_ta[r_npoints]=recv_buf[1];
      }
   }
   _braid_TFree(requests);
   _braid_TFree(statuses);

   /* If storing only C-points on the fine grid (and NOT using shell vectors),
    *  modify r_ca to mark only those coarse points that need to be sent to
    * initialize the C-points */

   _braid_GetCFactor(core, 0, &cfactor);
   if (_braid_CoreElt(core, storage) != 0 && _braid_CoreElt(core, useshell) != 1)
   {
      for (ii = 0; ii < npoints; ii++)
      {
         /* If the index for the next coarse point is not larger than the index
          * for the next C-point, then the coarse point is not needed */
         if ( !(r_fa[ii+1] > _braid_NextCPoint(r_fa[ii], cfactor)) )
         {
            r_ii = r_fa[ii] - r_ilower;
            r_ca[r_ii] = -1;
         }
      }
   }

   /* Move control of the refined grid maps to the balance structure */
   _braid_BalanceElt( bstruct, refined_ca) = r_ca;
   _braid_BalanceElt( bstruct, refined_fa) = r_fa;
   _braid_BalanceElt( bstruct, refined_ta) = r_ta;
   _braid_BalanceElt( bstruct, refined_taalloc) = r_ta_alloc;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------------
 * Build Phase 3: Get The new fine grid interval  
 *-------------------------------------------------------------------------------*/
braid_Int 
_braid_GetFineInterval(braid_Core            core,
                       braid_Int            *done,
                       _braid_BalanceStruct *bstruct)
{
    
    braid_Int fine_ilower, fine_iupper, fine_gupper, coarse_gupper;
    
    if ( _braid_BalanceElt(bstruct, lbalance) )
    {
       /* If user defined weighted load balancing is turned on, use the 
        * weighted fine interval to calculate new interval */
       _braid_GetWeightedFineInterval(core, done, bstruct);
    }
    else if ( _braid_BalanceElt(bstruct, refine ) )
    {
         /* If user defined load balancing is not in use, then 
          * redistribution is completed using the load balancing
          * threshold flag. Load balancing (block) is only completed
          * if the percentage increase in the number of total time
          * steps exceeds this value. The default is -1 indicating block
          * load balancing always occurs. */
          fine_gupper = _braid_BalanceElt( bstruct, refined_gupper );
          coarse_gupper = _braid_BalanceElt( bstruct, coarse_gupper); 
          
          braid_Real percent_increase = (double) ( fine_gupper - coarse_gupper )/ ( (double) ( fine_gupper + 1 ) ) ; 
          braid_Real load_bal_threshold = LOAD_BAL_THRES_EXPERIMENTAL;

          if ( load_bal_threshold < percent_increase )
          {
             _braid_GetBlockDistInterval(core, fine_gupper + 1, &fine_ilower, &fine_iupper);         
          }
          else
          { 
             // Copy the values from the refined grid. 
             fine_ilower = _braid_BalanceElt( bstruct, refined_ilower );
             fine_iupper = _braid_BalanceElt( bstruct, refined_iupper );
          }

         _braid_BalanceElt( bstruct, fine_ilower ) = fine_ilower;
         _braid_BalanceElt( bstruct, fine_iupper ) = fine_iupper;
         _braid_BalanceElt( bstruct, fine_gupper ) = fine_gupper;
    }
    else
    {
       //No refinement and no load balancing, so nothing needs to be done. 
       *done = 1; 
    }
    
    return _braid_error_flag;
}

/*--------------------------------------------------------------------------------
 * Build Phase 3b: Get Map from coarse grid to the new fine grid 
 *--------------------------------------------------------------------------------*/
braid_Int
_braid_GetFineGridMap(braid_Core            core,
                      braid_Int            *done,
                      _braid_BalanceStruct *bstruct)
{
   
   MPI_Comm     comm     =  _braid_CoreElt(core, comm );
   
   braid_Int    f_iupper = _braid_BalanceElt( bstruct, fine_iupper );
   braid_Int    f_ilower = _braid_BalanceElt( bstruct, fine_ilower );
   braid_Int    f_gupper = _braid_BalanceElt( bstruct, fine_gupper );
   braid_Int    r_iupper = _braid_BalanceElt( bstruct, refined_iupper );
   braid_Int    r_ilower = _braid_BalanceElt( bstruct, refined_ilower );
   braid_Int    *r_ca     = _braid_BalanceElt( bstruct, refined_ca);
   braid_Int    *r_fa     = _braid_BalanceElt( bstruct, refined_fa);
   braid_Real   *r_ta     = _braid_BalanceElt( bstruct, refined_ta);
   braid_Int    *send_map = _braid_BalanceElt( bstruct, send_map);
   
   MPI_Request request, *requests;
   MPI_Status  status, *statuses;

   braid_Real *f_ta_alloc, *f_ta, *send_buffer, *bptr, *recv_buffer;
   braid_Int *f_ca, *send_procs, *send_sizes, *iptr, *recv_map; 
   braid_Int size, isize, ii, r_i, proc, m, recv_size, r_ii, prevproc;
   braid_Int f_next, f_first, f_npoints, nsends, r_npoints, nreceived, j, f_i, f_ii; 

  
   /* Post f_next receive */
   f_next = -1;
   f_npoints = f_iupper - f_ilower + 1; 
   r_npoints = r_iupper - r_ilower + 1; 
   if (f_npoints > 0)
   {
      f_next = f_gupper+1;
      if (f_iupper < f_gupper)
      {
         MPI_Irecv(&f_next, 1, braid_MPI_INT, MPI_ANY_SOURCE, 3, comm, &request);
      }
   }

   recv_map = _braid_CTAlloc( braid_Int, f_npoints);
   f_ca = _braid_CTAlloc( braid_Int, f_npoints);
   f_ta_alloc = _braid_CTAlloc( braid_Real, f_npoints+2);
   f_ta = f_ta_alloc + 1;
   
   
   /* Compute send information and send f_next info */
   size = 2*sizeof(braid_Int);         /* size of two integers */
   isize = 0;
   _braid_NBytesToNReals(size, isize); /* convert to units of braid_Real */
   send_procs  = _braid_CTAlloc(braid_Int,  r_npoints);
   send_sizes  = _braid_CTAlloc(braid_Int,  r_npoints);
   send_buffer = _braid_CTAlloc(braid_Real, r_npoints*(1+isize+1));
   bptr = send_buffer;
   nsends = -1;

   prevproc = send_map[-1];
   ii = 0;
   for (r_ii = 0; r_ii < r_npoints; r_ii++)
   {
      r_i = r_ilower + r_ii;
      proc = send_map[r_ii];
      if ((proc != prevproc) || (nsends < 0))
      {
         nsends++;
         send_procs[nsends] = proc;
         bptr++; /* leave room for size value */

         if ((proc != prevproc) && (prevproc > -1))
         {
            /* Send f_next info */
            MPI_Send(&r_fa[ii], 1, braid_MPI_INT, prevproc, 3, comm);
         }
         prevproc = proc;
      }
      send_sizes[nsends] += (isize+1);

      iptr = (braid_Int *) bptr;
      iptr[0] = r_i;
      iptr[1] = r_ca[r_ii];
      bptr += isize;
      bptr[0] = r_ta[r_ii];
      bptr++;

      /* Update f_next info */
      if (r_fa[ii] == r_i)
      {
         ii++;
      }
   }
   nsends++;

   requests = _braid_CTAlloc(MPI_Request, nsends);
   statuses = _braid_CTAlloc(MPI_Status,  nsends);

   /* Post sends (do this first, since we will poll on receives) */
   bptr = send_buffer;
   for (m = 0; m < nsends; m++)
   {
      size = send_sizes[m];
      proc = send_procs[m];
      bptr[0] = (braid_Real) size; /* insert size at the beginning */
      MPI_Isend(bptr, (1+size), braid_MPI_REAL, proc, 4, comm, &requests[m]);
      bptr += (1+size);
   }

   /* Post receives */
   recv_size = f_npoints*(1+isize+1); /* max receive size */
   recv_buffer = _braid_CTAlloc(braid_Real, recv_size);
   nreceived = 0;
   while (nreceived < f_npoints)
   {
      /* post receive from arbitrary process (should always get at least one) */
      bptr = recv_buffer;
      size = recv_size;
      MPI_Recv(bptr, size, braid_MPI_REAL, MPI_ANY_SOURCE, 4, comm, &status);

      size = (braid_Int) bptr[0];
      bptr++;
      for (j = 0; j < size; j += (isize+1))
      {
         iptr = (braid_Int *) bptr;
         f_i = iptr[0];
         f_ii = f_i - f_ilower;
         recv_map[f_ii] = status.MPI_SOURCE;
         f_ca[f_ii] = iptr[1];
         bptr += isize;
         f_ta[f_ii] = bptr[0];
         bptr++;
         nreceived++;
      }

   }

   /* Finish sends and f_next receive */
   MPI_Waitall(nsends, requests, statuses);
   if (f_npoints > 0)
   {
      if (f_iupper < f_gupper)
      {
         MPI_Wait(&request, &status);
      }
   }

   /* Compute f_first */
   f_first = f_next;
   for (f_ii = 0; f_ii < f_npoints; f_ii++)
   {
      if (f_ca[f_ii] > -1)
      {
         f_first = f_ilower + f_ii;
         break;
      }
   }
   
   /* Transfer pointer ownership to the balance struct */
   _braid_BalanceElt( bstruct , fine_first  )     = f_first;
   _braid_BalanceElt( bstruct , fine_next   )     = f_next;
   _braid_BalanceElt( bstruct , fine_ca)          = f_ca;
   _braid_BalanceElt( bstruct , fine_ta)          = f_ta;
   _braid_BalanceElt( bstruct , fine_ta_alloc)    = f_ta_alloc;
   _braid_BalanceElt( bstruct , recv_map)         = recv_map;

   /* Free up memory */
   _braid_TFree(requests);
   _braid_TFree(statuses);
   _braid_TFree(send_procs);
   _braid_TFree(send_sizes);
   _braid_TFree(send_buffer);
   _braid_TFree(recv_buffer);

   return _braid_error_flag;
}
/*--------------------------------------------------------------------------------
 * Build Phase 4: Get the communication map   
 *-------------------------------------------------------------------------------*/

braid_Int 
_braid_BuildCommunicationMap(braid_Core            core,
                             braid_Int           *done,
                            _braid_BalanceStruct *bstruct)
{
    
   braid_Int cf, max_levels, min_coarse, g_lower, g_upper, nlevels;
   braid_Int rank, comm_size, i ,j, k, index, level;
   braid_Int haves, needs, bals, sends, response, low, high, sent;
   braid_Int  r_npoints, assumed_ilower, assumed_iupper, rilower_m1;
   braid_Int num_messages, num_needs, num_haves, curr_proc, new_proc, num_recvs;

   braid_Int *gupper, *cfactor, *iupper, *ilower, *owners;
   braid_Int *need_procs, *need_indices,  *have_procs, *have_indices,   *load_sends;
   braid_Int *send_buffer_1, *recv_buffer, *resp_buffer_1, *recv_indices;
   braid_Int *left_procs, *right_procs, *send_map_alloc, *send_map;
   braid_Int **send_buffer, **resp_buffer;

   /* Extract from the Balacne Structure */
   braid_Int fine_gupper    = _braid_BalanceElt( bstruct, fine_gupper );
   braid_Int fine_ilower    = _braid_BalanceElt( bstruct, fine_ilower );
   braid_Int fine_iupper    = _braid_BalanceElt( bstruct, fine_iupper );
   braid_Int refined_iupper = _braid_BalanceElt( bstruct, refined_iupper );   
   braid_Int refined_gupper = _braid_BalanceElt( bstruct, refined_gupper);
   braid_Int refined_ilower = _braid_BalanceElt( bstruct, refined_ilower );
   braid_Int coarse_gupper  = _braid_BalanceElt( bstruct, coarse_gupper);


   MPI_Request *send_requests, *resp_requests, request;
   MPI_Status status;
   MPI_Comm comm = _braid_CoreElt( core, comm );
   MPI_Comm_size( comm, &comm_size );
   MPI_Comm_rank( comm, &rank );

   /* Get the number of levels on the new fine grid */
   max_levels = _braid_CoreElt( core, max_levels );
   min_coarse = _braid_CoreElt( core, min_coarse );
   g_lower = 0;
   g_upper = fine_gupper;
   nlevels = 0;

   while (  nlevels + 1 < max_levels )
   {
      _braid_GetCFactor(core, nlevels, &cf);
      _braid_ProjectInterval( g_lower, g_upper, 0, cf, &g_lower, &g_upper );
      _braid_MapFineToCoarse( g_upper, cf, g_upper );
      if ( g_upper - g_lower >= min_coarse )
      {
         nlevels++;
      }
      else
      {
         break;
      }
   }
   nlevels++;

   /* Allocate the communication maps inside the structure */
   right_procs = _braid_CTAlloc(  braid_Int , nlevels );
   left_procs = _braid_CTAlloc(  braid_Int , nlevels );
   send_map_alloc = _braid_CTAlloc( braid_Int, refined_iupper - refined_ilower + 2 );
   send_map = send_map_alloc + 1;
   
   _braid_BalanceElt( bstruct, right_procs ) = right_procs;
   _braid_BalanceElt( bstruct, left_procs ) = left_procs;
   _braid_BalanceElt( bstruct, send_map_alloc ) = send_map_alloc;
   _braid_BalanceElt( bstruct, send_map ) = send_map;

   /* Get the new grid heirarchy information. This could maybe be done while getting
    * the number of levels above */
   gupper =  _braid_CTAlloc( braid_Int, nlevels );
   cfactor = _braid_CTAlloc( braid_Int, nlevels );
   ilower =  _braid_CTAlloc( braid_Int, nlevels );
   iupper =  _braid_CTAlloc( braid_Int, nlevels );
   gupper[0] = fine_gupper;
   iupper[0] = fine_iupper;
   ilower[0] = fine_ilower;
   
   _braid_GetCFactor( core, 0, &cfactor[0] );
   for ( i = 1; i < nlevels; i++ )
   {
      _braid_ProjectInterval( ilower[i-1], iupper[i-1] ,0, cfactor[i-1], &ilower[i], &iupper[i] );
      _braid_MapFineToCoarse( ilower[i], cfactor[i-1], ilower[i] );
      _braid_MapFineToCoarse( iupper[i], cfactor[i-1], iupper[i]);
      _braid_ProjectInterval( 0, gupper[i-1], 0, cfactor[i-1], &k, &gupper[i] );
      _braid_MapFineToCoarse( gupper[i], cfactor[i-1], gupper[i] );
      _braid_GetCFactor(core, i, &cfactor[i] );
   }

   /* Get a list of indices that are my neigbours on each level 
    * need_indices[l] is the index processor will send messages to on level l. 
    * recv indices are the indicies processor will recv messages from on level l */
   need_indices = _braid_CTAlloc( braid_Int, nlevels );
   need_procs = _braid_CTAlloc( braid_Int, nlevels + 1 );
   recv_indices = _braid_CTAlloc( braid_Int, nlevels );
   
   num_needs = 0;
   num_recvs = 0;
   for ( i = 0; i < nlevels; i++ )
   {
      index = iupper[i] + 1;
      if ( index <= gupper[i] && iupper[i] >= ilower[i])
      {
         for ( j = i; j > 0; j-- )
         {
            _braid_MapCoarseToFine( index, cfactor[j-1], index );
         }
         _braid_GetBlockDistProc( core, gupper[0] + 1, index, &need_procs[num_needs] );
         need_indices[num_needs++] = index;
      }
      if ( !_braid_BalanceElt( bstruct, lbalance ) )
      {
         index = ilower[i] - 1;
         if ( index >= 0 && iupper[i] >= ilower[i] )
         {
            for ( j = i; j > 0; j--)
            {
               _braid_MapCoarseToFine( index, cfactor[j-1], index );
            }
            recv_indices[num_recvs++] = index;
         }

      }
   }

   /* If not load balencing using weights then fill maps and return. Note that
    * the send and recv procs arrays are not allocated. The InitHeirarchy function
    * sets neighbours based on the weighted dist in the advent that these vectors 
    * are null. */
   
   /* Experimental flag turns of weighted dist load balancing in the case that the number 
    * of new time steps is "small".*/
   braid_Real load_bal_thres = LOAD_BAL_THRES_EXPERIMENTAL; 
   braid_Real percent_increase = (refined_gupper - coarse_gupper ) / ( coarse_gupper + 1 );
   braid_Int weighted_lbal = ( percent_increase > load_bal_thres ) ? 1 : 0; 
   
   if (  !_braid_BalanceElt( bstruct, lbalance ) && weighted_lbal && 1 )
   {
      /* Set Defualt value */
      for ( i = 0; i < nlevels; i++ )
      {
         left_procs[i] = -1;
         right_procs[i] = -1;
      }

      /* Get the processor that owns index to right on each level */
      for ( i = 0; i < num_needs; i++ )
      {
         _braid_GetBlockDistProc( core, fine_gupper + 1, need_indices[i], &right_procs[i]);
      }

      /* Get the processor that owns index to the left on each level */
      for ( i = 0; i < num_recvs; i++ )
      {
         _braid_GetBlockDistProc( core, fine_gupper + 1, recv_indices[i], &left_procs[i]);
      }

      /* Get the processor that owns, on the new fine grid, each index on the refined grid */
      for ( i = refined_ilower; i <= refined_iupper ; i++ )
      {
         _braid_GetBlockDistProc( core, fine_gupper + 1, i, &send_map[i-refined_ilower]) ;
      }

      /* Get the owner of the index to the left of refined_ilower ( needed for FRefine ? ) */
      _braid_GetBlockDistProc( core, fine_gupper + 1, refined_ilower - 1 , &rilower_m1 );
      send_map[-1] = rilower_m1;

      _braid_TFree( recv_indices );
      _braid_TFree( gupper );
      _braid_TFree( ilower );
      _braid_TFree( iupper );
      _braid_TFree( cfactor );
      _braid_TFree( need_procs );
      _braid_TFree( need_indices );
       return _braid_error_flag; 
   }

   /* The distribution is no longer uniform. Set flag to make this permanent */
   //BRAID_UNIFORM_DISTRIBUTION = 0;

   /* A recv map is built later on using information gathered as 
    * part of the load balancing proceedure. There is no need to explicitly 
    * communicate these values. */
   _braid_TFree( recv_indices );

   /* Get this processors ownership on the assumed partition (communication). owners[i] indicates the processor
    * that owns the index " i-assumed_ilower " */
   owners = NULL;
   _braid_GetPartition( core, fine_ilower, fine_iupper, fine_gupper , &assumed_ilower, &assumed_iupper, &owners );

   /* Calculate how many messages will be recieved. 
    *
    * Each processor will recieve one message for every point it owns on the every grid of 
    * the assumed partition ( expect index = 0 ). This allows us to send nearest neighbour information
    * for the new fine grid. 
    * In addition each processor recieves one message for every fine grid point in the assumed partition. This
    * information is needed for load balancing. Note that messages are combined into buffers whereever possible
    * to reduce communication */

   num_messages = 0;
   for ( i = _braid_max( assumed_ilower, 1 ); i <= assumed_iupper ; i++ )
   {
      index = i;
      level = 0;
      num_messages++;
      while ( level < nlevels -1 && !(index % cfactor[level])  )
      {
         index = index/cfactor[level++];
         num_messages++;
      }
   }
   num_messages += assumed_iupper - assumed_ilower + 1;

   /* Find the processor that owns each refined grid element on the assumed aprtition. The assumed 
    * partition will reply with the actual owner of these indicies */
   r_npoints = refined_iupper - refined_ilower + 1 ;
   load_sends = _braid_CTAlloc( braid_Int, r_npoints + 1);
   for ( i = refined_ilower; i <= refined_iupper; i++ )
   {
      _braid_GetBlockDistProc( core, gupper[0] + 1, i , &load_sends[i - refined_ilower] );
   }
   load_sends[r_npoints] = -1;

   /* Get a list of indices that this processor owns on both the coarse and fine grid. For these indicies we 
    * send a message to the owner on the assumed partition, but do not wait for a response. */
   have_indices = _braid_CTAlloc( braid_Int, iupper[0] - ilower[0] + 10 );
   have_procs = _braid_CTAlloc( braid_Int, iupper[0] - ilower[0] + 10);
   num_haves = -1;
   curr_proc = -1;
   new_proc = -1 ;
   for ( i = ilower[0] + 1 ; i <= iupper[0]; i++ )
   {
      level = 0;
      index = i;
      _braid_GetBlockDistProc( core, gupper[0] + 1, index, &new_proc );
      if ( new_proc != curr_proc )
      {
         curr_proc = new_proc;
         have_procs[++num_haves] = new_proc;
      }
      have_indices[num_haves]++;
      while ( level < nlevels -1 && !( index % cfactor[level] ) )
      {
         index = index/cfactor[level++];
         if ( index >  ilower[level] )
         {
            have_indices[num_haves]++;
         }
      }
   }
   num_haves++;


   /* Get the lowest and highest processors this processor sends to. */
   low = braid_Int_Max;
   high = -1;
   if ( num_haves > 0 )
   {
      low = _braid_min( have_procs[0], low );
      high = _braid_max( have_procs[num_haves-1], high );
   }
   if ( num_needs > 0 )
   {
      low = _braid_min( need_procs[0], low );
      high = _braid_max( need_procs[num_needs-1], high );
   }
   if ( r_npoints > 0 )
   {
      low = _braid_min( load_sends[0], low );
      high = _braid_max( load_sends[r_npoints -1], high );
   }
   braid_Int num_send_procs = _braid_max( high - low + 1 , 0 );

   /* Build the send buffer, send it, and post a recv for the response if needed */
   haves = 0;
   needs = 0;
   bals  = 0;
   response = 0;
   have_procs[ num_haves ] = -1;
   need_procs[ num_needs ] = -1;
   send_buffer = _braid_CTAlloc( braid_Int*, num_send_procs  );
   resp_buffer = _braid_CTAlloc( braid_Int*, num_send_procs );
   send_requests = _braid_CTAlloc( MPI_Request, num_send_procs );
   resp_requests = _braid_CTAlloc( MPI_Request, num_send_procs );

   sends = -1;
   sent = 1;
   braid_Int *bptr;

   for ( i = low; i <= high; i++ )
   {
      if ( sent )
      {
         /* Allocate a send buffer for the next processor */
         sends++;
         send_buffer[sends] = _braid_CTAlloc( braid_Int, 11 + nlevels + r_npoints  );
         sent = 0;
         bptr = &(send_buffer[sends][3]);
      }
      if ( i == have_procs[haves] )
      {
         /* Add the indicies that do not require a response */
         send_buffer[sends][0] = have_indices[haves++];
         sent = 1;
      }
      while( i == need_procs[needs] )
      {
         /* Add the indicies that require a response. These are the indicies that this
          * processor needs to know who owns for neighest neighbour information */         
         *bptr = need_indices[needs++];
         send_buffer[sends][1]++;
         bptr++;
         sent = 1;
      }
      while ( i == load_sends[bals] )
      {
         /* Add the indicies for the load balancing portion. These are the indicies the processor 
          * needs to know who owns for load balancing inforation */
         *bptr = refined_ilower + bals;
         send_buffer[sends][2]++;
         bals++;
         bptr++;
         sent = 1;
      }
      if ( sent )
      {
         /* Send it */
         j =  3 + send_buffer[sends][1] + send_buffer[sends][2];
         MPI_Isend( send_buffer[sends] , j , braid_MPI_INT, i, 20, comm, &send_requests[sends] );
         
         /* Post a recv if we need one */
         if (send_buffer[sends][1] + send_buffer[sends][2] > 0 )
         {
            /* Recv will one index for each of the ones asked for ?? */
            resp_buffer[response] = _braid_CTAlloc( braid_Int, j );
            MPI_Irecv( resp_buffer[response] , j , braid_MPI_INT, i, 21, comm, &resp_requests[response] );
            response++;
         }
      }
   }

   _braid_TFree ( need_procs );
   _braid_TFree ( have_procs );
   _braid_TFree ( have_indices );
   _braid_TFree ( need_indices );
   _braid_TFree( load_sends );

   /* Recieve all messages and respond as needed. This is assumed partition */
   i = 0;
   braid_Int recv_buffer_size = 10 + nlevels + assumed_iupper - assumed_ilower ;
   recv_buffer = _braid_CTAlloc( braid_Int, recv_buffer_size );
   resp_buffer_1 = _braid_CTAlloc( braid_Int, recv_buffer_size );
   braid_Int *rptr;
   while (num_messages > 0 )
   {
      MPI_Recv( recv_buffer, recv_buffer_size , braid_MPI_INT, MPI_ANY_SOURCE, 20, comm, &status );
      
      /* First, calculate the number of messages remaining after this recv  */
      num_messages -= (recv_buffer[0] + recv_buffer[1] + recv_buffer[2]);
      
      /* Respond if requested */
      if ( recv_buffer[1] + recv_buffer[2] > 0 )
      {

         /* Response buffer */
         resp_buffer_1[0] = recv_buffer[1]; 
         resp_buffer_1[1] = recv_buffer[2]; 
         bptr = &(resp_buffer_1[2]); 
         rptr = &(recv_buffer[3]);

         for ( i = 0; i < recv_buffer[1] + recv_buffer[2] ; i++ )
         {
            *bptr = owners[ *rptr - assumed_ilower ];
            bptr++; rptr++;
         }
         
         /* Ideally this would be an non-blocking send, but i am strugling to determine a reasonable
          * memory allocation for the resp_buffer */ 
         j = 2 + recv_buffer[1] + recv_buffer[2];
         MPI_Send( resp_buffer_1, j, braid_MPI_INT, status.MPI_SOURCE, 21, comm );
      }
   }

   if ( sends >= 0 )
   {
      MPI_Waitall( sends , send_requests , MPI_STATUS_IGNORE );
   }
   if ( response > 0 )
   {
      MPI_Waitall( response, resp_requests, MPI_STATUS_IGNORE );
   }

   _braid_TFree( owners );
   _braid_TFree( recv_buffer );
   _braid_TFree( resp_buffer_1 );
   _braid_TFree( send_requests );
   _braid_TFree( resp_requests );
   for( i = 0; i <= sends; i++ )
   {
      _braid_TFree( send_buffer[i] );
   }
   _braid_TFree( send_buffer );

   /*Extract the response information. */
   braid_Int *right_procs_ptr = right_procs;
   braid_Int *send_map_ptr = send_map;
   for ( i = 0; i < response; i++ )
   {
      bptr = resp_buffer[i];
      rptr = &(bptr[2]);
      
      memcpy( right_procs_ptr , rptr, sizeof(int)*bptr[0] );
      rptr += bptr[0]; 
      right_procs_ptr += bptr[0] ;

      memcpy( send_map_ptr, rptr, sizeof(int)*bptr[1] );
      send_map_ptr += bptr[1] ;
      _braid_TFree( resp_buffer[i] );
   }
   _braid_TFree( resp_buffer );

   
   /* Finish the send map */
   rilower_m1 = -1;
   if ( refined_ilower > 0 && refined_ilower <= refined_iupper)
   {
      _braid_GetProcLeftOrRight( core, 0, -1 , &curr_proc );
      MPI_Irecv( &rilower_m1 , 1, braid_MPI_INT, curr_proc, 23, comm, &request);
   }
   if ( refined_iupper < gupper[0] && refined_ilower <= refined_iupper )
   {
      _braid_GetProcLeftOrRight( core, 0, 1, &curr_proc );
      MPI_Send( &send_map[ refined_iupper - refined_ilower ], 1, braid_MPI_INT, curr_proc, 23, comm );
   }
   if ( refined_ilower > 0 && refined_ilower <= refined_iupper  )
   {
      MPI_Waitall ( 1, &request, MPI_STATUS_IGNORE );
   }
   send_map[-1] = rilower_m1;

   
   /*Nearest neighbour communication. As it stands, XBraid needs to know both who
    * it will recv from and who it sends to. Above, we discoved the processor to our
    * right on all levels. In this section, we send a message to those procs to let 
    * them know who they will be recving from */

   sends = 0;
   send_buffer_1 = _braid_CTAlloc( braid_Int, nlevels );
   send_requests = _braid_CTAlloc( MPI_Request, nlevels );
   for ( i = 0; i < nlevels; i++ )
   {
      if ( right_procs[i] > 0 && rank < comm_size - 1 && ilower[0] <= iupper[0] )
      {
         send_buffer_1[sends] = i;
         MPI_Isend( &send_buffer_1[sends] , 1, braid_MPI_INT, right_procs[i], 22, comm , &send_requests[sends]);
         sends++;
      }
   }

   for ( i = 0; i < nlevels; i++ )
   {
      if ( iupper[i] - ilower[i] + 1 > 0 && rank > 0 )
      {
         MPI_Recv( &index, 1, braid_MPI_INT, MPI_ANY_SOURCE, 22, comm , &status);
         left_procs[index] = status.MPI_SOURCE;
      }
      else
         left_procs[i] = -1;
   }

   if ( sends > 0 )
   {
      MPI_Waitall( sends, send_requests , MPI_STATUS_IGNORE );
   }

   _braid_TFree( send_buffer_1 );
   _braid_TFree( send_requests );

   _braid_TFree( gupper );
   _braid_TFree( iupper );
   _braid_TFree( ilower );
   _braid_TFree( cfactor );

   _braid_TFree( send_buffer );
   _braid_TFree( recv_buffer );

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------------
 * Destroy the Balance structure
 * -------------------------------------------------------------------------------*/

braid_Int
_braid_DestroyBalanceStruct( _braid_BalanceStruct *bstruct )
{
   braid_Int *left               = _braid_BalanceElt( bstruct, left_procs );
   braid_Int *right              = _braid_BalanceElt( bstruct, right_procs );
   braid_Int *send               = _braid_BalanceElt( bstruct, send_map_alloc );
   braid_Int *recv               = _braid_BalanceElt( bstruct, recv_map);
   braid_Int  *refined_ca        = _braid_BalanceElt( bstruct, refined_ca );           
   braid_Int  *refined_fa        = _braid_BalanceElt( bstruct, refined_fa);
   braid_Real *refined_taalloc   = _braid_BalanceElt( bstruct, refined_taalloc);
   braid_Int  *fine_ca           = _braid_BalanceElt( bstruct, fine_ca);
   braid_Real *fine_ta_alloc     = _braid_BalanceElt( bstruct, fine_ta_alloc);
    
   _braid_TFree( left );
   _braid_TFree( right );
   _braid_TFree( send );
   _braid_TFree( recv );
   _braid_TFree( refined_ca);
   _braid_TFree( refined_fa);
   _braid_TFree( refined_taalloc);
   _braid_TFree( fine_ca);
   _braid_TFree( fine_ta_alloc);

   _braid_TFree( bstruct );

   return _braid_error_flag;
}

braid_Int
_braid_GetBlockDistInterval(braid_Core core,
                            braid_Int  npoints,
                            braid_Int *ilower_ptr,
                            braid_Int *iupper_ptr)
{
    MPI_Comm comm = _braid_CoreElt( core, comm );    
    braid_Int comm_min  = _braid_CoreElt(core, comm_min);
    braid_Int myid      = _braid_CoreElt(core, myid);
    
    braid_Int  ilower,iupper, nprocs, cfactor, ncpoints;    
    MPI_Comm_size(comm, &nprocs);

    if ( !comm_min ) 
    {
       /*Distribute points evenly across the time processors */
       _braid_GetBlockDistInterval_basic( npoints, nprocs, myid, &ilower, &iupper);             
    }
    else
    {
       /*Distribute the C points evenly across the processors */
       _braid_GetCFactor(core,0,&cfactor); 
       ncpoints = (npoints - 1)/cfactor + 1 ;  
       _braid_GetBlockDistInterval_basic( ncpoints, nprocs, myid, &ilower, &iupper);
       
       if ( ilower == ncpoints )
       {
          /* Empty Processor */
          ilower = npoints;
          iupper = npoints - 1; 
       }
       else
       {
            /* Map back to the fine grid */
            _braid_MapCoarseToFine( ilower, cfactor, ilower );       
            _braid_MapCoarseToFine( iupper, cfactor, iupper );
            iupper = _braid_min( iupper + cfactor - 1, npoints-1 ); 
       }
    }

    *ilower_ptr = ilower;
    *iupper_ptr = iupper;

#if DEBUG
    printf( " PROC %d ILOWER %d IUPPER %d \n ", myid, ilower, iupper );
#endif

    return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns the index interval for 'proc' in a blocked data distribution
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetBlockDistInterval_basic(braid_Int   npoints,
                                  braid_Int   nprocs,
                                  braid_Int   proc,
                                  braid_Int  *ilower_ptr,
                                  braid_Int  *iupper_ptr)
{
   braid_Int  ilower, iupper, quo, rem, p;

   quo = npoints/nprocs;
   rem = npoints%nprocs;

   p = proc;
   ilower = p*quo + (p < rem ? p : rem);
   p = proc+1;
   iupper = p*quo + (p < rem ? p : rem) - 1;

   *ilower_ptr = ilower;
   *iupper_ptr = iupper;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns the processor that owns 'index' in a blocked data distribution
 * (returns -1 if index is out of range)
 *----------------------------------------------------------------------------*/

braid_Int 
_braid_GetBlockDistProc(braid_Core  core,
                        braid_Int   npoints,
                        braid_Int   index,
                        braid_Int   *proc_ptr)
{
   MPI_Comm comm = _braid_CoreElt(core, comm);
   braid_Int comm_min = _braid_CoreElt(core, comm_min );
   
   braid_Int nprocs, proc, index_copy, cfactor, ncpoints;
   MPI_Comm_size(comm, &nprocs);

   if ((index < 0) || (index > (npoints-1)))
   {
      proc = -1;
   }
   else if (!comm_min)
   {
      _braid_GetBlockDistProc_basic( npoints, nprocs, index, &proc );
   }
   else
   {
      //Find the index of the last C point
      _braid_GetCFactor(core,0,&cfactor);
      index_copy = index;
      while( _braid_IsFPoint(index_copy, cfactor) )
      {
         index_copy--;
      }
      _braid_MapFineToCoarse(index_copy,cfactor,index_copy);
      ncpoints = (npoints-1)/cfactor + 1;      
      _braid_GetBlockDistProc_basic( ncpoints, nprocs, index_copy, &proc);
   }

#if DEBUG
   printf(" INDEX %d PROC %d \n ", index, proc );
#endif
   *proc_ptr = proc;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns the processor that owns 'index' in a blocked data distribution
 * (returns -1 if index is out of range)
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetBlockDistProc_basic(braid_Int   npoints,
                        braid_Int   nprocs,
                        braid_Int   index,
                        braid_Int  *proc_ptr)
{
   braid_Int      proc, quo, rem, p, q;

   /* Compute processor number */
   if ((index < 0) || (index > (npoints-1)))
   {
      proc = -1;
   }
   else
   {
      quo = npoints/nprocs;
      rem = npoints%nprocs;

      if (quo > 0)
      {
         p = index/(quo+1);
         q = (index - rem*(quo+1))/quo;
         proc = (p < rem ? p : rem+q);
      }
      else
      {
         proc = index;
      }
   }

   *proc_ptr = proc;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Initialize the data structure used to store the global wfactor information
 *----------------------------------------------------------------------------*/

braid_Int
_braid_WeightedStructInit( _braid_WeightedStruct *wstruct )
{
    _braid_WeightedElt( wstruct, local_sum ) = 0;
    _braid_WeightedElt( wstruct, local_min ) = -1;
    _braid_WeightedElt( wstruct, local_max ) = -1;
    _braid_WeightedElt( wstruct, global_max ) = -1;
    _braid_WeightedElt( wstruct, global_min) = -1;
    _braid_WeightedElt( wstruct, global_sum ) = 0;
    _braid_WeightedElt( wstruct, local_start ) = -1;
    _braid_WeightedElt( wstruct, local_stop ) = -1;

    return _braid_error_flag;
}
/*----------------------------------------------------------------------------
 * Destroy the weighted data structure
 *----------------------------------------------------------------------------*/

braid_Int
_braid_WeightedStructDestroy( _braid_WeightedStruct *wstruct )
{
    _braid_TFree( wstruct );
    return _braid_error_flag;
}

/*----------------------------------------------------------------------------
  * Build the correct structure to get load balancing working 
* *----------------------------------------------------------------------------*/

braid_Int 
_braid_BuildWeightedStruct(braid_Core             core,
                           braid_Int             *done,
                          _braid_BalanceStruct  *bstruct,
                          _braid_WeightedStruct *wstruct)
{

    MPI_Comm    comm            = _braid_CoreElt( core, comm );    
    braid_Real *wfactors        = _braid_CoreElt(core, wfactors);    
    braid_Int   refined_npoints = _braid_BalanceElt( bstruct, refined_npoints) ;
    braid_Int   refined_ilower      = _braid_BalanceElt(bstruct, refined_ilower);
    braid_Int   refined_iupper      = _braid_BalanceElt(bstruct, refined_iupper);
    braid_Int   refined_gupper      = _braid_BalanceElt(bstruct, refined_gupper);

    braid_Int  myid, recving, sending, cfactor, i,j ;
    braid_Real recv_sum, send_sum, comm_min_wmax, comm_min_wmin;
    braid_Real local_max, local_sum, local_min, *send_data, *recv_data;
    MPI_Request send_request, recv_request;
    
    MPI_Comm_rank(comm, &myid);
    _braid_GetCFactor(core, 0, &cfactor);     
       
    /* If using the comm_min function then the weights must represent the
     * weights of the C-point interval. This section sums the weights across
     * each C-point interval, then sets the weight at the F points to zero. 
     * TODO Might be a better way to do this 
     */
    
    comm_min_wmin = -1;
    comm_min_wmax = -1;
    if (_braid_CoreElt(core, comm_min) && refined_iupper >= refined_ilower  )
    { 
       /* If the first point on the next processor is not a Cpoint then we will recieve 
        * a message regarding the sum of the weights on the interval */
       recv_sum = 0;
       recving = 0;
       if ( !_braid_IsCPoint( refined_iupper + 1 , cfactor ) && refined_iupper != refined_gupper )
       { 
          MPI_Irecv( &recv_sum, 1, MPI_DOUBLE, myid+1 , 30, comm, &recv_request );
          recving = 1;
       }

       /* If the first point on this processor is not a C point, then send a message 
        * regarding the sum on this interval. 
        * TODO This only works if comm_min is set true, and remains true, the entire 
        * time. If comm_min is set true after x iterations, this wont work.
        * */
       send_sum = 0;
       sending = 0;
       i = 0;       
       if ( !_braid_IsCPoint( refined_ilower , cfactor ) )
       {
          while ( !_braid_IsCPoint( refined_ilower + i, cfactor ) )
          {
            send_sum += wfactors[i];
            wfactors[i] = 0;
            comm_min_wmin = (  wfactors[i] < comm_min_wmin || comm_min_wmin < 0 ) ? wfactors[i] : comm_min_wmin;
            comm_min_wmax = ( wfactors[i] > comm_min_wmax ) ? wfactors[i] : comm_min_wmax;
            i++;
          }
          MPI_Isend( &send_sum, 1, MPI_DOUBLE, myid - 1, 30, comm, &send_request );
          sending = 1;
       }

       i = j = 0;
       while ( i < refined_npoints )
       {
          if (!_braid_IsCPoint( i + refined_ilower, cfactor) )
          {
             comm_min_wmin = (  wfactors[i] < comm_min_wmin || comm_min_wmin < 0 ) ? wfactors[i] : comm_min_wmin;
             comm_min_wmax = ( wfactors[i] > comm_min_wmax ) ? wfactors[i] : comm_min_wmax;
             wfactors[j] += wfactors[i]; 
             wfactors[i] = 0;
          }
          else
             j = i;
          i++;
       }

       if (recving)
       {
          MPI_Waitall( 1 , &recv_request , MPI_STATUS_IGNORE);
          wfactors[j] += recv_sum;
       }
       if (sending)
       {
          MPI_Waitall( 1, &send_request, MPI_STATUS_IGNORE);
       }
    }

    /* Get the local max and min information. If comm_min is used, these
     * will be the updated, C-interval weighted values. */
    local_max = 0;
    local_min = -1;
    local_sum = 0;

    for ( i = 0; i < refined_npoints; i++ )
    {
        if ( local_max < wfactors[i] )
        {
           local_max = wfactors[i];
        }
        if ( ( local_min > wfactors[i] && wfactors[i] > 0 ) || local_min < 0 )
        {
           local_min = wfactors[i];
        }
        local_sum += wfactors[i];
    }

    /* Do the allreduce */
    send_data = _braid_CTAlloc( braid_Real, 5 );
    recv_data = _braid_CTAlloc( braid_Real, 5 );
    send_data[0] = local_sum;
    send_data[1] = comm_min_wmax;
    send_data[2] = comm_min_wmin;
    send_data[3] = local_max;
    send_data[4] = local_min;

    /* sum_max_min is a user defined allreduce function. that takes the
     * sum[0],max[1],min[2],max[3],min[4] over all processors */ 
    MPI_Op sum3_max_min;
    MPI_Op_create( (MPI_User_function*) Sum3MaxMin , 1 , &sum3_max_min );
    MPI_Allreduce(send_data, recv_data, 5, braid_MPI_REAL, sum3_max_min, comm);
    MPI_Op_free( &sum3_max_min );

    /* Build the structure */
    _braid_WeightedStructInit( wstruct );        
    _braid_WeightedElt( wstruct, local_sum ) = local_sum;
    _braid_WeightedElt( wstruct, local_max ) = local_max;
    _braid_WeightedElt( wstruct, local_min ) = local_min;
    _braid_WeightedElt( wstruct, global_sum ) = recv_data[0];
    _braid_WeightedElt( wstruct, global_max ) = recv_data[3];
    _braid_WeightedElt( wstruct, global_min ) = recv_data[4];

    /* Weighted balancing is not needed if min weight == max weight */ 
    braid_Int weighted_lbal = _braid_BalanceElt(bstruct, lbalance);
    if (( !_braid_CoreElt(core,comm_min) && recv_data[3] == recv_data[4] ) ||
       ( _braid_CoreElt(core,comm_min) && recv_data[1] == recv_data[2] ))
    {
       weighted_lbal = 0;
    }
    
    _braid_TFree(send_data);
    _braid_TFree(recv_data);

    /* If weighted load balacning, get the cumulative sum */
    *done = 0;   
    if ( weighted_lbal )
    {
       braid_Real csum;
       MPI_Scan(&local_sum, &csum , 1, braid_MPI_REAL, MPI_SUM, comm);
       _braid_WeightedElt( wstruct, local_start )  = csum - local_sum;
       _braid_WeightedElt( wstruct, local_stop )   = csum;
    } 
    else if ( _braid_BalanceElt(bstruct, refine) && !weighted_lbal )
    {
        _braid_BalanceElt( bstruct, lbalance ) = 0;
    }
    else
    {  
      *done = 1; //No refine and no load balacne -- return!
    }

   return _braid_error_flag;
}


/*------------------------------------------------------------------------------------
 * Get the new fine grid interval using a weighted data distribution.
 *----------------------------------------------------------------------------------*/

braid_Int
_braid_GetWeightedFineInterval(braid_Core      core,
                               braid_Int      *done,
                               _braid_BalanceStruct  *bstruct )
{
    MPI_Comm      comm     = _braid_CoreElt( core, comm );
    braid_Real   *wfactors = _braid_CoreElt(core, wfactors);      
   
    braid_Int rank, comm_size, i, j;
    braid_Int *send_buffer_ilower, *send_buffer_iupper, num_recvs, num_sends;
    braid_Int comm_min, refined_gupper, refined_ilower, refined_ncpoints, cfactor;
    braid_Int partition_low, partition_high, num_partitions, fine_iupper, fine_ilower;
    braid_Real goal_load, global_sum, local_sum, c_sum_begin, c_sum_end;
    MPI_Request *send_requests, *recv_requests;
        
    MPI_Comm_rank( comm, &rank );
    MPI_Comm_size( comm, &comm_size );

    /* Build the weighted structure */ 
    _braid_WeightedStruct wstruct;
    _braid_WeightedStructInit( &wstruct );
    _braid_BuildWeightedStruct( core, done, bstruct, &wstruct );    
    
    if ( *done )
    {
       return _braid_error_flag;
    }
    else if ( _braid_BalanceElt(bstruct, lbalance ) == 0 )
    {
         /* If refinement is taking place, without load balancing, then
          * use the block data distrobution to get new interval */ 
          refined_gupper = _braid_BalanceElt( bstruct, refined_gupper );
         _braid_GetBlockDistInterval(core, refined_gupper + 1, &fine_ilower, &fine_iupper);         
         _braid_BalanceElt( bstruct, fine_ilower ) = fine_ilower;
         _braid_BalanceElt( bstruct, fine_iupper ) = fine_iupper;
         _braid_BalanceElt( bstruct, fine_gupper ) = refined_gupper;
         return _braid_error_flag;
    }

    /* Pull everything from the weighted struct for convienence */
    comm_min       = _braid_CoreElt(core, comm_min );
    refined_gupper = _braid_BalanceElt( bstruct, refined_gupper );
    refined_ilower = _braid_BalanceElt( bstruct, refined_ilower );
    global_sum =   _braid_WeightedElt( &wstruct, global_sum );
    local_sum =    _braid_WeightedElt( &wstruct, local_sum );
    c_sum_begin =  _braid_WeightedElt( &wstruct, local_start );
    c_sum_end =    _braid_WeightedElt( &wstruct, local_stop );

    _braid_GetCFactor(core,0,&cfactor);
    refined_ncpoints = refined_gupper/cfactor + 1;
    
    /* num_partitions is the number of interval processor boundaries. For example
     * 3 processors have two interval processor boundaries. */
    if (comm_min)
      num_partitions = _braid_min( comm_size, refined_ncpoints ); 
    else   
      num_partitions = _braid_min( comm_size, refined_gupper + 1 );
    num_partitions--;

    /* Goal load if the goal weight on each processor. The goal load is
     * always bigger than the max weight on any single time point */
    goal_load = global_sum / (num_partitions+1);       

    //Post recvs for the new ilower and iupper. 
    num_recvs = 0;
    recv_requests = _braid_CTAlloc( MPI_Request, 2 );
    if ( rank <= num_partitions )
    {
        fine_ilower = 0;
        fine_iupper = refined_gupper;
        if ( rank < num_partitions )
        {
           MPI_Irecv( &fine_iupper, 1, MPI_INT, MPI_ANY_SOURCE, 31, comm, &recv_requests[num_recvs++] );
        }
        if ( rank > 0 )
        {
           MPI_Irecv( &fine_ilower, 1, MPI_INT, MPI_ANY_SOURCE, 32, comm, &recv_requests[num_recvs++] );
        }
    }
    else
    {
        fine_ilower = refined_gupper + 1;
        fine_iupper = refined_gupper;
        
    }
 
    /* Find any partition boundaries that lie on my current grid. If so, send
       the ilower and iupper information to the corresponding processors.  */
    num_sends = 0;
    j = 0;
    if ( local_sum > 0 )
    {
        partition_low = (braid_Int) ceil( c_sum_begin/goal_load );
        partition_high = _braid_min( num_partitions , (braid_Int) floor(c_sum_end/goal_load ) );
      
        /* If a lower parttion boundary is an integer make sure it only gets accounted for by one processor  */
        if ( partition_low == 0 || fabs( ceil(c_sum_begin/goal_load) - floor(c_sum_begin/goal_load)) < 1e-13 )
        {
           partition_low++;
        }
                
        send_buffer_ilower  = _braid_CTAlloc( braid_Int , (partition_high - partition_low + 1) );
        send_buffer_iupper  = _braid_CTAlloc( braid_Int , (partition_high - partition_low + 1) );
        send_requests = _braid_CTAlloc( MPI_Request, 2*(partition_high-partition_low + 1) ) ;
        
        /* Send the ilower and iupper information */
        for ( i = partition_low; i <= partition_high; i++ )
        {
            while ( c_sum_begin < goal_load*i )
            {
               c_sum_begin += wfactors[j++];
            }

            if (comm_min)
            {
               send_buffer_ilower[(i-partition_low)] = j + cfactor -1 + refined_ilower;
               send_buffer_iupper[(i-partition_low)] = j + cfactor -2 + refined_ilower;
            }
            else
            {
               send_buffer_ilower[(i-partition_low)] = j + refined_ilower;
               send_buffer_iupper[(i-partition_low)] = j-1 + refined_ilower;
            }
            MPI_Isend( &send_buffer_ilower[i-partition_low], 1, braid_MPI_INT, i, 32, comm, &send_requests[num_sends++] );
            MPI_Isend( &send_buffer_iupper[i-partition_low], 1, braid_MPI_INT, i - 1, 31, comm, &send_requests[num_sends++] );
        }
    }

    /* Complete communication */
    if ( num_recvs > 0 )
    {
       MPI_Waitall( num_recvs , recv_requests , MPI_STATUS_IGNORE);
    }
    if ( num_sends > 0 )
    {
        MPI_Waitall( num_sends , send_requests , MPI_STATUS_IGNORE);
    }
#if DEBUG
    printf ("RANK %d ILOWR %d IUPPER %d \n", rank, fine_ilower, fine_iupper );
#endif
    /* Set the fine grid information */
    _braid_BalanceElt( bstruct, fine_ilower ) = fine_ilower;
    _braid_BalanceElt( bstruct, fine_iupper ) = fine_iupper;
    _braid_BalanceElt( bstruct, fine_gupper ) = refined_gupper;

    /* Free mallocs */
    if ( local_sum > 0 )
    {
        _braid_TFree( send_requests );
        _braid_TFree( send_buffer_ilower );
        _braid_TFree( send_buffer_iupper );
    }
    _braid_TFree( recv_requests );

    return _braid_error_flag;
}

/*------------------------------------------------------------------------------
 * Get the assumed partion, returning an array of the owners.
 *-----------------------------------------------------------------------------*/

braid_Int
_braid_GetPartition(braid_Core core,
                    braid_Int ilower,
                    braid_Int iupper,
                    braid_Int gupper,
                    braid_Int *assumed_ilower,
                    braid_Int *assumed_iupper,
                    braid_Int **owner)
{

   MPI_Comm comm = _braid_CoreElt( core, comm );
   
   braid_Int comm_size, rank, i;
   braid_Int sbuf_0, sbuf_1, proc, proc1, num_sends, num_messages;
   braid_Int *send_proc, **send_buffer, *recv_buffer;
   MPI_Request *send_request;
   MPI_Status recv_status;

   MPI_Comm_size( comm, &comm_size );
   MPI_Comm_rank( comm, &rank );

   /* get my assumed partition*/
   _braid_GetBlockDistInterval(core, gupper + 1, assumed_ilower, assumed_iupper );
   braid_Int *owners = _braid_CTAlloc( braid_Int, *assumed_iupper - *assumed_ilower + 1 );

   /* Fill in send buffers that claim ownership of an index in the assumed partition */
   num_sends = iupper - ilower + 1;
   send_buffer = _braid_CTAlloc( braid_Int*, num_sends + 1 );
   send_proc = _braid_CTAlloc( braid_Int, num_sends + 1 );
   sbuf_0 = -1;
   sbuf_1 = 0;
   proc1 = -1;
   for ( i = ilower; i <= iupper; i++ )
   {
      /* Each send buffer is [ num indices in message, index, index,...] */
      _braid_GetBlockDistProc( core, gupper +1, i , &proc );      
      if ( proc != proc1 )
      {
         sbuf_0++;
         sbuf_1 = 2;
         send_buffer[sbuf_0] = _braid_CTAlloc( braid_Int, num_sends + 1 );
         send_buffer[sbuf_0][1] = i ;
         send_proc[sbuf_0] = proc;
         proc1 = proc;
      }
      else
      {
         send_buffer[sbuf_0][sbuf_1] = i;
         sbuf_1++;
      }
      send_buffer[sbuf_0][0]++;
   }

   /* Send the buffers to the processors. */
   send_request = _braid_CTAlloc( MPI_Request , sbuf_0 + 1 );
   for ( i = 0; i <= sbuf_0; i++ )
   {
      MPI_Isend( &send_buffer[i][0] , 1 + send_buffer[i][0] , braid_MPI_INT, send_proc[i] , 33, comm, &send_request[i] );
   }

   /* Wait for messages from other processors claiming to own points in my assumed partition */
   num_messages = *assumed_iupper - *assumed_ilower + 1 ;
   recv_buffer = _braid_CTAlloc( braid_Int, num_messages + 1 );
   while ( num_messages > 0 )
   {
      MPI_Recv( recv_buffer, num_messages + 1, braid_MPI_INT, MPI_ANY_SOURCE, 33, comm, &recv_status);
      for ( i = 0 ; i < recv_buffer[0] ; i++ )
      {
         /* Each Recv buffer is [ num indices in message, index, index, ...] */
         owners[ recv_buffer[i+1] - *assumed_ilower] = recv_status.MPI_SOURCE;
      }
      num_messages -= recv_buffer[0];
   }

   /* Finish up the communication */
   if (sbuf_0 >= 0  )
   {
      MPI_Waitall( sbuf_0 + 1, send_request, MPI_STATUS_IGNORE );
   }

   /* Free up some memory */
   for ( i = 0; i <= sbuf_0; i++ )
   {
      _braid_TFree( send_buffer[i] );
   }
   _braid_TFree ( send_buffer );
   _braid_TFree( send_proc );
   _braid_TFree( send_request );
   _braid_TFree( recv_buffer );

   /* return */
   *owner = owners;
   return _braid_error_flag;
}


/*--------------------------------------------------------------------------------------
 * Used in the weighted data distribution all reduce call, this function combines all 
 * neccesary reductions into a single call.  
 *------------------------------------------------------------------------------------*/

void Sum3MaxMin(braid_Real *in, braid_Real *inout, braid_Int *len, MPI_Datatype *datatype)
{
    /* Sum, Sum, Max, Min */
    if ( *len != 5 )
    {
       abort();
    }

    inout[0] = in[0] + inout[0];
    inout[1] = _braid_max( in[1] , inout[1] );
    inout[2] = _braid_min( in[2] , inout[2] );
    inout[3] = _braid_max( in[3] , inout[3] );
    inout[4] = _braid_min( in[4] , inout[4] );

    return;
}
