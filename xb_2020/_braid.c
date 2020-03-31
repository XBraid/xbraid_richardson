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

/** \file _braid.c
 * \brief Source code for developer routines.  See braid.h for more information.
 *
 */

#include "_braid.h"
#include "braid_defs.h"
#include "_util.h"
#include "_braid_status.h"
#include "_braid_dist.h"

#define DEBUG 0

#if DEBUG
braid_Int  FRefine_count = 0;
#endif

braid_Int _braid_error_flag = 0;
FILE    *_braid_printfile  = NULL;



/*----------------------------------------------------------------------------
 * Returns the index interval for my processor on the finest grid level using 
 * a block data distribution. 
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetInitDistribution(braid_Core   core,
                           braid_Int   *ilower_ptr,
                           braid_Int   *iupper_ptr)
{
   MPI_Comm   comm    = _braid_CoreElt(core, comm);
   braid_Int  gupper = _braid_CoreElt(core, gupper);
   braid_Int  npoints, nprocs, proc;

   npoints = gupper + 1;
   MPI_Comm_size(comm, &nprocs);
   MPI_Comm_rank(comm, &proc);

   _braid_GetBlockDistInterval(core, npoints, ilower_ptr, iupper_ptr);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Sets the neighbouring processors for the initial distribution assuming the 
 * initial grid uses a block data distribution. 
 */

braid_Int
_braid_SetInitNeighbours(braid_Core     core,
                         braid_Int      nlevels,
                         braid_Int      gupper,
                         _braid_Grid  **grids)
{
   braid_Int i, level, cfactor, comm_size, proc;
   braid_Int ilower, iupper;
   _braid_Grid *grid;

   MPI_Comm comm = _braid_CoreElt( core, comm );
   MPI_Comm_size( comm, &comm_size );

   for ( level = 0; level < nlevels; level++ )
   {
      grid = grids[level];
      ilower = _braid_GridElt(grid, ilower);
      iupper = _braid_GridElt(grid, iupper);
      if ( ilower <= iupper )
      {
         braid_Int index = ilower - 1;
         for ( i = level-1; i > -1; i-- )
         {
            _braid_GetCFactor( core, i, &cfactor );
            _braid_MapCoarseToFine( index, cfactor, index );
         }
         _braid_GetBlockDistProc( core, gupper + 1 , index, &proc );
         _braid_GridElt( grids[level] , left_proc ) = proc;

         index = iupper + 1;
         for ( i = level-1; i > -1; i-- )
         {
            _braid_GetCFactor( core, i, &cfactor );
            _braid_MapCoarseToFine( index, cfactor, index );
         }
         _braid_GetBlockDistProc( core, gupper + 1 , index, &proc );
         _braid_GridElt( grids[level] , right_proc ) = proc;
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns the processor that owns the index to the left ( direction = -1 )
 * or the right ( direction = 1 ) on level l.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetProcLeftOrRight(braid_Core   core,
                          braid_Int    level,
                          braid_Int    direction,
                          braid_Int   *proc_ptr)
{
   _braid_Grid  **grids = _braid_CoreElt( core, grids );
   if ( direction < 0 )
   {
      *proc_ptr = _braid_GridElt( grids[level] , left_proc );
   }
   else
   {
      *proc_ptr = _braid_GridElt( grids[level] , right_proc );
   }
   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns the coarsening factor to use on grid 'level'
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetCFactor(braid_Core   core,
                  braid_Int    level,
                  braid_Int   *cfactor_ptr)
{
   braid_Int     *cfactors  = _braid_CoreElt(core, cfactors);
   braid_Int      cfdefault = _braid_CoreElt(core, cfdefault);
   braid_Int      cfactor;

   if (cfactors[level] != 0)
   {
      cfactor = cfactors[level];
   }
   else
   {
      cfactor = cfdefault;
   }
   *cfactor_ptr = cfactor;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_CommRecvInit(braid_Core           core,
                    braid_Int            level,
                    braid_Int            index,
                    braid_Vector        *vector_ptr,
                    _braid_CommHandle  **handle_ptr)
{
   MPI_Comm            comm = _braid_CoreElt(core, comm);
   braid_App           app  = _braid_CoreElt(core, app);
   _braid_CommHandle  *handle = NULL;
   void               *buffer;
   MPI_Request        *requests;
   MPI_Status         *status;
   braid_Int           proc, size, num_requests;
   braid_BufferStatus bstatus = (braid_BufferStatus)core;

   _braid_GetProcLeftOrRight( core, level, -1 , &proc );

   if (proc > -1)
   {
      handle = _braid_TAlloc(_braid_CommHandle, 1);

      /* Allocate buffer through user routine */
      _braid_BufferStatusInit( 0, 0, bstatus );
      _braid_CoreFcn(core, bufsize)(app, &size, bstatus);
      buffer = malloc(size);

      num_requests = 1;
      requests = _braid_CTAlloc(MPI_Request, num_requests);
      status   = _braid_CTAlloc(MPI_Status, num_requests);
      MPI_Irecv(buffer, size, MPI_BYTE, proc, 0, comm, &requests[0]);

      _braid_CommHandleElt(handle, request_type) = 1; /* recv type = 1 */
      _braid_CommHandleElt(handle, num_requests) = num_requests;
      _braid_CommHandleElt(handle, requests)     = requests;
      _braid_CommHandleElt(handle, status)       = status;
      _braid_CommHandleElt(handle, buffer)       = buffer;
      _braid_CommHandleElt(handle, vector_ptr)   = vector_ptr;
   }

   *handle_ptr = handle;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_CommSendInit(braid_Core           core,
                    braid_Int            level,
                    braid_Int            index,
                    braid_Vector         vector,
                    _braid_CommHandle  **handle_ptr)
{
   MPI_Comm            comm = _braid_CoreElt(core, comm);
   braid_App           app  = _braid_CoreElt(core, app);
   _braid_CommHandle  *handle = NULL;
   void               *buffer;
   MPI_Request        *requests;
   MPI_Status         *status;
   braid_Int           proc, size, num_requests;
   braid_BufferStatus  bstatus   = (braid_BufferStatus)core;
   

   _braid_GetProcLeftOrRight(core, level, 1, &proc);
   if (proc > -1)
   {
      handle = _braid_TAlloc(_braid_CommHandle, 1);

      /* Allocate buffer through user routine */
      _braid_BufferStatusInit( 0, 0, bstatus );
      _braid_CoreFcn(core, bufsize)(app, &size, bstatus);
      buffer = malloc(size);

      /* Note that bufpack may return a size smaller than bufsize */ 
      _braid_StatusElt(bstatus, size_buffer) = size;
      _braid_CoreFcn(core, bufpack)(app, vector, buffer, bstatus);
      size = _braid_StatusElt( bstatus, size_buffer );

      num_requests = 1;
      requests = _braid_CTAlloc(MPI_Request, num_requests);
      status   = _braid_CTAlloc(MPI_Status, num_requests);
      MPI_Isend(buffer, size, MPI_BYTE, proc, 0, comm, &requests[0]);

      _braid_CommHandleElt(handle, request_type) = 0; /* send type = 0 */
      _braid_CommHandleElt(handle, num_requests) = num_requests;
      _braid_CommHandleElt(handle, requests)     = requests;
      _braid_CommHandleElt(handle, status)       = status;
      _braid_CommHandleElt(handle, buffer)       = buffer;
   }

   *handle_ptr = handle;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_CommWait(braid_Core          core,
                _braid_CommHandle **handle_ptr)
{
   braid_App           app    = _braid_CoreElt(core, app);
   _braid_CommHandle  *handle = *handle_ptr;

   if (handle != NULL)
   {
      braid_Int      request_type = _braid_CommHandleElt(handle, request_type);
      braid_Int      num_requests = _braid_CommHandleElt(handle, num_requests);
      MPI_Request   *requests     = _braid_CommHandleElt(handle, requests);
      MPI_Status    *status       = _braid_CommHandleElt(handle, status);
      void          *buffer       = _braid_CommHandleElt(handle, buffer);
      braid_BufferStatus bstatus  = (braid_BufferStatus)core;

      MPI_Waitall(num_requests, requests, status);
      
      if (request_type == 1) /* recv type */
      {
         _braid_BufferStatusInit( 0, 0, bstatus );
         braid_Vector  *vector_ptr = _braid_CommHandleElt(handle, vector_ptr);
         _braid_CoreFcn(core, bufunpack)(app, buffer, vector_ptr, bstatus);
      }
      
      _braid_TFree(requests);
      _braid_TFree(status);
      _braid_TFree(handle);
      _braid_TFree(buffer);

      *handle_ptr = NULL;
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns an index into the local u-vector for grid 'level' at point 'index'.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UGetIndex(braid_Core   core,
                 braid_Int    level,
                 braid_Int    index,
                 braid_Int   *uindex_ptr,
                 braid_Int   *store_flag_ptr)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Int            ilower      = _braid_GridElt(grids[level], ilower);
   braid_Int            iupper      = _braid_GridElt(grids[level], iupper);
   braid_Int            clower      = _braid_GridElt(grids[level], clower);
   braid_Int            cfactor     = _braid_GridElt(grids[level], cfactor);
   braid_Int            uindex, ic, iclo, store_flag;

   uindex = -1;
   store_flag = -2;
   if ((index >= ilower) && (index <= iupper))
   {
      if ( _braid_CoreElt(core, useshell) == 1)
      {
         uindex = index-ilower;
         store_flag = 0;
         // If we are not on a fully-stored point
         // then we only have a shell, the store_flag should be -1
         if ( (_braid_CoreElt(core, storage) < 0) ||
              (level < _braid_CoreElt(core, storage)) )
         {
            if ( !_braid_IsCPoint(index, cfactor) )
            {
               store_flag = -1;
            }
         }
      }
      else
      {
         // If on level that only stores C-points
         if ( (_braid_CoreElt(core, storage) < 0) ||
              (level < _braid_CoreElt(core, storage)) )
         {
            if ( _braid_IsCPoint(index, cfactor) )
            {
               _braid_MapFineToCoarse(index, cfactor, ic);
               _braid_MapFineToCoarse(clower, cfactor, iclo);
               uindex = ic-iclo;
               store_flag = 0;
            }
         }
         else
         {
            uindex = index-ilower;
            store_flag = 0;
         }
      }
   }

   *uindex_ptr = uindex;
   *store_flag_ptr = store_flag;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns a reference to the local u-vector on grid 'level' at point 'index'.
 * If the u-vector is not stored, returns NULL. The referred u-vector might
 * just be a shell if that feature is used.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UGetVectorRef(braid_Core     core,
                     braid_Int      level,
                     braid_Int      index,
                     braid_Vector  *u_ptr)
{
   _braid_Grid        **grids = _braid_CoreElt(core, grids);
   braid_Vector        *ua    = _braid_GridElt(grids[level], ua);
   braid_Int            iu, sflag;
   braid_Vector         u = NULL;

   _braid_UGetIndex(core, level, index, &iu, &sflag);
   if (sflag>-2) // We have a full point or a shell (iu>=0)
   {
      u = ua[iu];
   }

   *u_ptr = u;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Stores a reference to the local u-vector on grid 'level' at point 'index'.
 * If the shellvector feature is used, the u-vector might be emptied so that
 * only the shell is stored. Otherwise, if the u-vector is not stored, nothing
 * is done.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_USetVectorRef(braid_Core    core,
                     braid_Int     level,
                     braid_Int     index,
                     braid_Vector  u)
{
   _braid_Grid        **grids = _braid_CoreElt(core, grids);
   braid_Vector        *ua    = _braid_GridElt(grids[level], ua);
   braid_Int            iu, sflag;

   _braid_UGetIndex(core, level, index, &iu, &sflag);
   // If sflag ==0, we have a full point, if sflag == -1, we have a shell
   if (sflag == 0)
   {
      ua[iu] = u;
   }
   else if (sflag == -1)
   {
      braid_App    app = _braid_CoreElt(core, app);
      _braid_CoreFcn(core, sfree)(app, u);
      ua[iu] = u;
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns a copy of the u-vector on grid 'level' at point 'index'.  If 'index'
 * is my "receive index" (as set by UCommInit(), for example), the u-vector will
 * be received from a neighbor processor.  If the u-vector is not stored, NULL
 * is returned. The copy might just be a shell if this feature is used.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UGetVector(braid_Core     core,
                  braid_Int      level,
                  braid_Int      index,
                  braid_Vector  *u_ptr)
{
   braid_App            app         = _braid_CoreElt(core, app);
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Vector        *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            recv_index  = _braid_GridElt(grids[level], recv_index);
   _braid_CommHandle   *recv_handle = _braid_GridElt(grids[level], recv_handle);
   braid_Vector         u = NULL;
   braid_Int            iu, sflag;

   if (index == recv_index)
   {
      /* If a recv was initiated, receive u value from neighbor processor */
      if (recv_index > -1)
      {
         _braid_CommWait(core, &recv_handle);
         _braid_GridElt(grids[level], recv_index)  = -1;
         _braid_GridElt(grids[level], recv_handle) = recv_handle;
         u = ua[-1];
      }
   }
   else
   {
      _braid_UGetIndex(core, level, index, &iu, &sflag);
      if (sflag == 0)
      {
         _braid_CoreFcn(core, clone)(app, ua[iu], &u);
      }
      else if (sflag == -1)
      {
         // In this case, sclone != NULL
         _braid_CoreFcn(core, sclone)(app, ua[iu], &u);
      }
   }

   *u_ptr = u;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Stores the u-vector on grid 'level' at point 'index'.  If 'index' is my "send
 * index", a send is initiated to a neighbor processor.  If 'move' is true, the
 * u-vector is moved into core storage instead of copied.  If the u-vector is
 * not stored, nothing is done or only the shell is copied/moved when the shellvector
 * feature is used.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_USetVector(braid_Core    core,
                  braid_Int     level,
                  braid_Int     index,
                  braid_Vector  u,
                  braid_Int     move)
{
   braid_App            app         = _braid_CoreElt(core, app);
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Vector        *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            send_index  = _braid_GridElt(grids[level], send_index);
   _braid_CommHandle   *send_handle = _braid_GridElt(grids[level], send_handle);
   braid_Int            iu, sflag;

   if (index == send_index)
   {
      /* Post send to neighbor processor */
      _braid_CommSendInit(core, level, index, u, &send_handle);
      _braid_GridElt(grids[level], send_index)  = -1;
      _braid_GridElt(grids[level], send_handle) = send_handle;
   }

   _braid_UGetIndex(core, level, index, &iu, &sflag);
   if (sflag == 0) // We have a full point
   {
      if (ua[iu] != NULL)
      {
         _braid_CoreFcn(core, free)(app, ua[iu]);
      }
      if (move)
      {
         ua[iu] = u;                                   /* move the vector */
      }
      else
      {
         _braid_CoreFcn(core, clone)(app, u, &ua[iu]); /* copy the vector */
      }
   }
   else if (sflag == -1) // We have a shell
   {
      if (ua[iu] != NULL)
      {
         _braid_CoreFcn(core, free)(app, ua[iu]);
      }
      if (move)
      {
         // We are on an F-point, with shellvector option. We only keep the shell.
         _braid_CoreFcn(core, sfree)(app, u);
         ua[iu] = u;                                   /* move the vector */
      }
      else
      {
         _braid_CoreFcn(core, sclone)(app, u, &ua[iu]); /* copy the vector */
      }
   }
   else if (move) // We store nothing
   {
      _braid_CoreFcn(core, free)(app, u);              /* free the vector */
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Basic communication (from the left, to the right)
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UCommInitBasic(braid_Core  core,
                      braid_Int   level,
                      braid_Int   recv_msg,
                      braid_Int   send_msg,
                      braid_Int   send_now)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Int            ilower      = _braid_GridElt(grids[level], ilower);
   braid_Int            iupper      = _braid_GridElt(grids[level], iupper);
   braid_Vector        *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            recv_index  = -1;
   braid_Int            send_index  = -1;
   _braid_CommHandle   *recv_handle = NULL;
   _braid_CommHandle   *send_handle = NULL;
   braid_Int            iu, sflag;

   /* Post receive */
   if (recv_msg)
   {
      recv_index = ilower-1;
      _braid_CommRecvInit(core, level, recv_index, &ua[-1], &recv_handle);
   }

   /* Post send */
   if (send_msg)
   {
      send_index = iupper;
      if (send_now)
      {
         _braid_UGetIndex(core, level, send_index, &iu, &sflag);
         if (sflag < 0)
         {
            // We should never get here : we do not communicate shells...
            abort();
         }
         _braid_CommSendInit(core, level, send_index, ua[iu], &send_handle);
         send_index = -1;
      }
   }

   _braid_GridElt(grids[level], recv_index)  = recv_index ;
   _braid_GridElt(grids[level], send_index)  = send_index ;
   _braid_GridElt(grids[level], recv_handle) = recv_handle;
   _braid_GridElt(grids[level], send_handle) = send_handle;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Working on all intervals
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UCommInit(braid_Core  core,
                 braid_Int   level)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Int            ilower      = _braid_GridElt(grids[level], ilower);
   braid_Int            iupper      = _braid_GridElt(grids[level], iupper);
   braid_Int            cfactor     = _braid_GridElt(grids[level], cfactor);
   braid_Vector        *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            recv_index  = -1;
   braid_Int            send_index  = -1;
   _braid_CommHandle   *recv_handle = NULL;
   _braid_CommHandle   *send_handle = NULL;
   braid_Int            iu, sflag;
   
   /* Note that this routine works for the case of all points being C-points, 
    * i.e., cfactor = 1.  A send and receive are always posted. */

   if (ilower <= iupper)
   {
      /* Post receive */
      _braid_CommRecvInit(core, level, ilower-1, &ua[-1], &recv_handle);
      recv_index = ilower-1;
      
      /* Only post send if iupper is a C-point, otherwise compute and send later */
      if ( _braid_IsCPoint(iupper, cfactor) )
      {
         _braid_UGetIndex(core, level, iupper, &iu, &sflag);
         _braid_CommSendInit(core, level, iupper, ua[iu], &send_handle);
         send_index = -1;
      }
      else
      {
         send_index = iupper;
      }
   }

   _braid_GridElt(grids[level], recv_index)  = recv_index ;
   _braid_GridElt(grids[level], send_index)  = send_index ;
   _braid_GridElt(grids[level], recv_handle) = recv_handle;
   _braid_GridElt(grids[level], send_handle) = send_handle;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Working only on F-pt intervals
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UCommInitF(braid_Core  core,
                  braid_Int   level)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Int            ilower      = _braid_GridElt(grids[level], ilower);
   braid_Int            iupper      = _braid_GridElt(grids[level], iupper);
   braid_Int            cfactor     = _braid_GridElt(grids[level], cfactor);
   braid_Vector        *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            recv_index  = -1;
   braid_Int            send_index  = -1;
   _braid_CommHandle   *recv_handle = NULL;
   _braid_CommHandle   *send_handle = NULL;
   braid_Int            iu, sflag;

   if (ilower <= iupper)
   {
      /* Only post receive if ilower is an F-point */
      if ( _braid_IsFPoint(ilower, cfactor) )
      {
         _braid_CommRecvInit(core, level, ilower-1, &ua[-1], &recv_handle);
         recv_index = ilower-1;
      }

      /* Only post send if iupper is a C-point and iupper+1 is an F-point.  This
       * check allows for the case of cfactor=1, i.e., all C-points.  Otherwise, 
       * if iupper+1 is an F-point, set send_index, so that when that point is 
       * computed later, it is sent. */
      if ( _braid_IsCPoint(iupper, cfactor) && _braid_IsFPoint(iupper+1, cfactor))
      {
         _braid_UGetIndex(core, level, iupper, &iu, &sflag);
         _braid_CommSendInit(core, level, iupper, ua[iu], &send_handle);
         send_index = -1;
      }
      else if ( _braid_IsFPoint(iupper+1, cfactor) )
      {
         send_index = iupper;
      }
   }

   _braid_GridElt(grids[level], recv_index)  = recv_index ;
   _braid_GridElt(grids[level], send_index)  = send_index ;
   _braid_GridElt(grids[level], recv_handle) = recv_handle;
   _braid_GridElt(grids[level], send_handle) = send_handle;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Finish up communication
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UCommWait(braid_Core  core,
                 braid_Int   level)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   _braid_CommHandle   *recv_handle = _braid_GridElt(grids[level], recv_handle);
   _braid_CommHandle   *send_handle = _braid_GridElt(grids[level], send_handle);

   _braid_CommWait(core, &recv_handle);
   _braid_CommWait(core, &send_handle);
   _braid_GridElt(grids[level], recv_index)  = -1;
   _braid_GridElt(grids[level], send_index)  = -1;
   _braid_GridElt(grids[level], recv_handle) = recv_handle;
   _braid_GridElt(grids[level], send_handle) = send_handle;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetInterval(braid_Core   core,
                   braid_Int    level,
                   braid_Int    interval_index,
                   braid_Int   *flo_ptr,
                   braid_Int   *fhi_ptr,
                   braid_Int   *ci_ptr)
{
   _braid_Grid  **grids   = _braid_CoreElt(core, grids);
   braid_Int      ilower  = _braid_GridElt(grids[level], ilower);
   braid_Int      iupper  = _braid_GridElt(grids[level], iupper);
   braid_Int      clower  = _braid_GridElt(grids[level], clower);
   braid_Int      cupper  = _braid_GridElt(grids[level], cupper);
   braid_Int      cfactor = _braid_GridElt(grids[level], cfactor);
   braid_Int      flo, fhi, ci;

   flo = ilower;
   fhi = iupper;
   ci  = -1;

   if ( _braid_IsCPoint(clower, cfactor) )
   {
      flo = clower + (interval_index-1)*cfactor + 1;
      fhi = clower + (interval_index  )*cfactor - 1;
      if (flo < ilower)
      {
         flo = ilower;
      }
      if (fhi > iupper)
      {
         fhi = iupper;
      }

      ci = clower + interval_index*cfactor;
      if (ci > cupper)
      {
         ci = -1;  /* return -1 if no C-points */
      }
   }

   *flo_ptr = flo;
   *fhi_ptr = fhi;
   *ci_ptr  = ci;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_AccessVector(braid_Core          core,
                    braid_AccessStatus  status,
                    braid_Vector        u)
{
   braid_App      app    = _braid_CoreElt(core, app);

   _braid_CoreFcn(core, access)(app, u, status);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Get an initial guess for ustop to use in the step routine (implicit schemes)
 * This vector may just be a shell. User should be able to deal with it
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetUInit(braid_Core     core,
                braid_Int      level,
                braid_Int      index,
                braid_Vector   u,
                braid_Vector  *ustop_ptr)
{
   _braid_Grid    **grids    = _braid_CoreElt(core, grids);
   braid_Int        ilower   = _braid_GridElt(grids[level], ilower);
   braid_Int        storage  = _braid_CoreElt(core, storage);
   braid_Vector    *va       = _braid_GridElt(grids[level], va);

   braid_Vector     ustop = *ustop_ptr;
   braid_Int        ii;

   ii = index-ilower;

   _braid_UGetVectorRef(core, level, index, &ustop);

   /* If ustop is NULL, then storage is only at C-points on this level and this
    * is an F-point.  See the comment block around FRestrict() for the fixed-point
    * logic behind our choices in ustop. */
   if( ustop == NULL)
   {
      if( (level == 0) || ( storage == -2 ) )
      {
         ustop = u;
      }
      else
      {
         ustop = va[ii];
      }
   }

   /* If you have storage at this point, use it, unless you're in compatibility mode (-2). */
   else if( storage == -2 )
   {
      if ( _braid_CoreElt(core, useshell) == 1)
      {
         // Should not happen, ustop is never NULL with useshell option
         // unless there are inconsistent options (i.e. useshell && storage==-2)
         abort();
      }
      ustop = u;
   }

   *ustop_ptr = ustop;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Integrate one time step
 *----------------------------------------------------------------------------*/

braid_Int
_braid_Step(braid_Core     core,
            braid_Int      level,
            braid_Int      index,
            braid_Vector   ustop,
            braid_Vector   u)
{
   braid_App        app      = _braid_CoreElt(core, app);
   braid_Real       tol      = _braid_CoreElt(core, tol);
   braid_Int        iter     = _braid_CoreElt(core, niter);
   braid_Int       *rfactors = _braid_CoreElt(core, rfactors);
   _braid_Grid    **grids    = _braid_CoreElt(core, grids);
   braid_StepStatus status   = (braid_StepStatus)core;
   braid_Int        nrefine  = _braid_CoreElt(core, nrefine);
   braid_Int        gupper   = _braid_CoreElt(core, gupper);
   braid_Int        ilower   = _braid_GridElt(grids[level], ilower);
   braid_Real      *ta       = _braid_GridElt(grids[level], ta);
   braid_Vector    *fa       = _braid_GridElt(grids[level], fa);

   braid_Int        ii;

   ii = index-ilower;
   _braid_StepStatusInit(ta[ii-1], ta[ii], index-1, tol, iter, level, nrefine, gupper, status);

   /* If ustop is set to NULL, use a default approach for setting it */
   if (ustop == NULL)
   {
      _braid_GetUInit(core, level, index, u, &ustop);
   }

   if (level == 0)
   {
      _braid_CoreFcn(core, step)(app, ustop, NULL, u, status);
      rfactors[ii] = _braid_StatusElt(status, rfactor); 
   }
   else
   {
      if ( _braid_CoreElt(core, residual) == NULL )
      {
         _braid_CoreFcn(core, step)(app, ustop, NULL, u, status);
         if(fa[ii] != NULL)
         {
            _braid_CoreFcn(core, sum)(app, 1.0, fa[ii], 1.0, u);
         }
      }
      else
      {
         _braid_CoreFcn(core, step)(app, ustop, fa[ii], u, status);
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Compute residual
 *----------------------------------------------------------------------------*/

braid_Int
_braid_Residual(braid_Core     core,
                braid_Int      level,
                braid_Int      index,
                braid_Vector   ustop,
                braid_Vector   r)
{
   braid_App        app      = _braid_CoreElt(core, app);
   braid_Real       tol      = _braid_CoreElt(core, tol);
   braid_Int        iter     = _braid_CoreElt(core, niter);
   braid_Int       *rfactors = _braid_CoreElt(core, rfactors);
   _braid_Grid    **grids    = _braid_CoreElt(core, grids);
   braid_StepStatus status   = (braid_StepStatus)core;
   braid_Int        nrefine  = _braid_CoreElt(core, nrefine);
   braid_Int        gupper   = _braid_CoreElt(core, gupper);
   braid_Int        ilower   = _braid_GridElt(grids[level], ilower);
   braid_Real      *ta       = _braid_GridElt(grids[level], ta);

   braid_Vector     rstop;
   braid_Int        ii;

   ii = index-ilower;
   _braid_StepStatusInit(ta[ii-1], ta[ii], index-1, tol, iter, level, nrefine, gupper, status);
   if ( _braid_CoreElt(core, residual) == NULL )
   {
      /* By default: r = ustop - \Phi(ustart)*/
      _braid_GetUInit(core, level, index, r, &rstop);
      _braid_CoreFcn(core, step)(app, rstop, NULL, r, status);
      _braid_CoreFcn(core, sum)(app, 1.0, ustop, -1.0, r);
      if (level == 0)
      {
         /*TODO Remove this line after modifing the _braid_StatusSetRFactor to set the rfactor in the array directly */
         rfactors[ii] = _braid_StatusElt(status, rfactor);
      }
   }
   else
   {
      /* Call the user's residual routine */
      _braid_CoreFcn(core, residual)(app, ustop, r, status);
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Compute FAS residual = f - residual
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FASResidual(braid_Core     core,
                   braid_Int      level,
                   braid_Int      index,
                   braid_Vector   ustop,
                   braid_Vector   r)
{
   braid_App        app    = _braid_CoreElt(core, app);
   _braid_Grid    **grids  = _braid_CoreElt(core, grids);
   braid_Int        ilower = _braid_GridElt(grids[level], ilower);
   braid_Vector    *fa     = _braid_GridElt(grids[level], fa);

   braid_Int        ii;

   _braid_Residual(core, level, index, ustop, r);
   if (level == 0)
   {
      _braid_CoreFcn(core, sum)(app, 0.0, r, -1.0, r);
   }
   else
   {
      ii = index-ilower;
      if(fa[ii] == NULL)
      {
         _braid_CoreFcn(core, sum)(app, 0.0, r, -1.0, r);
      }
      else
      {
         _braid_CoreFcn(core, sum)(app, 1.0, fa[ii], -1.0, r);
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Coarsen in space
 *----------------------------------------------------------------------------*/

braid_Int
_braid_Coarsen(braid_Core     core,
               braid_Int      level,    /* coarse level */
               braid_Int      f_index,  /* fine index */
               braid_Int      c_index,  /* coarse index */
               braid_Vector   fvector,
               braid_Vector  *cvector)
{
   braid_App      app             = _braid_CoreElt(core, app);
   _braid_Grid  **grids           = _braid_CoreElt(core, grids);
   braid_CoarsenRefStatus cstatus = (braid_CoarsenRefStatus)core;
   braid_Int      nrefine         = _braid_CoreElt(core, nrefine);
   braid_Int      gupper          = _braid_CoreElt(core, gupper);
   braid_Int      c_ilower        = _braid_GridElt(grids[level], ilower);
   braid_Int      f_ilower        = _braid_GridElt(grids[level-1], ilower);
   braid_Real    *c_ta            = _braid_GridElt(grids[level], ta);
   braid_Real    *f_ta            = _braid_GridElt(grids[level-1], ta);

   braid_Int      c_ii = c_index-c_ilower;
   braid_Int      f_ii = f_index-f_ilower;
   
   if ( _braid_CoreElt(core, scoarsen) == NULL )
   {
      /* No spatial coarsening needed, just clone the fine vector.*/
      _braid_CoreFcn(core, clone)(app, fvector, cvector);
   }
   else
   {
      /* Call the user's coarsening routine */
      _braid_CoarsenRefStatusInit(f_ta[f_ii], f_ta[f_ii-1], f_ta[f_ii+1], 
                                  c_ta[c_ii-1], c_ta[c_ii+1],
                                  level-1, nrefine, gupper, cstatus);
      _braid_CoreFcn(core, scoarsen)(app, fvector, cvector, cstatus);
   }
   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Refine in space (basic routine)
 *----------------------------------------------------------------------------*/

braid_Int
_braid_RefineBasic(braid_Core     core,
                   braid_Int      level,    /* fine level */
                   braid_Real    *f_ta,     /* pointer into fine time array */
                   braid_Real    *c_ta,     /* pointer into coarse time array */
                   braid_Vector   cvector,
                   braid_Vector  *fvector)
{
   braid_App              app     = _braid_CoreElt(core, app);
   braid_CoarsenRefStatus cstatus = (braid_CoarsenRefStatus)core;
   braid_Int              nrefine = _braid_CoreElt(core, nrefine);
   braid_Int              gupper  = _braid_CoreElt(core, gupper);

   if ( _braid_CoreElt(core, scoarsen) == NULL )
   {
      /* No spatial refinement needed, just clone the fine vector.*/
      _braid_CoreFcn(core, clone)(app, cvector, fvector);
   }
   else
   {
      /* Call the user's refinement routine */
      _braid_CoarsenRefStatusInit(f_ta[0], f_ta[-1], f_ta[+1], c_ta[-1], c_ta[+1],
                                  level, nrefine, gupper, cstatus);
      _braid_CoreFcn(core, srefine)(app, cvector, fvector, cstatus);
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Refine in space
 *----------------------------------------------------------------------------*/

braid_Int
_braid_Refine(braid_Core     core,
              braid_Int      level,    /* fine level */
              braid_Int      f_index,  /* fine index */
              braid_Int      c_index,  /* coarse index */
              braid_Vector   cvector,
              braid_Vector  *fvector)
{
   _braid_Grid  **grids           = _braid_CoreElt(core, grids);
   braid_Int      f_ilower        = _braid_GridElt(grids[level], ilower);
   braid_Int      c_ilower        = _braid_GridElt(grids[level+1], ilower);
   braid_Real    *f_ta            = _braid_GridElt(grids[level], ta);
   braid_Real    *c_ta            = _braid_GridElt(grids[level+1], ta);

   braid_Int      c_ii = c_index-c_ilower;
   braid_Int      f_ii = f_index-f_ilower;

   _braid_RefineBasic(core, level, &f_ta[f_ii], &c_ta[c_ii], cvector, fvector);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Create a new grid object
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GridInit(braid_Core     core,
                braid_Int      level,
                braid_Int      ilower,
                braid_Int      iupper,
                _braid_Grid  **grid_ptr)
{
   _braid_Grid   *grid;
   braid_Real    *ta;

   grid = _braid_CTAlloc(_braid_Grid, 1);
   
   _braid_GridElt(grid, level)  = level;
   _braid_GridElt(grid, ilower) = ilower;
   _braid_GridElt(grid, iupper) = iupper;
   _braid_GridElt(grid, recv_index) = -1;
   _braid_GridElt(grid, send_index) = -1;
   _braid_GridElt(grid, left_proc)  = -1;
   _braid_GridElt(grid, right_proc)  = -1;
   
   /* Store each processor's time slice, plus one time value to the left 
    * and to the right */
   ta = _braid_CTAlloc(braid_Real, iupper-ilower+3);
   _braid_GridElt(grid, ta_alloc) = ta;
   _braid_GridElt(grid, ta)       = ta+1;  /* shift */
   
   *grid_ptr = grid;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GridClean(braid_Core    core,
                 _braid_Grid  *grid)
{
   braid_App      app      = _braid_CoreElt(core, app);
   braid_Int      ilower   = _braid_GridElt(grid, ilower);
   braid_Int      iupper   = _braid_GridElt(grid, iupper);
   braid_Int      nupoints = _braid_GridElt(grid, nupoints);
   braid_Vector  *ua       = _braid_GridElt(grid, ua);
   braid_Vector  *va       = _braid_GridElt(grid, va);
   braid_Vector  *fa       = _braid_GridElt(grid, fa);
   braid_Vector  *ua_alloc = _braid_GridElt(grid, ua_alloc);
   braid_Vector  *va_alloc = _braid_GridElt(grid, va_alloc);
   braid_Vector  *fa_alloc = _braid_GridElt(grid, fa_alloc);
   
   braid_Int      ii;

   if (ua_alloc)
   {
      for (ii = 0; ii < nupoints; ii++)
      {
         if (ua[ii] != NULL)
         {
            _braid_CoreFcn(core, free)(app, ua[ii]);
            ua[ii] = NULL;
         }
      }
   }
   if (va_alloc)
   {
      for (ii = -1; ii <= (iupper-ilower); ii++)
      {
         if (va[ii] != NULL)
         {
            _braid_CoreFcn(core, free)(app, va[ii]);
            va[ii] = NULL;
         }
      }
   }
   if (fa_alloc)
   {
      for (ii = -1; ii <= (iupper-ilower); ii++)
      {
         if (fa[ii] != NULL)
         {
            _braid_CoreFcn(core, free)(app, fa[ii]);
            fa[ii] = NULL;
         }
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GridDestroy(braid_Core    core,
                   _braid_Grid  *grid)
{
   if (grid)
   {
      braid_Vector  *ua_alloc = _braid_GridElt(grid, ua_alloc);
      braid_Real    *ta_alloc = _braid_GridElt(grid, ta_alloc);
      braid_Vector  *va_alloc = _braid_GridElt(grid, va_alloc);
      braid_Vector  *fa_alloc = _braid_GridElt(grid, fa_alloc);

      _braid_GridClean(core, grid);

      if (ua_alloc)
      {
         _braid_TFree(ua_alloc);
      }
      if (ta_alloc)
      {
         _braid_TFree(ta_alloc);
      }
      if (va_alloc)
      {
         _braid_TFree(va_alloc);
      }
      if (fa_alloc)
      {
         _braid_TFree(fa_alloc);
      }

      _braid_TFree(grid);
   }
   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Set initial guess at C-points, and initialize shell at F-points when using
 * shell vectors.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_InitGuess(braid_Core  core,
                 braid_Int   level)
{
   braid_App      app      = _braid_CoreElt(core, app);
   braid_Int      seq_soln = _braid_CoreElt(core, seq_soln);
   _braid_Grid  **grids    = _braid_CoreElt(core, grids);
   braid_Int      ilower   = _braid_GridElt(grids[level], ilower);
   braid_Int      iupper   = _braid_GridElt(grids[level], iupper);
   braid_Int      clower   = _braid_GridElt(grids[level], clower);
   braid_Int      cupper   = _braid_GridElt(grids[level], cupper);
   braid_Int      cfactor  = _braid_GridElt(grids[level], cfactor);
   braid_Real    *ta       = _braid_GridElt(grids[level], ta);
   braid_Vector  *va       = _braid_GridElt(grids[level], va);

   braid_Vector   u;
   braid_Int      i, iu, sflag;

   /* Set the calling function */
   _braid_CoreElt(core, calling_function) = braid_Caller_InitGuess;
   
   if ( (level == 0) && (seq_soln == 1) )
   {
      /* If first processor, grab initial condition */
      if(ilower == 0)
      {
         _braid_CoreFcn(core, init)(app, ta[0], &u);
         _braid_USetVector(core, 0, 0, u, 0);
         ilower += 1;
      }
      /* Else, receive point to the left */
      else
      {
         _braid_UCommInitBasic(core, 0, 1, 0, 0);     /* Post receive to the left*/
         _braid_UGetVector(core, 0, ilower-1, &u);    /* Wait on receive */
      }

      /* Set Flag so that USetVector initiates send when iupper is available */
      _braid_GridElt(grids[level], send_index)  = iupper;

      /* Initialize all points on the finest grid with sequential time marching */
      for(i = ilower; i <= iupper; i++)
      {
         _braid_Step(core, 0, i, NULL, u);       /* Step forward */
         _braid_USetVector(core, 0, i, u, 0);    /* Store: copy u into core,
                                                    sending to left if needed */
      }

      _braid_UCommWait(core, 0);                 /* Wait on comm to finish */
   }
   else if (level == 0)
   {
      if (_braid_CoreElt(core, useshell)==1)
      {
         for (i = ilower; i <=iupper; i++)
         {
            if (_braid_IsCPoint(i,cfactor))
            {
               // We are on a C-point, init full vector
               _braid_CoreFcn(core, init)(app, ta[i-ilower], &u);
            }
            else
            {
               // We are on a F-point, init shell only
               _braid_CoreFcn(core, sinit)(app, ta[i-ilower], &u);
            }
            _braid_USetVectorRef(core, level, i, u);
         }
      }
      else
      {
         /* Only initialize the C-points on the finest grid */
         for (i = clower; i <= cupper; i += cfactor)
         {
            _braid_CoreFcn(core, init)(app, ta[i-ilower], &u);
            _braid_USetVectorRef(core, level, i, u);
         }
      }
   }
   else
   {
      for (i = ilower; i <= iupper; i++)
      {
         _braid_UGetIndex(core, level, i, &iu, &sflag);
         if (sflag == 0) // Full point
         {
            _braid_CoreFcn(core, clone)(app, va[i-ilower], &u);
            _braid_USetVectorRef(core, level, i, u);
         }
         else if (sflag == -1) // Shell
         {
            _braid_CoreFcn(core, sclone)(app, va[i-ilower], &u);
            _braid_USetVectorRef(core, level, i, u);
         }
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_ComputeFullRNorm(braid_Core  core,
                        braid_Int   level,
                        braid_Real *return_rnorm)
{
   MPI_Comm           comm        = _braid_CoreElt(core, comm);
   braid_App          app         = _braid_CoreElt(core, app);
   braid_Real         tol         = _braid_CoreElt(core, tol);
   braid_Int          iter        = _braid_CoreElt(core, niter);
   _braid_Grid      **grids       = _braid_CoreElt(core, grids);
   braid_StepStatus   status      = (braid_StepStatus)core;
   braid_Int          nrefine     = _braid_CoreElt(core, nrefine);
   braid_Int          gupper      = _braid_CoreElt(core, gupper);
   braid_Int          ncpoints    = _braid_GridElt(grids[level], ncpoints);
   braid_Int          tnorm       = _braid_CoreElt(core, tnorm);
   braid_Real        *ta          = _braid_GridElt(grids[level], ta);
   braid_Int          ilower      = _braid_GridElt(grids[level], ilower);
   _braid_CommHandle *send_handle;
   braid_Int          send_index;

   braid_Int        flo, fhi, fi, ci, ii, interval;
   braid_Real       rnorm_temp, rnorm = 0, global_rnorm = 0;
   braid_Vector     u, r;
   
   /* Set the calling function */
   _braid_CoreElt(core, calling_function) = braid_Caller_ComputeFullRNorm;

   _braid_UCommInit(core, level);

   /* Start from the right-most interval. */
   for (interval = ncpoints; interval > -1; interval--)
   {
      _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);
      if (flo <= fhi)
      {
         _braid_UGetVector(core, level, flo-1, &u);
      }
      else if (ci > 0)
      {
         _braid_UGetVector(core, level, ci-1, &u);
      }

      /* Generate F-points and get residual. */
      for (fi = flo; fi <= fhi; fi++)
      {
         _braid_CoreFcn(core, clone)(app, u, &r);
         _braid_Step(core, level, fi, NULL, u);

         /* Update local processor norm. */
         ii = fi-ilower;
         _braid_StepStatusInit(ta[ii-1], ta[ii], fi-1, tol, iter, level, nrefine, gupper, status);
         _braid_CoreFcn(core, full_rnorm_res)(app, u, r, status);
         _braid_CoreFcn(core, spatialnorm)(app, r, &rnorm_temp); 
         if(tnorm == 1)       /* one-norm */ 
         {  
            rnorm += rnorm_temp;
         }
         else if(tnorm == 3)  /* inf-norm */
         {  
            rnorm = (((rnorm_temp) > (rnorm)) ? (rnorm_temp) : (rnorm));
         }
         else                 /* default two-norm */
         {  
            rnorm += (rnorm_temp*rnorm_temp);
         }

         /* Communicate w/ neighbor nodes. */
         send_handle = _braid_GridElt(grids[level], send_handle);
         send_index  = _braid_GridElt(grids[level], send_index);
         if (fi == send_index)
         {
            /* Post send to neighbor processor */
            _braid_CommSendInit(core, level, fi, u, &send_handle);
            _braid_GridElt(grids[level], send_index)  = -1;
            _braid_GridElt(grids[level], send_handle) = send_handle;
         }
         _braid_CoreFcn(core, free)(app, r);
      }
      /* Residual from C-point. */
      if (ci > 0)
      {
         /* Update local processor norm. */
         ii = ci-ilower;
         _braid_StepStatusInit(ta[ii-1], ta[ii], ci-1, tol, iter, level, nrefine, gupper, status);
         _braid_UGetVector(core, level, ci, &r);
         _braid_CoreFcn(core, full_rnorm_res)(app, r, u, status);
         _braid_CoreFcn(core, spatialnorm)(app, u, &rnorm_temp);

         if(tnorm == 1)       /* one-norm */ 
         {  
            rnorm += rnorm_temp;
         }
         else if(tnorm == 3)  /* inf-norm */
         {  
            rnorm = (((rnorm_temp) > (rnorm)) ? (rnorm_temp) : (rnorm));
         }
         else                 /* default two-norm */
         {  
            rnorm += (rnorm_temp*rnorm_temp);
         }
         _braid_CoreFcn(core, free)(app, r);
      }

      if ((flo <= fhi) || (ci > 0))
      {
         _braid_CoreFcn(core, free)(app, u);
      }
   }

   _braid_UCommWait(core, level);

   /* Compute global residual norm. */
   if(tnorm == 1)       /* one-norm reduction */
   {  
      MPI_Allreduce(&rnorm, &global_rnorm, 1, MPI_DOUBLE, MPI_SUM, comm);
   }
   else if(tnorm == 3)  /* inf-norm reduction */
   {  
      MPI_Allreduce(&rnorm, &global_rnorm, 1, MPI_DOUBLE, MPI_MAX, comm);
   }
   else                 /* default two-norm reduction */
   {  
      MPI_Allreduce(&rnorm, &global_rnorm, 1, MPI_DOUBLE, MPI_SUM, comm);
      global_rnorm = sqrt(global_rnorm);
   }

   *return_rnorm = global_rnorm;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Do nu sweeps of F-then-C relaxation
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FCRelax(braid_Core  core,
               braid_Int   level)
{
   braid_App       app      = _braid_CoreElt(core, app);
   braid_Int      *nrels    = _braid_CoreElt(core, nrels);
   _braid_Grid   **grids    = _braid_CoreElt(core, grids);
   braid_Int       nlevels  = _braid_CoreElt(core, nlevels);
   braid_Int       iter     = _braid_CoreElt(core, niter);
   braid_Int       nrefine  = _braid_CoreElt(core, nrefine);
   braid_Int       order    = _braid_CoreElt(core, order);
   braid_Real      tol      = _braid_CoreElt(core, tol );
   braid_Int       ncpoints = _braid_GridElt(grids[level], ncpoints);  
   braid_Int       ilower   = _braid_GridElt(grids[level], ilower);
   braid_Int       iupper   = _braid_GridElt(grids[level], iupper);
   braid_Int       cupper   = _braid_GridElt(grids[level], cupper);
   braid_Int       clower   = _braid_GridElt(grids[level], clower);
   braid_Int       cfactor  = _braid_GridElt(grids[level], cfactor);
   braid_Int       gupper   = _braid_GridElt(grids[level], gupper);
   braid_Real     *ta       = _braid_GridElt(grids[level], ta);
   braid_Real     *dtk      = _braid_CoreElt(core, dtk);
   braid_Vector   *ua       = _braid_GridElt(grids[level], ua);
   braid_Real     *ta_c, time_left, a, b;   
   MPI_Comm        comm     = _braid_CoreElt(core, comm);
   
   MPI_Request     send_request, recv_request;
   void *recv_buff, *send_buff;
   braid_BufferStatus bstatus = (braid_BufferStatus)core;
   braid_StepStatus status = (braid_StepStatus)core;
   braid_Int send_flag;  

   braid_Vector    u, bigstep;
   braid_Int       flo, fhi, fi, ci, size, proc;
   braid_Int       nu, nrelax, interval;
   
   /* Set the calling function */
   _braid_CoreElt(core, calling_function) = braid_Caller_FCRelax;
   
   braid_Int richardson = _braid_CoreElt(core, richardson);
   if ( level > 0 || ncpoints == 0 || nlevels <= 1 || iter + nrefine == 0 )
   {
      /* In this case, nothing needs to be done regarding richardson */
      richardson = 0;
   }

   nrelax  = nrels[level];
   
   send_flag = 0;
   for (nu = 0; nu < nrelax; nu++)
   {
      _braid_UCommInit(core, level);
      braid_Int dtk_index = ncpoints;

      if ( richardson )
      {
         if (  ilower > 0 )
         {
            //Need to post a recv for a cpoint from the left. 
            
            _braid_GetProcLeftOrRight(core, 1, -1, &proc);
            _braid_BufferStatusInit( 0, 0, bstatus );
            _braid_CoreFcn(core, bufsize)(app, &size, bstatus);
             recv_buff = malloc(size);
             
             MPI_Irecv(recv_buff, size, MPI_BYTE, proc, 84, comm, &recv_request);
         }
         if ( iupper < gupper && cupper + cfactor <= gupper   )            
         {
            //Need to post a send of ciupper         
            _braid_GetProcLeftOrRight(core, 1, 1, &proc);
            _braid_BufferStatusInit( 0 ,0 ,bstatus );
            _braid_CoreFcn(core,bufsize)(app, &size, bstatus);
            send_buff = malloc(size);
           
            braid_Int iu, is_stored;
            _braid_UGetIndex(core, level, cupper, &iu, &is_stored);
            _braid_CoreFcn(core, bufpack)(app, ua[iu], send_buff, bstatus);
            size = _braid_StatusElt( bstatus, size_buffer );
            MPI_Isend(send_buff, size, MPI_BYTE, proc, 84, comm, &send_request);
            send_flag = 1;
         }
      }      
      
      /* Start from the right-most interval */
      for (interval = ncpoints; interval > -1; interval--)
      {
         _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);

         if (flo <= fhi)
         {
            _braid_UGetVector(core, level, flo-1, &u);
         }
         else if (ci > 0)
         {
            _braid_UGetVector(core, level, ci-1, &u);
         }
         
         /* Get the solution at previous coarse grid point */
         if ( richardson && ci > 0)  
         {
            if ( ci == clower )
            {
              /* The previous coarse grid solution was sent
               * as a message. Wait and unpack */
               MPI_Wait( &recv_request, MPI_STATUS_IGNORE);
               _braid_BufferStatusInit(0,0,bstatus);
               _braid_CoreFcn(core, bufunpack)(app, recv_buff, &bigstep, bstatus);
               _braid_TFree(recv_buff);                  
              
               ta_c = _braid_GridElt(grids[1], ta );
               time_left = ta_c[-1];
            }
            else
            {
               /* We own the solution at previous coarse grid point, so clone it */
                braid_Int iu, is_stored;
               _braid_UGetIndex(core, level, ci - cfactor , &iu, &is_stored);

               _braid_CoreFcn(core, clone)(app, ua[iu], &bigstep );
               time_left = ta[ci-cfactor-ilower];
            }
         }

          /* F-relaxation */
         for (fi = flo; fi <= fhi; fi++)
         {
            _braid_Step(core, level, fi, NULL, u);
            _braid_USetVector(core, level, fi, u, 0);
         }

         /* C-relaxation */
         if (ci > 0)
         {
            _braid_Step(core, level, ci, NULL, u);
            _braid_USetVector(core, level, ci, u, 1);

            if (  richardson )
             {
               dtk_index--;
               braid_Real DTK = pow( ta[ci-ilower] - time_left, order );
               braid_Real dtk_temp = dtk[dtk_index];
               a = DTK / ( DTK - dtk_temp );
               b = - dtk_temp / ( DTK - dtk_temp );

               _braid_StepStatusInit(time_left, ta[ci-ilower], ci-cfactor, tol, 
                                     iter, level, nrefine, gupper, status);            
               _braid_CoreFcn(core, step)(app, u, NULL, bigstep, status);              
               _braid_CoreFcn(core, sum)(app, a, u, b, bigstep );

                braid_Int iu, is_stored;
               _braid_UGetIndex(core, level, ci, &iu, &is_stored);

               _braid_CoreFcn(core, free)(app, ua[iu]);
               _braid_CoreFcn(core, clone)(app, bigstep, &ua[iu]);
               _braid_CoreFcn(core, free)(app,bigstep);
            }

         }

         /* if ((flo <= fhi) && (interval == ncpoints)) */
         if ((flo <= fhi) && !(ci > 0))
         {
            _braid_CoreFcn(core, free)(app, u);
         }
      }

      _braid_UCommWait(core, level);

      if (send_flag) 
      {
          MPI_Wait( &send_request, MPI_STATUS_IGNORE);
          _braid_TFree( send_buff );
      }

   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * F-Relax on level and restrict to level+1
 *
 * At the bottom of FRestrict(), the coarse-grid right-hand side is computed and
 * stored in c_fa[].  The result stored in c_fa[] must be considered carefully
 * to maintain a fixed point nature of the algorithm.  The FAS coarse-grid
 * equation is
 *
 *   A_c(u_c) = R(f - A(u)) + A_c(R(u))
 *
 * where R() is restriction, and A() and A_c() are the fine and coarse
 * operators.  We desire our algorithm to be a fixed-point method, in that the
 * exact solution (from sequential time stepping) should yield a zero initial
 * residual and that the residual should stay zero as Braid iterates.  This is
 * actually a subtle point.  In this case, the solution u_c should be equivalent
 * to the solution u.  Since R(f - A(u)) = 0, this implies that
 *
 *   A_c(u_c) = A_c(R(u))
 *
 * This further implies that the initial guess provided to the implicit solve
 * routines in A_c() must be consistent on the right an left hand sides.
 *
 * The function (@ref _braid_GetUInit()) provides this functionality by
 * encapsulating all decisions about how to set this initial guess (called
 * ustop) The ustop values for A_c(R(u)) are set in (@ _braid_Residual) and for
 * A_c(u_c), they are set in (@ _braid_step).
 *
 * storage == -2 implies that (@ref _braid_GetUInit()) will always use the
 * previous time step value as the initial guess, even if you have better
 * information, which is the case at C-points on the fine grid, and at all
 * points on coarse-levels.
 *
 * storage == -1 implies that (@ref _braid_GetUInit()) will always use the best
 * information.  The value at the previous time step is used as the initial
 * guess only at F-points on the fine-grid.  The va[] values provide all initial
 * guess at F- and C-points on coarse-grids.
 *
 * storage == 0 is the same as -1, except that storage also exists at fine-grid
 * F-points.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FRestrict(braid_Core   core,
                 braid_Int    level)
{
   MPI_Comm             comm        = _braid_CoreElt(core, comm);
   braid_App            app         = _braid_CoreElt(core, app);
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_AccessStatus   astatus     = (braid_AccessStatus)core;
   braid_Int            iter        = _braid_CoreElt(core, niter);
   braid_Int            print_level = _braid_CoreElt(core, print_level);
   braid_Int            access_level= _braid_CoreElt(core, access_level);
   braid_Int            tnorm       = _braid_CoreElt(core, tnorm);
   braid_Real          *tnorm_a     = _braid_CoreElt(core, tnorm_a);
   braid_Int            nrefine     = _braid_CoreElt(core, nrefine);
   braid_Int            gupper      = _braid_CoreElt(core, gupper);
   braid_Int            cfactor     = _braid_GridElt(grids[level], cfactor);
   braid_Int            ncpoints    = _braid_GridElt(grids[level], ncpoints);
   braid_Real          *ta          = _braid_GridElt(grids[level], ta);
   braid_Int            f_ilower    = _braid_GridElt(grids[level], ilower);
   
   braid_Int            richardson  = _braid_CoreElt(core, richardson);
   braid_Int            est_error   = _braid_CoreElt(core, est_error);
   braid_Int            order       = _braid_CoreElt(core, order);
   braid_Real          *ta_c        = _braid_GridElt(grids[1], ta );
   braid_Real          *dtk_core    = _braid_CoreElt(core, dtk);
   _braid_CommHandle   *recv_handle = NULL;
   _braid_CommHandle   *send_handle = NULL;

   braid_Int            c_level, c_ilower, c_iupper, c_index, c_i, c_ii;
   braid_Vector         c_u, *c_va, *c_fa;

   braid_Vector         u, r;
   braid_Int            interval, flo, fhi, fi, ci;
   braid_Real           rnorm, grnorm, rnorm_temp, rnm, factor, dtk, DTK;
   braid_Real          *estimate;    

   c_level  = level+1;
   c_ilower = _braid_GridElt(grids[c_level], ilower);
   c_iupper = _braid_GridElt(grids[c_level], iupper);
   c_va     = _braid_GridElt(grids[c_level], va);
   c_fa     = _braid_GridElt(grids[c_level], fa);

   rnorm = 0.0;

   /* Set the calling function */
   _braid_CoreElt(core, calling_function) = braid_Caller_FRestrict;
   
   _braid_UCommInit(core, level);

   /* Start from the right-most interval.
    * 
    * Do an F-relax and then a C-relax.  These relaxations are needed to compute
    * the residual, which is needed for the coarse-grid right-hand-side and for
    * convergence checking on the finest grid.  This loop updates va and fa. */
   for (interval = ncpoints; interval > -1; interval--)
   {
      _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);

      if (flo <= fhi)
      {
         _braid_UGetVector(core, level, flo-1, &r);
      }
      else if (ci > 0)
      {
         _braid_UGetVector(core, level, ci-1, &r);
      }

      /* F-relaxation */
      _braid_GetRNorm(core, -1, &rnm);
      for (fi = flo; fi <= fhi; fi++)
      {
         _braid_Step(core, level, fi, NULL, r);
         
         /* Allow user to process current vector, note that r here is
          * temporarily holding the state vector */
         if( (access_level >= 3) )
         {
            _braid_AccessStatusInit(ta[fi-f_ilower], fi, rnm, iter, level, nrefine, gupper,
                                    0, 0, astatus);
            _braid_AccessVector(core, astatus, r);
         }
         _braid_USetVector(core, level, fi, r, 0); /*BO: moved so that access modifes the vector */
      }

      /* Allow user to process current C-point */
      if( (access_level>= 3) && (ci > -1) )
      {
         _braid_UGetVectorRef(core, level, ci, &u);
         _braid_AccessStatusInit(ta[ci-f_ilower], ci, rnm, iter, level, nrefine, gupper,
                                 0, 0, astatus);
         _braid_AccessVector(core, astatus, u);
      }
      
      /* Compute residual and restrict */
      if (ci > 0)
      {
         /* Compute FAS residual */
         _braid_UGetVectorRef(core, level, ci, &u);
         _braid_FASResidual(core, level, ci, u, r);

         /* Compute rnorm (only on level 0). Richardson computes the rnorm later ) */
         if (level == 0 && !richardson )
         {
            _braid_CoreFcn(core, spatialnorm)(app, r, &rnorm_temp);
            tnorm_a[interval] = rnorm_temp;       /* inf-norm uses tnorm_a */
            if(tnorm == 1) 
            {  
               rnorm += rnorm_temp;               /* one-norm combination */ 
            }
            else if(tnorm == 2)
            {  
               rnorm += (rnorm_temp*rnorm_temp);  /* two-norm combination */
            }
         }

         /* Restrict u and residual, coarsening in space if needed */
         _braid_MapFineToCoarse(ci, cfactor, c_index);
         _braid_Coarsen(core, c_level, ci, c_index, u, &c_va[c_index-c_ilower]);
         _braid_Coarsen(core, c_level, ci, c_index, r, &c_fa[c_index-c_ilower]);
      }
      else if (ci == 0)
      {
         /* Restrict initial condition, coarsening in space if needed */
         _braid_UGetVectorRef(core, level, 0, &u);
         _braid_Coarsen(core, c_level, 0, 0, u, &c_va[0]);
      }

      if ((flo <= fhi) || (ci > 0))
      {
         _braid_CoreFcn(core, free)(app, r);
      }
   }
   _braid_UCommWait(core, level);

   
   /* Now apply coarse residual to update fa values */
   /* Set initial guess on coarse level */
   _braid_InitGuess(core, c_level);

   /* Initialize update of c_va[-1] boundary */
   if (c_ilower <= c_iupper)
   {
      _braid_CommRecvInit(core, c_level, c_ilower-1, &c_va[-1],
                          &recv_handle);
      _braid_CommSendInit(core, c_level, c_iupper, c_va[c_iupper-c_ilower],
                          &send_handle);
   }

   /* Start with rightmost point */

   if ( est_error )
   {
        estimate = _braid_CTAlloc(braid_Real, c_iupper-c_ilower + 1 );
   }

   for (c_i = c_iupper; c_i >= c_ilower; c_i--)
   {
      if (c_i > 0)
      {
         c_ii = c_i - c_ilower;
         if (c_ii == 0)
         {
            /* Finalize update of c_va[-1] */
            _braid_CommWait(core, &recv_handle);
         }
         _braid_CoreFcn(core, clone)(app, c_va[c_ii-1], &c_u);
         _braid_Residual(core, c_level, c_i, c_va[c_ii], c_u);
         
         if ( level == 0 && richardson  ) 
         {    
               dtk = dtk_core[c_ii];
               DTK = pow( ta_c[c_ii] - ta_c[c_ii-1], order );
               factor = DTK / ( DTK - dtk );
               _braid_CoreFcn(core, sum)(app, factor, c_u, factor, c_fa[c_ii]);   

               /* Compute the rnorm */
               _braid_CoreFcn(core, sum)(app, 1.0, c_fa[c_ii], -1.0, c_u );               
               _braid_CoreFcn(core, spatialnorm)(app, c_u, &rnorm_temp);
               
               tnorm_a[c_ii] = rnorm_temp;       /* inf-norm uses tnorm_a */
               if(tnorm == 1) 
               {  
                   rnorm += rnorm_temp;               /* one-norm combination */ 
               }
               else if(tnorm == 2)
               {  
                  rnorm += (rnorm_temp*rnorm_temp);  /* two-norm combination */
               }
         }    
         else
         {
            _braid_CoreFcn(core, sum)(app, 1.0, c_u, 1.0, c_fa[c_ii]);        
         }


         if ( level == 0 && est_error )
         {

              braid_Real est_temp;
              dtk = dtk_core[c_ii]  ;
              DTK = pow( ta_c[c_ii] - ta_c[c_ii-1], order );
             _braid_CoreFcn(core, spatialnorm)(app, c_fa[c_ii], &est_temp);
             if ( richardson )
             {
                est_temp = est_temp / DTK;
             }
             else
             {
                est_temp = est_temp / ( DTK - dtk );
             }
             estimate[ c_ii ] = est_temp; 
         }
              
         _braid_CoreFcn(core, free)(app, c_u);
      }
   }

   _braid_CommWait(core, &send_handle);

   /* Compute rnorm (only on level 0) */
   if (level == 0 )
   {
      if(tnorm == 1)          /* one-norm reduction */
      {  
         MPI_Allreduce(&rnorm, &grnorm, 1, braid_MPI_REAL, MPI_SUM, comm);
      }
      else if(tnorm == 3)     /* inf-norm reduction */
      {  
         _braid_Max(tnorm_a, ncpoints, &rnorm); 
         MPI_Allreduce(&rnorm, &grnorm, 1, braid_MPI_REAL, MPI_MAX, comm);
      }
      else                    /* default two-norm reduction */
      {  
         MPI_Allreduce(&rnorm, &grnorm, 1, braid_MPI_REAL, MPI_SUM, comm);
         grnorm = sqrt(grnorm);
      }

      /* Store new rnorm */
      _braid_SetRNorm(core, -1, grnorm); 
   }

   /* If debug printing, print out tnorm_a for this interval. This
    * should show the serial propagation of the exact solution */
   if ((print_level >= 2) && (level == 0) )
   {
      _braid_PrintSpatialNorms(core, tnorm_a, ncpoints);
   }

   /* Need to finalize the error estimates at the F-points */ 
   if ( level == 0 && est_error )
   {
      _braid_FinalizeErrorEstimates( core, estimate , c_iupper-c_ilower + 1 );
      _braid_TFree(estimate);
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * F-Relax on level and interpolate to level-1
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FInterp(braid_Core  core,
               braid_Int   level)
{
   braid_App          app          = _braid_CoreElt(core, app);
   _braid_Grid      **grids        = _braid_CoreElt(core, grids);
   braid_AccessStatus astatus      = (braid_AccessStatus)core;
   braid_Int          iter         = _braid_CoreElt(core, niter);
   braid_Int          access_level = _braid_CoreElt(core, access_level);
   braid_Int          nrefine      = _braid_CoreElt(core, nrefine);
   braid_Int          gupper       = _braid_CoreElt(core, gupper);
   braid_Int          ilower       = _braid_GridElt(grids[level], ilower);
   braid_Int          ncpoints     = _braid_GridElt(grids[level], ncpoints);
   braid_Vector      *va           = _braid_GridElt(grids[level], va);
   braid_Real        *ta           = _braid_GridElt(grids[level], ta);
   
   braid_Real         rnorm;
   braid_Int          f_level, f_cfactor, f_index;
   braid_Vector       f_u, f_e;

   braid_Vector       u, e;
   braid_Int          flo, fhi, fi, ci;
   braid_Int          interval;

   f_level   = level-1;
   f_cfactor = _braid_GridElt(grids[f_level], cfactor);

   /* Set the calling function */
   _braid_CoreElt(core, calling_function) = braid_Caller_FRestrict;
   
   _braid_GetRNorm(core, -1, &rnorm);
   
   _braid_UCommInitF(core, level);

   /**
    * Start from the right-most interval 
    *
    * First, generate the coarse-grid F-points through F-relaxation and
    * interpolate them to the fine grid, where they are C-points.  Second,
    * interpolate the coarse-grid C-points to the fine-grid.  The user-defined
    * spatial refinement (if set) is also called.  
    **/
   for (interval = ncpoints; interval > -1; interval--)
   {
      _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);

      /* Relax and interpolate F-points, refining in space if needed */
      if (flo <= fhi)
      {
         _braid_UGetVector(core, level, flo-1, &u);
      }
      for (fi = flo; fi <= fhi; fi++)
      {
         _braid_Step(core, level, fi, NULL, u);
         /* Allow user to process current vector */
         if( (access_level >= 3) )
         {
            _braid_AccessStatusInit(ta[fi-ilower], fi, rnorm, iter, level, nrefine, gupper,
                                    0, 0, astatus);
            _braid_AccessVector(core, astatus, u);
         }
         _braid_USetVector(core, level, fi, u, 0); /* BO: moved so access modifies */
         e = va[fi-ilower];
         _braid_CoreFcn(core, sum)(app, 1.0, u, -1.0, e);
         _braid_MapCoarseToFine(fi, f_cfactor, f_index);
         _braid_Refine(core, f_level, f_index, fi, e, &f_e);
         _braid_UGetVectorRef(core, f_level, f_index, &f_u);
         _braid_CoreFcn(core, sum)(app, 1.0, f_e, 1.0, f_u);
         _braid_USetVectorRef(core, f_level, f_index, f_u);
         _braid_CoreFcn(core, free)(app, f_e);
      }
      if (flo <= fhi)
      {
         _braid_CoreFcn(core, free)(app, u);
      }

      /* Interpolate C-points, refining in space if needed */
      if (ci > 0)
      {
         _braid_UGetVectorRef(core, level, ci, &u);
         /* Allow user to process current C-point */
         if( (access_level >= 3) )
         {
            _braid_AccessStatusInit(ta[ci-ilower], ci, rnorm, iter, level, nrefine, gupper,
                                    0, 0, astatus);
            _braid_AccessVector(core, astatus, u);
         }
         e = va[ci-ilower];
         _braid_CoreFcn(core, sum)(app, 1.0, u, -1.0, e);
         _braid_MapCoarseToFine(ci, f_cfactor, f_index);
         _braid_Refine(core, f_level, f_index, ci, e, &f_e);
         _braid_UGetVectorRef(core, f_level, f_index, &f_u);
         _braid_CoreFcn(core, sum)(app, 1.0, f_e, 1.0, f_u);
         _braid_USetVectorRef(core, f_level, f_index, f_u);
         _braid_CoreFcn(core, free)(app, f_e);
      }
   }

   _braid_UCommWait(core, level);

   /* Clean up */
   _braid_GridClean(core, grids[level]);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Refine in Space at every point on processor if r_space is set
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FRefineSpace(braid_Core   core,
                    braid_Int   *refined_ptr)
{

   MPI_Comm      comm    = _braid_CoreElt(core, comm);
   _braid_Grid **grids   = _braid_CoreElt(core, grids);   
   
   braid_Int     r_space = _braid_CoreElt(core, r_space);
   braid_Int     ilower  = _braid_GridElt(grids[0], ilower);
   braid_Int     iupper  = _braid_GridElt(grids[0], iupper);
   braid_Real   *ta      = _braid_GridElt(grids[0], ta); 

   braid_Int     global_r_space;
   braid_Int     i, ii;
   braid_Vector  c_vec, f_vec;  

   if ( _braid_CoreElt(core, scoarsen) != NULL )
   {
       
      if ( r_space > 0 )  
      {         
         for ( i = ilower; i <= iupper; i++ )
         {
            ii = i-ilower;
            _braid_UGetVectorRef(core, 0, i , &c_vec);
             
            if ( c_vec != NULL )
            {
               _braid_RefineBasic(core, -1, &ta[ii], &ta[ii], c_vec, &f_vec);
               _braid_USetVectorRef(core, 0, i, f_vec);
            }
         }                       
      }
    
      /* Check if any refinment was completed globally. If true then refine the 
         initial time point if not done already, increase nrefine, and return 2 */
      MPI_Allreduce(&r_space, &global_r_space, 1, braid_MPI_INT, MPI_MAX, comm);
      if (global_r_space > 0)
      { 
         /* Need to make sure to refine the first point. If it is a lone C point
            on a processor r_space then can never be set for that processor */ 
         if ( ilower == 0 && r_space == 0 )
         {
            _braid_UGetVectorRef(core, 0, 0, &c_vec);
            _braid_RefineBasic(core, -1, &ta[0], &ta[0], c_vec, &f_vec);
            _braid_USetVectorRef(core, 0, 0, f_vec);     
         }       
                
         *refined_ptr = 2;
         _braid_CoreElt(core, nrefine) += 1;
      }
      else
      {
         *refined_ptr = 0;
      }
   }
   else
   {
      *refined_ptr = 0;
   }
    
   /* Reset r_space */
   _braid_CoreElt(core, r_space) = 0;
   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Create a new fine grid (level 0) and corresponding grid hierarchy by refining
 * the current fine grid based on user-provided refinement factors.  Return the
 * boolean 'refined_ptr' to indicate whether grid refinement was actually done.
 * To simplify the algorithm, refinement factors are automatically capped to be
 * no greater than the coarsening factor (for level 0).  The grid data is also
 * redistributed to achieve good load balance in the temporal dimension.  If the
 * refinement factor is 1 in each time interval, no refinement is done.
 *
 * This routine is somewhat complex, but an attempt was made to use consistent
 * terminology throughout.  We refer to the initial level 0 grid as the "coarse"
 * grid and the new level 0 grid as the "fine" grid.  The routine starts with
 * the coarse grid and creates an intermediate "refined" grid that is co-located
 * with the parent coarse grid data distribution.  On this refined grid, a
 * number of things are constructed: the mapping between the coarse and fine
 * indexes; the new fine time values; and the injected and (possibly) spatially
 * refined coarse u-vectors.  This data is then redistributed by first sending
 * the mapping and time value information to the appropriate processors (the
 * receiving side polls for arbitrary messages until each fine index is
 * accounted for).  The u-vector values are communicated in a second phase
 * (without polling).  Finally, a new hierarchy is created and the fine grid
 * values are initialized by integrating the communicated u-vector values to the
 * next C-point to the right.  Note that in the case of C-point storage, some
 * u-vector values may not need to be communicated.
 *
 * The variable names use certain conventions as well.  No prefix (usually)
 * indicates a coarse variable, and the prefixes 'r_' and 'f_' indicate data on
 * the "refined" and "fine" grids, respectively.  Characters 'c', 'f', and 'a'
 * usually mean "coarse", "fine", and "array".  Indexes 'i', 'r_i', and 'f_i'
 * refer to global indexes, while 'ii', 'r_ii', and 'f_ii' are local indexes.
 *
 * Here are some important variables along with a brief description (they are
 * illustrated further in the example below):
 *
 *   gupper, f_gupper - global upper index for coarse and fine grids
 *   ilower,   iupper,   npoints   - extents and size of local coarse interval
 *   r_ilower, r_iupper, r_npoints - extents and size of local refined interval
 *   f_ilower, f_iupper, f_npoints - extents and size of local fine interval
 *
 *   r_ca - index map from fine to coarse on the refined grid  (size 'r_npoints')
 *   r_ta - time values on the refined grid                    (size 'r_npoints+2')
 *   r_fa - index map from coarse to fine on the refined grid  (size 'npoints+1')
 *          (note the extra value)
 *   f_ca - index map from fine to coarse on the fine grid     (size 'f_npoints')
 *
 *   send_ua - array of u-vectors to send to new processors    (size 'npoints')
 *   recv_ua - array of u-vectors received from old processors (size 'f_npoints')
 *
 *   f_first - fine grid index for the first coarse point in my fine interval
 *   f_next  - fine grid index for the first coarse point in the next fine interval
 *
 * Example: Some processor p owns the coarse interval, ilower = 29, iupper = 33.
 * The coarsening factor is 3 and 'rfactors' indicates the refinement factor for
 * the coarse interval to the left.  From this, an intermediate refined grid and
 * a final fine grid are formed as described above.
 *
 *   Coarse            |-----------c-----------|-----------|-----------c
 *     Grid           29          30          31          32          33
 * 
 * rfactors            3           1           2           3           2
 * 
 *  Refined ---c---|---|-----------c-----|-----|---c---|---|-----c-----|
 *     Grid   57  58  59          60    61    62  63  64  65    66    67
 * 
 * r_ilower   57
 *     r_ca   -1  -1  29          30    -1    31  -1  -1  32    -1    33
 *     r_ta    *   *   *           *     *     *   *   *   *     *     *
 *     r_fa           59          60          62          65          67          70
 *  send_ua            *           *           *           *           *
 * 
 *     Fine      |-----|---c---|---|-----c-----|---|---c---|-----|-----c
 *     Grid     61    62  63  64  65    66    67  68  69  70    71    72
 *              
 * f_ilower     61
 *     f_ca     -1    31  -1  -1  32    -1    33  -1  -1  34    -1    35
 *  f_first           62
 *   f_next                                                                       74
 *  recv_ua      0     *   0   0   *     0     *   0   0   *     0     *
 * 
 * When storing C-pts only, data from coarse indices 29 and 34 are not needed,
 * so we have the following differences from above:
 *     r_ca   -1  -1  -1          30    -1    31  -1  -1  32    -1    33
 *  send_ua            0           *           *           *           *
 *     f_ca     -1    31  -1  -1  32    -1    33  -1  -1  -1    -1    35
 *  recv_ua      0     *   0   0   *     0     *   0   0   0     0     *
 * 
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FRefine(braid_Core   core,
               braid_Int   *refined_ptr)
{
   MPI_Comm           comm            = _braid_CoreElt(core, comm);
   braid_App          app             = _braid_CoreElt(core, app);
   braid_Int          iter            = _braid_CoreElt(core, niter);
   braid_Int          nrefine         = _braid_CoreElt(core, nrefine);
   braid_AccessStatus astatus         = (braid_AccessStatus)core;
   braid_BufferStatus bstatus         = (braid_BufferStatus)core;
   braid_Int          access_level    = _braid_CoreElt(core, access_level);
   _braid_Grid      **grids           = _braid_CoreElt(core, grids);
   braid_Int          ncpoints        = _braid_GridElt(grids[0], ncpoints);
   braid_Int          myid = _braid_CoreElt(core, myid );
   braid_Real     rnorm;

   /* Use prefix 'r_' for the refined version of the current grid level 0, and
    * use prefix 'f_' for the fine grid in the new distribution. */

   braid_Int      npoints, ilower, iupper, gupper, ii;
   braid_Int      r_ilower, r_i, r_ii;
   braid_Int      f_npoints, f_ilower, f_iupper, f_gupper, f_i, f_j, f_ii;
   braid_Int     *r_ca, *r_fa, *f_ca, f_first, f_next, next;
   braid_Real    *ta, *r_ta, *f_ta;

   braid_Vector  *send_ua, *recv_ua, u;
   braid_Int     *send_procs, *send_unums;
   braid_Int     *recv_procs, *recv_unums, *recv_map;
   braid_Real    **send_buffers, **recv_buffers, *bptr;
   void          *buffer;
   braid_Int      send_size, recv_size, size, max_usize;
   braid_Int     nsends, nrecvs, proc, prevproc;
   braid_Int      unum, send_msg, recv_msg;
   MPI_Request   *requests;
   MPI_Status    *statuses;

   _braid_Grid   *f_grid;
   braid_Int      done, cfactor;
   braid_Int     *send_map, *left_procs, *right_procs;
   braid_Int      m, interval, flo, fhi, fi, ci, f_hi, f_ci;
   
   /* Set the calling function */
   _braid_CoreElt(core, calling_function) = braid_Caller_FRefine;
   _braid_GetCFactor(core,0, &cfactor);
#if DEBUG
   braid_Int  myproc;
   /*cfactor = 6;*/ /* RDF HACKED TEST */
   MPI_Comm_rank(comm, &myproc);
#endif

   /*-----------------------------------------------------------------------*/
   /* 1. Compute f_gupper and the local interval extents for both the refined
    * and fine grids.  The local refined interval contains the fine grid points
    * underlying the coarse interval (ilower-1, iupper]. */

   /* Get the new refined, balanced distrobution */
   done = 0;
   
   _braid_BalanceStruct *bstruct;
   bstruct = _braid_CTAlloc( _braid_BalanceStruct , 1 );
   _braid_BuildBalanceStructure( core, &done, bstruct);
   
   /* If done then no load balancing or time refinement is necessary */
   if ( done )
   {
      _braid_DestroyBalanceStruct( bstruct );
      _braid_FRefineSpace(core, refined_ptr);
      return _braid_error_flag;
   }
       
    *refined_ptr = _braid_BalanceElt(bstruct, refine);
   _braid_CoreElt(core, r_space) = 0 ;
      
   /*current grid */
   ta       = _braid_GridElt( grids[0], ta);
   gupper   = _braid_GridElt( grids[0], gupper);
   ilower   = _braid_GridElt( grids[0], ilower);
   iupper   = _braid_GridElt( grids[0], iupper);
   npoints  = iupper - ilower + 1;
    
   /* refined grid */
   r_ilower  = _braid_BalanceElt( bstruct, refined_ilower );
   r_ta      = _braid_BalanceElt( bstruct, refined_ta);

   /* New fine grid */
   f_ilower    = _braid_BalanceElt( bstruct, fine_ilower );
   f_iupper    = _braid_BalanceElt( bstruct, fine_iupper );
   f_gupper    = _braid_BalanceElt( bstruct, fine_gupper );
   f_ta        = _braid_BalanceElt( bstruct, fine_ta);
   f_npoints   = f_iupper - f_ilower + 1;

   /* Grid maps (see introduction) */
   r_ca      = _braid_BalanceElt( bstruct, refined_ca);
   r_fa      = _braid_BalanceElt( bstruct, refined_fa);
   f_ca      = _braid_BalanceElt( bstruct, fine_ca);
   f_next    = _braid_BalanceElt( bstruct, fine_next);
   f_first   = _braid_BalanceElt( bstruct, fine_first);

   /* Communication maps */
   recv_map    = _braid_BalanceElt( bstruct, recv_map );    /* Proc each time step comes from */
   send_map    = _braid_BalanceElt( bstruct, send_map );    /* Proc each time step goes to */
   left_procs  = _braid_BalanceElt( bstruct, left_procs );  /* Proc to my left on each level */
   right_procs = _braid_BalanceElt( bstruct, right_procs ); /* Proc to my right on each level */
   
#if DEBUG   
   braid_Int r_iupper  = _braid_BalanceElt( bstruct, refined_iupper );
   printf(" %d : ilower %d iupper %d r_ilower %d r_iupper %d f_ilower %d f_iupper %d \n ",
          myid,     ilower,   iupper,   r_ilower,   r_iupper,   f_ilower,   f_iupper); 
#endif

   /* 1. Initialize the new fine grid and copy over the time values */
   _braid_GridInit(core, 0, f_ilower, f_iupper, &f_grid);
   braid_Real *f_ta_grid = _braid_GridElt(f_grid, ta);
   memcpy(f_ta_grid, f_ta, sizeof(braid_Real)*f_npoints);
 
   /*-----------------------------------------------------------------------*/
   /* 2. Build u-vectors on the fine grid (send_ua) by first integrating on the
    * coarse grid, then injecting and refining spatially.  Redistribute these
    * u-vectors to the fine grid (recv_ua). */

   send_ua = _braid_CTAlloc(braid_Vector, npoints);
   send_procs = _braid_CTAlloc(braid_Int, npoints);
   send_unums = _braid_CTAlloc(braid_Int, npoints);
   send_buffers = _braid_CTAlloc(braid_Real *, npoints);

   recv_ua = _braid_CTAlloc(braid_Vector, f_npoints);
   recv_procs = _braid_CTAlloc(braid_Int, f_npoints);
   recv_unums = _braid_CTAlloc(braid_Int, f_npoints);
   recv_buffers = _braid_CTAlloc(braid_Real *, f_npoints);

   _braid_GetRNorm(core, -1, &rnorm);

   _braid_UCommInitF(core, 0);

   /* Start from the right-most interval */
   for (interval = ncpoints; interval > -1; interval--)
   {
      _braid_GetInterval(core, 0, interval, &flo, &fhi, &ci);

      /* Integrate F-points and refine in space */
      if (flo <= fhi)
      {
         _braid_UGetVector(core, 0, flo-1, &u);
         for (fi = flo; fi <= fhi; fi++)
         {
            _braid_Step(core, 0, fi, NULL, u);
            _braid_USetVector(core, 0, fi, u, 0); /* needed for communication */

            /* Set send_ua */
            ii = fi - ilower;
            r_ii = r_fa[ii] - r_ilower;
            if (r_ca[r_ii] > -1)
            {
               _braid_RefineBasic(core, -1, &r_ta[r_ii], &ta[ii], u, &send_ua[ii]);
            }

            /* Allow user to process current vector */
            if( (access_level >= 3) )
            {
               _braid_AccessStatusInit(ta[ii], fi, rnorm, iter, 0, nrefine, gupper,
                                       0, 0, astatus);
               _braid_AccessVector(core, astatus, u);
            }
         }
         _braid_CoreFcn(core, free)(app, u);
      }

      /* Refine C-points in space */
      if (ci > -1)
      {
         _braid_UGetVectorRef(core, 0, ci, &u);

         /* Set send_ua */
         ii = ci - ilower;
         r_ii = r_fa[ii] - r_ilower;
         if (r_ca[r_ii] > -1)
         {
            _braid_RefineBasic(core, -1, &r_ta[r_ii], &ta[ii], u, &send_ua[ii]);
         }

         /* Allow user to process current vector */
         if( (access_level >= 3) )
         {
            _braid_AccessStatusInit(ta[ii], ci, rnorm, iter, 0, nrefine, gupper,
                                    0, 0, astatus);
            _braid_AccessVector(core, astatus, u);
         }
      }
   }

   _braid_UCommWait(core, 0);

   /* Compute nsends, send_procs, and send_unums from send_ua array */
   nsends = -1;
   prevproc = -1;
   for (ii = 0; ii < npoints; ii++)
   {
      if (send_ua[ii] != NULL)
      {
         r_i = r_fa[ii];
         proc = send_map[r_i - r_ilower];
         if (proc != prevproc)
         {
            nsends++;
            send_procs[nsends] = proc;
            send_unums[nsends] = 0;
            prevproc = proc;
         }
         send_unums[nsends]++;
      }
   }
   nsends++;

   /* Compute nrecvs, recv_procs, and recv_unums from f_ca array */
   nrecvs = -1;
   prevproc = -1;
   for (f_ii = 0; f_ii < f_npoints; f_ii++)
   {
      if (f_ca[f_ii] > -1)
      {
         proc = recv_map[f_ii];
         if (proc != prevproc)
         {
            nrecvs++;
            recv_procs[nrecvs] = proc;
            recv_unums[nrecvs] = 0;
            prevproc = proc;
         }
         recv_unums[nrecvs]++;
      }
   }
   nrecvs++;

   requests = _braid_CTAlloc(MPI_Request, (nsends+nrecvs));
   statuses = _braid_CTAlloc(MPI_Status,  (nsends+nrecvs));

   _braid_BufferStatusInit( 1, 0, bstatus );
   _braid_CoreFcn(core, bufsize)(app, &max_usize, bstatus); /* max buffer size */
   _braid_NBytesToNReals(max_usize, max_usize);

   
   /* Post u-vector receives */
   for (m = 0; m < nrecvs; m++)
   {
         if (recv_procs[m] != myid )
         {
            unum = recv_unums[m]; /* Number of u-vectors being received */
            recv_size = unum*(1 + max_usize);
            recv_buffers[m] = _braid_CTAlloc(braid_Real, recv_size);
            MPI_Irecv(recv_buffers[m], recv_size, braid_MPI_REAL, recv_procs[m], 5, comm,
                      &requests[m]);
         }
         else
         {
            /* Dont need to post a recv from myself! */
            requests[m] = MPI_REQUEST_NULL;
            recv_buffers[m] = NULL;
         }
   }
   
   /* Post u-vector sends */
   ii = 0;
   for (m = 0; m < nsends; m++)
   {
      unum = send_unums[m]; /* Number of u-vectors being sent */
      
      if ( send_procs[m] != myid )
      {   
         send_size = unum*(1 + max_usize);
         send_buffers[m] = _braid_CTAlloc(braid_Real, send_size);
         send_size = 0; /* Recompute send_size and realloc buffer */
         bptr = send_buffers[m];
         while (unum > 0)
         {
            if (send_ua[ii] != NULL)
            {
               /* Pack u into buffer, adjust size, and put size into buffer */
               buffer = &bptr[1];
               _braid_CoreFcn(core, bufsize)(app, &size, bstatus);
               _braid_StatusElt( bstatus, size_buffer ) = size;
               _braid_CoreFcn(core, bufpack)(app, send_ua[ii], buffer, bstatus);
               size = _braid_StatusElt(bstatus, size_buffer);
               _braid_CoreFcn(core, free)(app, send_ua[ii]);
               send_ua[ii] = NULL;
               _braid_NBytesToNReals(size, size);
               bptr[0] = (braid_Int) size; /* insert size at the beginning */
               bptr += (1+size);
               send_size += (1+size);
               unum--;
            }
            ii++;
         }
      
         send_buffers[m] = _braid_TReAlloc(send_buffers[m], braid_Real, send_size);
         MPI_Isend(send_buffers[m], send_size, braid_MPI_REAL, send_procs[m], 5, comm,
                   &requests[ nrecvs + m]);
      }
      else
      {
         /* We do not need to send a message to ourselves. */
         requests[nrecvs+m] = MPI_REQUEST_NULL;
         send_buffers[m] = NULL;
         while ( unum > 0 )
         {
            if (send_ua[ii] != NULL)
            {
               unum--;
            }
            ii++;
         }  
      }

   }
   /* Finish communication */
   MPI_Waitall(nrecvs+nsends, requests, statuses);
   
   /* Unpack u-vectors */
   f_ii = 0;
   ii = 0;
   for (m = 0; m < nrecvs; m++)
   {
      unum = recv_unums[m];
      bptr = recv_buffers[m];
      if ( recv_procs[m] != myid ) 
      {
         while (unum > 0)
         {
            if (f_ca[f_ii] > -1)
            {
               /* Unpack buffer into u-vector */
               buffer = &bptr[1];
               _braid_CoreFcn(core, bufunpack)(app, buffer, &recv_ua[f_ii], bstatus);
               size = (braid_Int) bptr[0];
               bptr += (1+size);
               unum--;
            }
            f_ii++;
         }
      }
      else
      {
         /* These time steps already exist on this processor. Time steps were
          * deleted from send_ua after packing, so send_ua only contains time 
          * steps that are to remain on this processor. */
         while ( unum > 0 )
         {
            if ( f_ca[f_ii] > -1 )
            {
               while ( send_ua[ii] == NULL )
               {
                  ii++;
               }
               recv_ua[f_ii] = send_ua[ii];
               ii++;
               unum--;
            }
            f_ii++;
         }
      }
   }
   
   /* Free up some memory */
   _braid_TFree(send_procs);
   _braid_TFree(send_unums);
   for (m = 0; m < nsends; m++)
   {
      _braid_TFree(send_buffers[m]);
   }
   _braid_TFree(send_buffers);
   _braid_TFree(recv_procs);
   _braid_TFree(recv_unums);
   for (m = 0; m < nrecvs; m++)
   {
      _braid_TFree(recv_buffers[m]);
   }
   _braid_TFree(recv_buffers);
   _braid_TFree(requests);
   _braid_TFree(statuses);
   
   braid_Int  level, nlevels = _braid_CoreElt(core, nlevels);
   _braid_TFree(_braid_CoreElt(core, rfactors));
   _braid_TFree(_braid_CoreElt(core, wfactors));
      
   if ( _braid_CoreElt(core, estimate ) != NULL )
   {
      _braid_TFree(_braid_CoreElt(core, estimate));
   }
   if ( _braid_CoreElt(core, dtk) != NULL )
   {
      _braid_TFree(_braid_CoreElt(core, dtk));
   }
   _braid_TFree(_braid_CoreElt(core, tnorm_a));

   for (level = 0; level < nlevels; level++)
   {
      _braid_GridDestroy(core, grids[level]);
   }

   /*-----------------------------------------------------------------------*/
   /* 5. Build the new fine grid hierarchy, then use recv_ua to populate the
    * initial values on grid level 0.  This is done by integrating values to the
    * next C-point the right.  Because we require that rfactor <= cfactor, each
    * C-point either has a corresponding value or has a value in the F-interval
    * immediately to the left that can be integrated.  Communication from the
    * left processor may still be needed. */

#if DEBUG
   printf("%d %d: 5\n", FRefine_count, myproc);
#endif

   /* Initialize new hierarchy */
   _braid_CoreElt(core, gupper)  = f_gupper;

   /*braid_SetCFactor(core,  0, cfactor);*/ /* RDF HACKED TEST */
   _braid_InitHierarchy(core, f_grid, 1, left_procs, right_procs);

   /* Initialize communication */
   recv_msg = 0;
   send_msg = 0;
   if (f_first > _braid_NextCPoint(f_ilower, cfactor))
   {
      recv_msg = 1;
   }
   if (f_next > _braid_NextCPoint(f_iupper+1, cfactor))
   {
      send_msg = 1;
   }
   _braid_UCommInitBasic(core, 0, recv_msg, send_msg, 0);

#if DEBUG
   printf("%d %d: 6 recv_msg = %d, send_msg = %d\n",
          FRefine_count, myproc, recv_msg, send_msg);
#endif

   /* Start from the right-most point */
   f_i = f_iupper;
   next = f_next;
   while (f_i >= f_ilower)
   {
      /* Find the next value to the left */
      u = NULL;
      for ( ; f_i >= f_ilower; f_i--)
      {
         f_ii = f_i - f_ilower;
         if (recv_ua[f_ii] != NULL)
         {
            u = recv_ua[f_ii];
            break;
         }
      }
      if ((f_i < f_ilower) && (recv_msg))
      {
         f_i = f_ilower-1; /* receive value from left processor */
         _braid_UGetVector(core, 0, f_i, &u);
      }

      /* Integrate the value if needed and set */
      if (u != NULL)
      {
         f_j = f_i;
         f_ci = _braid_PriorCPoint(next-1, cfactor);
         if (f_i < f_ci)
         {
            /* integrate */
            f_hi = _braid_min(f_ci, f_iupper);
            for ( ; f_j < f_hi; f_j++)
            {
               _braid_USetVector(core, 0, f_j, u, 0);
               _braid_Step(core, 0, f_j+1, NULL, u);
            }
         }
         _braid_USetVector(core, 0, f_j, u, 1);
         next = f_i;
         f_i--;
      }
   }


    _braid_CoreElt(core, nrefine) += *refined_ptr;
#if DEBUG
   printf("%d %d: 7\n", FRefine_count, myproc);
   fflush(stdout);
#endif

   _braid_UCommWait(core, 0);

#if DEBUG
   printf("%d %d: 8\n", FRefine_count, myproc);
   fflush(stdout);
   FRefine_count++;
#endif

   /* Free up remaining memory */
   _braid_TFree(send_ua);
   _braid_TFree(recv_ua);
   _braid_DestroyBalanceStruct( bstruct );

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Access to XBraid on grid level
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FAccess(braid_Core     core,
               braid_Int      level,
               braid_Int      done)
{
   braid_App           app         = _braid_CoreElt(core, app);
   _braid_Grid       **grids       = _braid_CoreElt(core, grids);
   braid_AccessStatus  astatus     = (braid_AccessStatus)core;
   braid_Int           iter        = _braid_CoreElt(core, niter);
   braid_Int           nrefine     = _braid_CoreElt(core, nrefine);
   braid_Int           gupper      = _braid_CoreElt(core, gupper);
   braid_Int           access_level= _braid_CoreElt(core, access_level);
   braid_Int           ncpoints    = _braid_GridElt(grids[level], ncpoints);
   braid_Real          *ta         = _braid_GridElt(grids[level], ta);
   braid_Int           ilower      = _braid_GridElt(grids[level], ilower);

   braid_Real     rnorm;
   braid_Vector   u;
   braid_Int      interval, flo, fhi, fi, ci;

   /* Set the calling function */
   _braid_CoreElt(core, calling_function) = braid_Caller_FAccess;
   
   _braid_UCommInitF(core, level);
   
   _braid_GetRNorm(core, -1, &rnorm);

   /* Start from the right-most interval */
   for (interval = ncpoints; interval > -1; interval--)
   {
      _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);

      /* Give access at F-points */
      if (flo <= fhi)
      {
         _braid_UGetVector(core, level, flo-1, &u);
      }
      for (fi = flo; fi <= fhi; fi++)
      {
         _braid_Step(core, level, fi, NULL, u);
         _braid_USetVector(core, level, fi, u, 0);

         if (access_level >= 1)
         {
            _braid_AccessStatusInit( ta[fi-ilower], fi, rnorm, iter, level, nrefine, gupper,
                                     done, 0, astatus);
            _braid_AccessVector(core, astatus, u);
         }
      }
      if (flo <= fhi)
      {
         _braid_CoreFcn(core, free)(app, u);
      }

      /* Give access at C-points */
      if ((ci > -1) && (access_level >= 1))
      {
         _braid_UGetVectorRef(core, level, ci, &u);
         _braid_AccessStatusInit( ta[ci-ilower], ci, rnorm, iter, level, nrefine, gupper,
                                  done, 0, astatus);
         _braid_AccessVector(core, astatus, u);
      }
   }
   _braid_UCommWait(core, level);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Initialize grid hierarchy
 *----------------------------------------------------------------------------*/

braid_Int
_braid_InitHierarchy(braid_Core    core,
                     _braid_Grid  *fine_grid,
                     braid_Int     refined,
                     braid_Int    *recv_procs,
                     braid_Int    *send_procs)
{
   MPI_Comm       comm       = _braid_CoreElt(core, comm);
   braid_Int      max_levels = _braid_CoreElt(core, max_levels);
   braid_Int      min_coarse = _braid_CoreElt(core, min_coarse);
   braid_Int     *nrels      = _braid_CoreElt(core, nrels);
   braid_Int      nrdefault  = _braid_CoreElt(core, nrdefault);
   braid_Int      gupper     = _braid_CoreElt(core, gupper);
   braid_Int     *rfactors   = _braid_CoreElt(core, rfactors);
   braid_Real    *wfactors   = _braid_CoreElt(core, wfactors);
   
   braid_Int      nlevels    = _braid_CoreElt(core, nlevels);
   _braid_Grid  **grids      = _braid_CoreElt(core, grids);

   braid_Int      richardson = _braid_CoreElt(core, richardson);
   braid_Int      est_error  = _braid_CoreElt(core, est_error); 
   braid_Real    *dtk, *estimate;

   /**
    * These are some common index names used to refer to intervals and
    * time points.  Here's what they mean.
    *
    * cfactor            - coarsening factor, fixed on each level
    * ilower, iupper     - lowest and highest time indices on a level for one processor, 
    *                      could be C or F points
    * clower, cupper     - lowest and highest C-point indices on a level for one processor
    *                      analagous to ilower, iupper being projected onto C-points
    * clo, chi           - the coarse level indices for clower and cupper
    * f_iupper, f_ilower - lowest and highest F-point indices on a level
    * c_iupper, c_ilower - lowest and highest time indices on the coarse level 
    *                      (analagous to ilower and iupper, but on the next level down)
    * flo, fhi, ci       - describes an interval [ci, flo, flo+1, ..., fhi] 
    * gclower, gcupper   - global lowest and highest C-point indices on a level
    **/

   braid_Int      level;
   braid_Int      ilower, iupper;
   braid_Int      clower, cupper, cfactor, ncpoints, nupoints;
   braid_Real    *ta;
   braid_Vector  *ua;
   braid_Vector  *va;
   braid_Vector  *fa;

   _braid_Grid   *grid;
   braid_Real    *f_ta;
   braid_Int      i, f_i, f_ilower, clo, chi, gclower, gcupper;

   MPI_Request    request1, request2;
   MPI_Status     status;
   braid_Int      left_proc, right_proc;

   grids[0] = fine_grid;

   /* Do sequential time marching if min_coarse is already reached */
   if ( gupper <= min_coarse )
   {
      max_levels = 1;
   }

   /* Allocate space for rfactors (and initialize to zero) */
   ilower = _braid_GridElt(grids[0], ilower);
   iupper = _braid_GridElt(grids[0], iupper);
   rfactors = _braid_CTAlloc(braid_Int, iupper-ilower+2); /* Ensures non-NULL */
   wfactors = _braid_CTAlloc(braid_Real, iupper-ilower+2); /* Ensures non-NULL */
   for(i = 0; i < iupper-ilower+2; i++)
   {
      /* We need to ensure an rfactor of 1 for global index 0, and for
       * the case of a user-defined residual function and F-relaxation,
       * a default rfactor value at F-points */
      rfactors[i] = 1;
      wfactors[i] = 1;
      
   }
   _braid_CoreElt(core, rfactors) = rfactors;
   _braid_CoreElt(core, wfactors) = wfactors;
   

   /* Set up nrels array */
   for (level = 0; level < max_levels; level++)
   {
      if (nrels[level] < 0)
      {
         nrels[level] = nrdefault;
      }
   }

   /* Coarsen global grid to determine nlevels */
   gclower = 0;
   gcupper = gupper;
   for (level = 0; level < max_levels; level++)
   {
      grid = grids[level];
      ilower = _braid_GridElt(grid, ilower);
      iupper = _braid_GridElt(grid, iupper);
      if (level > 0)
      {
         /* Copy ta info from level-1 grid */
         ta       = _braid_GridElt(grid, ta);
         f_ilower = _braid_GridElt(grids[level-1], ilower);
         f_ta     = _braid_GridElt(grids[level-1], ta);
         cfactor  = _braid_GridElt(grids[level-1], cfactor);
         for (i = ilower; i <= iupper; i++)
         {
            _braid_MapCoarseToFine(i, cfactor, f_i);
            ta[i-ilower] = f_ta[f_i-f_ilower];
         }
      }

      _braid_GetCFactor(core, level, &cfactor);
      
      _braid_GridElt(grid, gupper) = gcupper;

      _braid_ProjectInterval(gclower, gcupper, 0, cfactor, &gclower, &gcupper);
      _braid_MapFineToCoarse(gclower, cfactor, gclower);
      _braid_MapFineToCoarse(gcupper, cfactor, gcupper);

      /* Coarsen */
      _braid_ProjectInterval(ilower, iupper, 0, cfactor, &clower, &cupper);
      _braid_MapFineToCoarse(clower, cfactor, clo);
      _braid_MapFineToCoarse(cupper, cfactor, chi);
      ncpoints = chi-clo+1;
      if (ncpoints < 0)
      {
         ncpoints = 0;
      }
      _braid_GridElt(grid, clower)   = clower;
      _braid_GridElt(grid, cupper)   = cupper;
      _braid_GridElt(grid, cfactor)  = cfactor;
      _braid_GridElt(grid, ncpoints) = ncpoints;
      if ( (gclower < gcupper) && (max_levels > level+1) &&
           ((gcupper - gclower) >= min_coarse) )
      {
         /* Initialize the coarse grid */
         _braid_GridInit(core, level+1, clo, chi, &grids[level+1]);
      }
      else
      {
         /* This is the coarsest level */
         if ( (level > 0) || (!refined) )
         {
            /* If this is a true coarse level (it has a fine grid above it in
             * the current hierarchy) or it is a fine level that was not built
             * by refining a coarser grid, then do serial time integration by
             * setting only one C-point and the rest F-points */
            if (ilower == 0)
            {
               ncpoints = 1;
            }
            else
            {
               ncpoints = 0;
            }
            /* clower > cupper indicates empty interval */
            _braid_GridElt(grid, clower)   = ilower;
            _braid_GridElt(grid, cupper)   = 0;
            _braid_GridElt(grid, cfactor)  = gupper+1;
            _braid_GridElt(grid, ncpoints) = ncpoints;
         }

         /* Stop coarsening */
         break;
      }
      
      if(level == 0)
      {   
         /* Allocate space for storage of residual norm at each C-point */
         _braid_CoreElt(core, tnorm_a)  = _braid_CTAlloc(braid_Real, ncpoints);
      }
   }
   nlevels = level+1;
   _braid_CoreElt(core, nlevels) = nlevels;

   /* Allocate ua, va, and fa here */
   for (level = 0; level < nlevels; level++)
   {
      grid = grids[level];
      ilower = _braid_GridElt(grid, ilower);
      iupper = _braid_GridElt(grid, iupper);
      if (level > 0)
      {
         va = _braid_CTAlloc(braid_Vector, iupper-ilower+2);
         fa = _braid_CTAlloc(braid_Vector, iupper-ilower+2);
         _braid_GridElt(grid, va_alloc) = va;
         _braid_GridElt(grid, fa_alloc) = fa;
         _braid_GridElt(grid, va)       = va+1;  /* shift */
         _braid_GridElt(grid, fa)       = fa+1;  /* shift */
      }

      // If on level that only stores C-points and not using the shell vector feature
      if ( ((_braid_CoreElt(core, storage) < 0) ||
            (level < _braid_CoreElt(core, storage))) &&
           (_braid_CoreElt(core, useshell)!=1) )
      {
         nupoints = _braid_GridElt(grid, ncpoints);   /* only C-points */
      }
      else
      {
         nupoints = iupper-ilower+1;                  /* all points */
      }
      ua = _braid_CTAlloc(braid_Vector, nupoints+1);
      _braid_GridElt(grid, nupoints)  = nupoints;
      _braid_GridElt(grid, ua_alloc)  = ua;
      _braid_GridElt(grid, ua)        = ua+1;  /* shift */
   }
   
   /* Set the neighbouring processor on each level */
   if (!(recv_procs && send_procs) )
   {
      _braid_SetInitNeighbours(core, nlevels, gupper, grids);
   }
   else
   {
      for (i = 0 ; i < nlevels; i++ )
      {
         _braid_GridElt( grids[i], left_proc ) = recv_procs[i] ;
         _braid_GridElt( grids[i], right_proc ) = send_procs[i] ;
      }
   }

   /* Communicate ta[-1] and ta[iupper-ilower+1] information */
   for (level = 0; level < nlevels; level++)
   {
      grid = grids[level];
      ilower = _braid_GridElt(grid, ilower);
      iupper = _braid_GridElt(grid, iupper);
      ta     = _braid_GridElt(grid, ta);

      if (ilower <= iupper)
      {
         _braid_GetProcLeftOrRight(core, level, -1, &left_proc);
         _braid_GetProcLeftOrRight(core, level, 1, &right_proc);

         /* Post receive to set ta[-1] on each processor*/
         if (left_proc > -1)
         {
            MPI_Irecv(&ta[-1], sizeof(braid_Real), MPI_BYTE,
                      left_proc, 1, comm, &request1);
         }
         else
         {
            /* Place a repeat value to indicate the start of the time-line for this level */
            ta[-1] = ta[0]; 
         }
         /* Post receive to set ta[iupper-ilower+1] on each processor */
         if ( _braid_CoreElt(core, scoarsen) != NULL )
         {
            if (right_proc > -1)
            {
               MPI_Irecv(&ta[iupper-ilower+1], sizeof(braid_Real), MPI_BYTE,
                         right_proc, 1, comm, &request2);
            }
            else
            {
               /* Place a repeat value to indicate the end the time-line for this level */
               ta[iupper-ilower+1] = ta[iupper-ilower];
            }
         }

         /* Post send that sets ta[-1] on each processor */
         if (right_proc > -1)
         {
            MPI_Send(&ta[iupper-ilower], sizeof(braid_Real), MPI_BYTE,
                     right_proc, 1, comm);
         }
         /* Post send that sets ta[iupper-ilower+1] on each processor */
         if ( (left_proc > -1) && ( _braid_CoreElt(core, scoarsen) != NULL ) )
         {
            MPI_Send(&ta[0], sizeof(braid_Real), MPI_BYTE, left_proc, 1, comm);
         }

         /* Finish receive */
         if (left_proc > -1)
         {
            MPI_Wait(&request1, &status);
         }
         if ( (right_proc > -1) && ( _braid_CoreElt(core, scoarsen) != NULL ) )
         {
            MPI_Wait(&request2, &status);
         }
      }
   }
   
   /* Allocate the dtk and estimate arrays if we are doing either richardson
    * or error estimation */
   if ( richardson || est_error )
   {
      grid = grids[0];
      ilower   = _braid_GridElt(grid, ilower);
      iupper   = _braid_GridElt(grid, iupper);
      ncpoints = _braid_GridElt(grid, ncpoints);

      dtk = _braid_CTAlloc(braid_Real, ncpoints);
      _braid_CoreElt(core, dtk) = dtk;
      if (est_error)
      {
         estimate = _braid_CTAlloc(braid_Real, iupper-ilower+1);
         _braid_CoreElt(core, estimate) = estimate;
      }   
      
      _braid_GetDtk(core);
   }
   return _braid_error_flag;
}


braid_Int
_braid_FinalizeErrorEstimates( braid_Core   core, 
                               braid_Real  *estimate,
                               braid_Int    length)
{

   braid_Int     order       = _braid_CoreElt(core, order );
   braid_Int     myid        = _braid_CoreElt(core, myid );
   MPI_Comm      comm        = _braid_CoreElt(core, comm );
  _braid_Grid  **grids       = _braid_CoreElt(core, grids);
   braid_Real   *error_est   = _braid_CoreElt(core, estimate);
   braid_Real   *ta       = _braid_GridElt(grids[0], ta);
   braid_Int     ilower   = _braid_GridElt(grids[0],ilower);
   braid_Int     iupper   = _braid_GridElt(grids[0],iupper);
   braid_Int     gupper   = _braid_GridElt(grids[0],gupper);
   braid_Int     cfactor  = _braid_GridElt(grids[0],cfactor);
   braid_Int     ncpoints = _braid_GridElt(grids[0], ncpoints);

   braid_Int recv_flag, send_flag_l, send_flag_r, i, last_cpoint, estimate_index ; 
   braid_Real recv_value, send_value_l, send_value_r, dtk, factor ;
   MPI_Request recv_request, send_request_l, send_request_r;

   if ( ilower <= iupper )
   {

      /* Recv situation 1: Recv from left if ilower is > last C point */
      send_flag_l = send_flag_r = recv_flag = 0; 
      last_cpoint = ( gupper / cfactor ) * cfactor ;
      if ( ilower > last_cpoint && ilower > 0 )
      { 
         recv_flag++;
         MPI_Irecv( &recv_value, 1, braid_MPI_REAL, myid -1, 98, comm, &recv_request );
      }
      
      /* Recv situation 2: recv from right if iupper is a f point */
      else if ( !_braid_IsCPoint(iupper, cfactor) && iupper != gupper )
      {
         recv_flag++;
         MPI_Irecv( &recv_value, 1, braid_MPI_REAL, myid + 1, 98, comm, &recv_request );
      }
      
      /* Send situation 1: Send to the right if proc owns a point in last interval, but not gupper */
      if ( iupper >= last_cpoint && iupper != gupper )
      {
         /* If we have no c points, wait to recv before sending */ 
         if ( ncpoints == 0 )
         {  
            recv_flag--;
            MPI_Wait( &recv_request , MPI_STATUS_IGNORE );
            send_value_r = recv_value;
         }
         else
         {
            send_value_r = estimate[length - 1];
         }

         send_flag_r++;
         MPI_Isend( &send_value_r, 1, braid_MPI_REAL, myid + 1 , 98, comm, &send_request_r );
      }
      
      /* Send situation 2: Send to the left if ilower -1 is not a C point */
      if ( ilower <= last_cpoint && ilower > 0 && !_braid_IsCPoint( ilower - 1, cfactor ) )
      {
         /* If we have no c points, wait to recv before sending */
         if ( ncpoints == 0 )
         {
            recv_flag--;
            MPI_Wait( &recv_request, MPI_STATUS_IGNORE );
            send_value_l = recv_value; 
         }
         else
         {
            send_value_l = estimate[0];
         }

         send_flag_l++;
         MPI_Isend( &send_value_l, 1, braid_MPI_REAL, myid - 1, 98, comm, &send_request_l );
      }

      /* Update the values */
      estimate_index = 0;
      factor = 0;   
      for ( i = ilower; i <= iupper; i++ )
      {
         /* If this is an F point interval at end of proc ownership */
         if (estimate_index >= length )
         {
            if (recv_flag)
            {
               recv_flag--;
               MPI_Wait( &recv_request , MPI_STATUS_IGNORE );
            }       
            factor = recv_value; 
         }
         else
         {
            factor = estimate[ estimate_index ] ;
         }

         /* Update the factor */
         dtk = pow( ta[i-ilower] - ta[i-ilower-1] , order );
         error_est[ i-ilower ] = factor * dtk ;
         
         /* increase the index at the Cpoints, unless this is the
          * last C point. In that case, use the same index */ 
         if ( _braid_IsCPoint( i, cfactor ) && i != last_cpoint )
         {
            estimate_index++;
         }
      }

      /* All recvieves should be completed. Finish up sends */
      if ( send_flag_r )
      {
         MPI_Wait( & send_request_r, MPI_STATUS_IGNORE );
      }
      if ( send_flag_l )
      {
         MPI_Wait( & send_request_l, MPI_STATUS_IGNORE );
      }

   }
   
   return _braid_error_flag;
}

braid_Int 
_braid_GetDtk(braid_Core core )
{
  
   _braid_Grid  **grids       = _braid_CoreElt(core, grids);
   _braid_Grid   *grid        = grids[0];
   MPI_Comm      comm         = _braid_CoreElt(core, comm);
   braid_Int     myid         = _braid_CoreElt(core, myid);
   braid_Int     order        = _braid_CoreElt(core, order); 
   braid_Int     ilower       = _braid_GridElt(grid, ilower);
   braid_Int     iupper       = _braid_GridElt(grid, iupper);
   braid_Int     gupper       = _braid_GridElt(grid, gupper);
   braid_Int     ncpoints     = _braid_GridElt(grid, ncpoints);
   braid_Int     cfactor      = _braid_GridElt(grid, cfactor);
   braid_Real   *dtk_core     = _braid_CoreElt(core, dtk ); 
   braid_Real   *ta           = _braid_GridElt(grid, ta);
   braid_Int fhi,flo,ci,interval,send_flag,recv_flag,i;
   braid_Real send_val, recv_val, dtk;
   MPI_Request recv_left, send_right ;
   
   int comm_size;
   MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
   
   if ( ilower <= iupper ) 
   {  
      send_flag = recv_flag = 0;
      send_val = recv_val = 0.0;
      
      if ( !_braid_IsCPoint( ilower - 1, cfactor) && ilower > 0 )
      {
         recv_flag = 1;
         MPI_Irecv( &recv_val, 1, braid_MPI_REAL, myid -1 , 56, comm, &recv_left );
      }

      braid_Int cpoint = ncpoints-1;
      for (interval = ncpoints; interval > -1; interval--)
      {
         _braid_GetInterval(core, 0, interval, &flo, &fhi, &ci);             
        
         dtk = 0;

         /* Wait for the value from the left */
         if ( interval == 0 && recv_flag )
         {
            MPI_Wait( &recv_left, MPI_STATUS_IGNORE);
            dtk = recv_val;

            recv_flag = 0;
         }
         
         //Finish up the interval 
         for ( i = flo; i <= fhi; i++ ) 
         {
            dtk += pow( ta[i-ilower] - ta[i-ilower-1] , order );
         }   
         /* Save dtk at the C points or send it on to the right is 
          * last point is not a C point */
         if ( ci > 0 ) 
         {
            dtk += pow( ta[ci-ilower] - ta[ci-ilower-1] , order );
            dtk_core[cpoint] = dtk;
            cpoint--;
         }
         
         if ( iupper != gupper && interval == ncpoints && !_braid_IsCPoint(iupper, cfactor) )
         {
            send_val  = dtk;
            send_flag = 1;
            MPI_Isend( &send_val, 1, braid_MPI_REAL, myid + 1, 56, comm, &send_right );
         }
      }   
      /* All recvieves should be complete. Finish up sends, and exit */
      if ( send_flag )
      {
         MPI_Wait( &send_right , MPI_STATUS_IGNORE);
      }
   }

   return _braid_error_flag;
}


/*----------------------------------------------------------------------------
 * Print the residual norm at ever C-point for debugging purposes 
 *----------------------------------------------------------------------------*/

braid_Int
_braid_PrintSpatialNorms(braid_Core    core,
                         braid_Real   *rnorms,     /* This processor's local residual norms at C-points */
                         braid_Int     n)          /* Length of the rnorms array */
{
   MPI_Comm       comm       = _braid_CoreElt(core, comm);
   MPI_Comm       comm_world = _braid_CoreElt(core, comm_world);
   _braid_Grid  **grids      = _braid_CoreElt(core, grids);
   braid_Int      cfactor    = _braid_GridElt(grids[0], cfactor);
   braid_Int      gupper     = _braid_CoreElt(core, gupper);

   braid_Int      g_ncpoints = ceil( ((braid_Real) (gupper + 1) ) / ((braid_Real) cfactor ));
   
   braid_Int     *recvcounts = NULL;
   braid_Real    *recvbuf;
   braid_Int     *displs;
   braid_Int      i, p, myid_t, myid_world, my_root_global_rank;

   MPI_Comm_size(comm, &p);
   MPI_Comm_rank(comm_world, &myid_world);
   MPI_Comm_rank(comm, &myid_t);

   /* We need to know all the processor's belonging to the temporal
    * communicator of global rank 0.  It is only these processors that are
    * involved with printing output. */
   my_root_global_rank = myid_world;
   MPI_Bcast(&my_root_global_rank, 1, braid_MPI_INT, 0, comm);

   if(my_root_global_rank == 0)
   {
      if(myid_t == 0)
      {
         recvbuf = _braid_CTAlloc(braid_Real, g_ncpoints);
         recvcounts = _braid_CTAlloc(braid_Int, p);
         displs = _braid_CTAlloc(braid_Int, p);
      }

      /* Rank 0 gather's every processor's number of C-points, which forms the
       * displacements (displs) for the Gatherv call below */ 
      MPI_Gather(&n, 1, braid_MPI_INT, recvcounts, 1, braid_MPI_INT, 0, comm);

      if(myid_t == 0)
      {
         displs[0] = 0;
         for(i = 1; i < p; i++)
         {
            displs[i] = displs[i-1] + recvcounts[i-1];
         }
      }

      /* Gather over comm */
      MPI_Gatherv(rnorms, n, braid_MPI_REAL, recvbuf, recvcounts, displs,
                  braid_MPI_REAL, 0, comm);

      if(myid_t == 0)
      {
         for(i = 0; i <= g_ncpoints; i++){
            _braid_printf("  Braid:  time step: %6d, rnorm: %1.2e\n", i*cfactor, recvbuf[i]);
         }
      }

      if(myid_t == 0)
      {
         _braid_TFree(recvbuf);
         _braid_TFree(recvcounts);
         _braid_TFree(displs);
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Copy the initialized C-points on the fine grid, to all coarse levels.
 * Allows first down cycle to be skipped, in FMG fashion.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_CopyFineToCoarse(braid_Core  core)
{
   braid_App      app     = _braid_CoreElt(core, app);
   _braid_Grid  **grids   = _braid_CoreElt(core, grids);
   braid_Int      nlevels = _braid_CoreElt(core, nlevels);
   
   braid_Int      f_index, index, iu, is_stored, level, f_cfactor;
   braid_Int      ilower, iupper;
   braid_Vector   u, *va;

   for(level = 1; level < nlevels; level++)
   {

      f_cfactor = _braid_GridElt(grids[level-1], cfactor);
      iupper    = _braid_GridElt(grids[level], iupper);
      ilower    = _braid_GridElt(grids[level], ilower);
      va        = _braid_GridElt(grids[level], va);

      /* Loop over all points belonging to this processor, and if a C-point,
       * then carry out spatial coarsening and copy to ua and va */
      for (index=ilower; index<=iupper; index++)
      {
         _braid_MapCoarseToFine(index, f_cfactor, f_index);
         _braid_UGetVector(core, level-1, f_index, &u);
         _braid_Coarsen(core, level, f_index, index, u, &va[index-ilower]);
         
         _braid_CoreFcn(core, free)(app, u);
         _braid_CoreFcn(core, clone)(app, va[index-ilower], &u);
         _braid_USetVectorRef(core, level, index, u);
         _braid_UGetIndex(core, level, index, &iu, &is_stored);
         if (is_stored == -2) /* Case where F-points are not stored, and we are not using shell vectors */
         {
            _braid_CoreFcn(core, free)(app, u);
         }
         else if (is_stored == -1) /* This is a shell vector */
         {
            // We free the data in u, keeping the shell
            _braid_CoreFcn(core, sfree)(app, u);
         }
 
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Set the residual norm for iteration iter.  If iter < 0, set the rnorm for the
 * last iteration minus |iter|-1.  Also set the initial residual norm.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_SetRNorm(braid_Core  core,
                braid_Int   iter,
                braid_Real  rnorm)
{
   braid_Real  *rnorms   = _braid_CoreElt(core, rnorms);
   braid_Int    max_iter = _braid_CoreElt(core, max_iter);
   braid_Int    skip     = _braid_CoreElt(core, skip);
   braid_Int    k;

   k = iter;
   if (iter < 0)
   {
      k = _braid_CoreElt(core, niter) + 1 + iter;
   }

   if ((k > -1) && (k <= max_iter)) 
   {
      rnorms[k] = rnorm;

      /* Set initial residual norm if not already set */
      if ((k == 0) || ((k==1) && (skip == 1)))
      {
         if ( _braid_CoreElt(core, rnorm0) == braid_INVALID_RNORM )
         {
            _braid_CoreElt(core, rnorm0) = rnorm;
         }
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Get the residual norm for iteration iter.  If iter < 0, get the rnorm for the
 * last iteration minus |iter|-1.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetRNorm(braid_Core  core,
                braid_Int   iter,
                braid_Real *rnorm_ptr)
{
   braid_Real  *rnorms   = _braid_CoreElt(core, rnorms);
   braid_Int    max_iter = _braid_CoreElt(core, max_iter);
   braid_Int    k;

   /* Initialize to invalid value */
   *rnorm_ptr = braid_INVALID_RNORM;

   k = iter;
   if (iter < 0)
   {
      k = _braid_CoreElt(core, niter) + 1 + iter;
   }

   if ((k > -1) && (k <= max_iter)) 
   {
      *rnorm_ptr = rnorms[k];
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Same as SetRNorm, but sets full residual norm
 *----------------------------------------------------------------------------*/

braid_Int
_braid_SetFullRNorm(braid_Core  core,
                    braid_Int   iter,
                    braid_Real  rnorm)
{
   braid_Real  *rnorms   = _braid_CoreElt(core, full_rnorms);
   braid_Int    max_iter = _braid_CoreElt(core, max_iter);
   braid_Int    skip     = _braid_CoreElt(core, skip);
   braid_Int    k;

   if (rnorms == NULL)
   {
      return _braid_error_flag;
   }

   k = iter;
   if (iter < 0)
   {
      k = _braid_CoreElt(core, niter) + 1 + iter;
   }

   if ((k > -1) && (k <= max_iter)) 
   {
      rnorms[k] = rnorm;

      /* Set initial residual norm if not already set */
      if ((k == 0) || ((k==1) && (skip == 1)))
      {
         if ( _braid_CoreElt(core, full_rnorm0) == braid_INVALID_RNORM )
         {
            _braid_CoreElt(core, full_rnorm0) = rnorm;
         }
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Same as GetRNorm, but gets full residual norm
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetFullRNorm(braid_Core  core,
                    braid_Int   iter,
                    braid_Real *rnorm_ptr)
{
   braid_Real  *rnorms   = _braid_CoreElt(core, full_rnorms);
   braid_Int    max_iter = _braid_CoreElt(core, max_iter);
   braid_Int    k;

   /* Initialize to invalid value */
   *rnorm_ptr = braid_INVALID_RNORM;

   if (rnorms == NULL)
   {
      return _braid_error_flag;
   }

   k = iter;
   if (iter < 0)
   {
      k = _braid_CoreElt(core, niter) + 1 + iter;
   }

   if ((k > -1) && (k <= max_iter)) 
   {
      *rnorm_ptr = rnorms[k];
   }

   return _braid_error_flag;
}

/*---------------------------------------------------------------------------
 * Allow access to the global access function following a Vcycle
 *---------------------------------------------------------------------------*/
braid_Int
_braid_GlobalAccess(braid_Core   core )
{

   /* Set the calling function */
   _braid_CoreElt(core, calling_function) = braid_Caller_GlobalAccess;
   
   if ( _braid_CoreElt( core, gaccess ) != NULL )
   {
      MPI_Comm      comm       = _braid_CoreElt(core, comm);
      braid_App     app        = _braid_CoreElt(core, app);
      _braid_Grid  **grids     = _braid_CoreElt(core, grids);
      braid_Int      ilower    = _braid_GridElt(grids[0], ilower);
      braid_Int      iupper    = _braid_GridElt(grids[0], iupper);
      braid_Int i, index;

      braid_Int npoints = iupper-ilower+1;
      braid_Vector  *ua = _braid_CTAlloc(braid_Vector, npoints);
      braid_Int *umask = _braid_CTAlloc(braid_Int, npoints);
      
      for ( i = ilower; i <= iupper; i++ )
      {
         _braid_UGetIndex(core, 0, i, &index, &umask[i-ilower]);
         _braid_UGetVectorRef(core, 0 , i , &ua[i-ilower] );
         if ( ua[i-ilower] == NULL )
         {
            umask[i-ilower] = -2;
         }
      }

      braid_GlobalAccessStatus gstatus   = (braid_GlobalAccessStatus)core;
      _braid_CoreFcn(core, gaccess)(app, comm, ua, umask, npoints, gstatus );

      _braid_TFree(umask);
      _braid_TFree(ua);
   }

   return _braid_error_flag;
}


/*----------------------------------------------------------------------------
 * Process the error with code ierr raised at the given line and source file
 *----------------------------------------------------------------------------*/
void
_braid_ErrorHandler(const char *filename,
                    braid_Int   line,
                    braid_Int   ierr,
                    const char *msg)
{
   _braid_error_flag |= ierr;

#ifdef braid_PRINT_ERRORS
   if (msg)
   {
      _braid_printf("braid error in file \"%s\", line %d, error code = %d - %s\n",
                    filename, line, ierr, msg);
   }
   else
   {
      _braid_printf("braid error in file \"%s\", line %d, error code = %d\n",
                    filename, line, ierr);
   }
#endif
}

