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

/** \file _braid_dist.h
 * \brief Header files for the data distribution functions and load balancing.
 */
#ifndef _braiddist_HEADER
#define _braiddist_HEADER
#include "_braid.h"
#include "_util.h"

/* The balance structure contains all the information needed to complete 
 * the new load balancing algorithm. It gets passed around until it has
 * been fully built. Then it can be passed to FRefine, to complete the 
 * load balancing and temporal refinement. 
 */

/* Experimental load balancing threshold. If this is of use, it could go into the 
 * core as "another" flag. The idea is that temporal load balancing is not nessesary
 * when only a small number of time steps are refined. For example, under the current 
 * scenario, refining the very first time step, and no other time steps, could potentially
 * require load balancing on every processor. The downside if that if we do this, then the 
 * partition is no longer uniform, so the assumed partition needs to be called. */
#define LOAD_BAL_THRES_EXPERIMENTAL -1

typedef struct
{
    //Set by InitBalanceStruct
    braid_Int  refine;          /* bool -- refinment in time */
    braid_Int  lbalance;        /* bool -- weigted load balance */ 
    
    braid_Int  coarse_ilower;   /* coarse grid ilower */
    braid_Int  coarse_iupper;   /* coarse grid iupper */
    braid_Int  coarse_gupper;   /* coarse grid gupper */
    
    // Set inside Get Refined Interval
    braid_Int  refined_ilower;  /* refined grid ilower */
    braid_Int  refined_iupper;  /* refined grid iupper */
    braid_Int  refined_gupper;  /* refined grid gupper */
    braid_Int  refined_npoints; /* number of points refined grid */   
    braid_Int  *refined_ca;           /* map from refined to coarse grid */
    braid_Int  *refined_fa;           /* index map from coarse to fine grid */ 
    braid_Real *refined_taalloc;     /* allocation for the refined time grid */
    braid_Real *refined_ta;           /* refined time grid (shifted from alloc) */
     
    //Set Inside 
    braid_Int  fine_ilower;     /* fine grid ilower */
    braid_Int  fine_iupper;     /* fine grid iupper */
    braid_Int  fine_gupper;     /* fine grid igupper */
    braid_Int  *fine_ca; 
    braid_Real *fine_ta;
    braid_Real *fine_ta_alloc; 
    braid_Int  fine_first;
    braid_Int  fine_next;
    
    //Communication Maps
    braid_Int  *right_procs;    /* maps to processor owning ilower -1 on lvel l. */
    braid_Int  *left_procs;     /* maps to processor owning iupper +1 on lvel l */
    braid_Int  *send_map_alloc; /* send_map allocation. */
    braid_Int  *send_map;       /* load balancing map */
    braid_Int  *recv_map;       /* load balancing recv map */ 

} _braid_BalanceStruct;



/**
 * Main function to build the structure. Calls GetRefined
 * interval, get Fine interval, and build comm map 
 **/
braid_Int 
_braid_BuildBalanceStructure(braid_Core            core,
                             braid_Int            *done,
                            _braid_BalanceStruct *bstruct);

/**
 * Build Phase 1: Initialize the Balance structure
 **/
braid_Int
_braid_InitBalanceStruct(_braid_BalanceStruct *bstruct,
                         braid_Int             refine,
                         braid_Int             lbalance,
                         braid_Int             coarse_ilower,
                         braid_Int             coarse_iupper,
                         braid_Int             coarse_gupper);


/** 
 * Buld Phase 2: This function takes the rfactors, and gets 
 * the new ilower and iupper for the refined grid. Also, if
 * load balancing will be completed, this function fixes the
 * wfactors array so that it matches the new refined grid. 
 **/
braid_Int
_braid_GetRefinedInterval(braid_Core            core,
                          braid_Int            *done,
                         _braid_BalanceStruct *bstruct);

/** 
 * Buld Phase 2b: Get the mapping for the indices of the 
 * refined grid **/
braid_Int
_braid_GetRefinedGridMap(braid_Core            core,
                         braid_Int            *done,
                         _braid_BalanceStruct *bstruct);

/**
 * Build Phase 3: Get the new ilower and iupper for the 
 * fine grid. 
 **/
braid_Int
_braid_GetFineInterval(braid_Core           core,
                       braid_Int           *done,
                      _braid_BalanceStruct *bstruct );

/**
 * Build Phase 3: Get the new ilower and iupper for the 
 * fine grid. 
 **/
braid_Int
_braid_GetFineGridMap(braid_Core           core,
                       braid_Int           *done,
                      _braid_BalanceStruct *bstruct );

/**
 * Build Phase 4: Build the required communication maps 
 **/
braid_Int
_braid_BuildCommunicationMap(braid_Core            core,
                             braid_Int            *done,
                             _braid_BalanceStruct *bstruct);
/**
 * Destroy the balance structure 
 **/
braid_Int
_braid_DestroyBalanceStruct(_braid_BalanceStruct *bstruct);


/**************************************************************************
 * Block Data distribution functions. These functions use an O(1) global
 * function to distribute time points evenly across the available time
 * processors. This is the defualt distribution used by xbraid, and is also
 * the distribution used to calculate the assumed partition.
 **************************************************************************/

/**
 * Returns the index interval for *proc* in a blocked data distribution.
 **/
braid_Int
_braid_GetBlockDistInterval(braid_Core  core,
                            braid_Int   npoints,
                            braid_Int  *ilower_ptr,
                            braid_Int  *iupper_ptr);

/**
 * Helper function for GetBlockDistInterval.
 **/
braid_Int
_braid_GetBlockDistInterval_basic(braid_Int   npoints,
                                  braid_Int   nprocs,
                                  braid_Int   proc,
                                  braid_Int  *ilower_ptr,
                                  braid_Int  *iupper_ptr);

/**
 * Returns the processor that owns *index* in a blocked data distribution
 * (returns -1 if *index* is out of range).
 **/

braid_Int
_braid_GetBlockDistProc(braid_Core  core,
                        braid_Int   npoints,
                        braid_Int   index,
                        braid_Int  *proc_ptr);

/**
 * Helper Function for the GetBlockDistProc function
 **/

braid_Int
_braid_GetBlockDistProc_basic(braid_Int   npoints,
                              braid_Int   nprocs,
                              braid_Int   index,
                              braid_Int  *proc_ptr);


/**************************************************************************
 * Weighted Data distribution functions. These functions use the user defined
 * weights to distribute time points across the interval. The assumed partition
 * algorithm is used to minimize global communication and storage.
 **************************************************************************/

/** _braid_WeightedStruct is used to pass information about the weights to the corresponding
 * functions. This is used for simplicicty, as it allows all communication needed during
 * the load balancing and temporal refinment steps to be completed at the same time.
 **/

#define _braid_WeightedElt( weighted, elt)  ( (weighted) -> elt )
typedef struct
{
    braid_Real local_max;   /* local max of the weights */
    braid_Real local_min;   /* local min of the weights */
    braid_Real local_sum;   /* local sum of the weights */
    braid_Real global_max;  /* global max of the weights */
    braid_Real global_min;  /* global min of the weights */
    braid_Real global_sum;  /* global sum of the weights */
    braid_Real local_start; /* cumulative sum of the weights exclusive */
    braid_Real local_stop;  /* cumulative sum of the weights inclusive */

} _braid_WeightedStruct;

/**
 * Init a _braid_WeightedStruct structure
 **/
braid_Int
_braid_WeightedStructInit( _braid_WeightedStruct *wstruct );


/**
 * Destroy a _braid_WeightedStruct structure
 **/
braid_Int
_braid_WeightedStructDestroy( _braid_WeightedStruct *wstruct );

/**
 * Build the Balance structure required to load balance on a (possibly refined )
 * new grid, using a weighted distribution.
 */
braid_Int
_braid_BuildWeightedStruct(braid_Core             core,
                           braid_Int             *done,
                          _braid_BalanceStruct  *bstruct, 
                          _braid_WeightedStruct *wstruct);

/**
 * Gets the new fine grid interval based on the user defined weights
 * and stores that information in the Balance structure.
 */
braid_Int
_braid_GetWeightedFineInterval(braid_Core            core,
                               braid_Int            *done,
                              _braid_BalanceStruct  *bstruct );


/**
 * Chooses an assumed partition based on a global O(1) function using the BlockDist
 * functions, and populates an array, owners, with the processor numbers that own the
 * corresponding indicies in the assumed partition.
 */
braid_Int
_braid_GetPartition(braid_Core    core,
                    braid_Int     ilower,
                    braid_Int     iupper,
                    braid_Int     gupper,
                    braid_Int    *assumed_ilower,
                    braid_Int    *assumed_iupper,
                    braid_Int   **owners);

/** 
 * User defined MPI function used to reduce communication. Function takes
 * a double of length 5 and returns the sum,max,min,max,min. 
 */
void 
Sum3MaxMin(braid_Real *in, 
           braid_Real *inout, 
           braid_Int *len, 
           MPI_Datatype *datatype);

#endif
