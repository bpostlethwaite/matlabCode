/* oneProjector_core.c
   $Id: oneProjector_core.c 272 2007-07-19 05:40:29Z mpf $

   ----------------------------------------------------------------------
   This file is part of SPGL1 (Spectral Projected Gradient for L1).

   Copyright (C) 2007 Ewout van den Berg and Michael P. Friedlander,
   Department of Computer Science, University of British Columbia, Canada.
   All rights reserved. E-mail: <{ewout78,mpf}@cs.ubc.ca>.

   SPGL1 is free software; you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1 of the
   License, or (at your option) any later version.

   SPGL1 is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General
   Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with SPGL1; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
   USA
   ----------------------------------------------------------------------
*/

#include "oneProjector_core.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>     /* provides DBL_EPSILON */
#include <sys/types.h>



/* ======================================================================= */
/*                     H E L P E R   F U N C T I O N S                     */
/* ======================================================================= */

#define swap_double(a,b) { double c; c = (a); (a) = (b); (b) = c; }

/*!
  \brief Perform the "sift" operation for the heap-sort algorithm.

  A heap is a collection of items arranged in a binary tree.  Each
  child node is greater than or equal to its parent.  If x[k] is the
  parent, than its children are x[2k+1] and x[2k+2].

  This routine promotes ("sifts up") children that are smaller than
  their parents.  Thus, this is a "inverse" heap, where the smallest
  element of the heap is the root node.

  \param[in]     root       The root index from which to start sifting.
  \param[in]     lastChild  The last child (largest node index) in the sift operation.
  \param[in,out] x          The array to be sifted.
*/
static void heap_sift( int root, int lastChild, double x[] )
{
    int child;

    for (; (child = (root * 2) + 1) <= lastChild; root = child) {

	if (child < lastChild)
	    if ( fabs(x[child]) < fabs(x[child+1]) )
		child++;
	
	if ( fabs(x[child]) <= fabs(x[root]) )
	    break;

	swap_double(  x[root],  x[child] );
    }
}


/*!
  \brief Perform the "sift" operation for the heap-sort algorithm.

  A heap is a collection of items arranged in a binary tree.  Each
  child node is greater than or equal to its parent.  If x[k] is the
  parent, than its children are x[2k+1] and x[2k+2].

  This routine promotes ("sifts up") children that are smaller than
  their parents.  Thus, this is a "inverse" heap, where the smallest
  element of the heap is the root node.

  Elements in y are associated with those in x and are reordered accordingly.

  \param[in]     root       The root index from which to start sifting.
  \param[in]     lastChild  The last child (largest node index) in the sift operation.
  \param[in,out] x          The array to be sifted.
  \param[in,out] y          The array to be sifted accordingly.
*/
static void heap_sift_2( int root, int lastChild, double x[], double y[] )
{
    int child;

    for (; (child = (root * 2) + 1) <= lastChild; root = child) {

	if (child < lastChild)
	    if ( fabs(x[child]) < fabs(x[child+1]) )
		child++;
	
	if ( fabs(x[child]) <= fabs(x[root]) )
	    break;

	swap_double( x[root], x[child] );
        swap_double( y[root], y[child] );
    }
}


/*!
  \brief Discard the smallest element and contract the heap.

  On entry, the numElems of the heap are stored in x[0],...,x[numElems-1],
  and the biggest element is x[0].  The following operations are performed:
    -# Swap the first and last elements of the heap
    -# Shorten the length of the heap by one.
    -# Restore the heap property to the contracted heap.
       This effectively makes x[0] the next smallest element
       in the list.  

  \param[in]     numElems   The number of elements in the current heap.
  \param[in,out] x          The array to be modified.

  \return  The number of elements in the heap after it has been contracted.
*/
static int heap_del_min(int numElems, double x[])
{
    int lastChild = numElems - 1;

    assert(numElems > 0);

    /* Swap the smallest element with the lastChild. */
    swap_double(x[0], x[lastChild]);

    /* Contract the heap size, thereby discarding the smallest element. */
    lastChild--;
    
    /* Restore the heap property of the contracted heap. */
    heap_sift(0, lastChild, x);

    return numElems - 1;
}


/*!
  \brief Discard the smallest element of x and contract the heaps.

  On entry, the numElems of the heap are stored in x[0],...,x[numElems-1],
  and the smallest element is x[0].  The following operations are performed:
    -# Swap the first and last elements of both heaps
    -# Shorten the length of the heaps by one.
    -# Restore the heap property to the contracted heap x.
       This effectively makes x[0] the next smallest element
       in the list.  

  \param[in]     numElems   The number of elements in the current heap.
  \param[in,out] x          The array to be modified.
  \param[in,out] y          The array to be modified accordingly

  \return  The number of elements in each heap after they have been contracted.
*/
static int heap_del_min_2( int numElems, double x[], double y[] )
{
    int lastChild = numElems - 1;

    assert(numElems > 0);

    /* Swap the smallest element with the lastChild. */
    swap_double( x[0], x[lastChild] );
    swap_double( y[0], y[lastChild] ); 

    /* Contract the heap size, thereby discarding the smallest element. */
    lastChild--;
    
    /* Restore the heap property of the contracted heap. */
    heap_sift_2( 0, lastChild, x, y );

    return numElems - 1;
}


/*!
  
  \brief  Build a heap by adding one element at a time.
  
  \param[in]      n   The length of x and ix.
  \param[in,out]  x   The array to be heapified.

*/
static void heap_build( int n, double x[] )
{    
    int i;

    for (i = n/2; i >= 0; i--) heap_sift( i, n-1, x );
}


/*!
  
  \brief  Build a heap by adding one element at a time.
  
  \param[in]      n   The length of x and ix.
  \param[in,out]  x   The array to be heapified.
  \param[in,out]  y   The array to be reordered in sync. with x.

*/
static void heap_build_2( int n, double x[], double y[] )
{
    int i;

    for (i = n/2; i >= 0; i--) heap_sift_2( i, n-1, x, y );
}




/* ----------------------------------------------------------------------- */
int projectI(double xPtr[], double bPtr[], double lambda, int n)
/* ----------------------------------------------------------------------- */
{  int
       i, j;
   double
       b,          /* Current element of vector b */
       csb,        /* Cumulative sum of b */
       alpha = 0,
       tau   = 0;

   /* The vector xPtr[] is initialized to bPtr[] prior to the function call */

   /* Check if lambda is essentially zero.  Exit with x = 0. */
   if (lambda < DBL_EPSILON) {
       for (i = 0; i < n; i++) xPtr[i] = 0;
       return 0;
   }

   /* Check if ||b||_1 <= lambda.  Exit with x = b. */
   for (csb = 0, i = 0; i < n; i++) csb += fabs(bPtr[i]);
   if (csb <= lambda)
       return 0;

   /* Set up the heap */
   heap_build(n, xPtr);

   /* Initialise csb with -lambda so we don't have to subtract this at every iteration */
   csb = -lambda;

   /* Determine threshold value tau */
   for (i = n, j = 0; j < n; tau = alpha)
   {  
      b = fabs(xPtr[0]); /* Get current minimum                      */
      j ++;              /* Give compiler some room for optimization */
      csb += b;          /* Update the cumulative sum of b           */

      /* Move heap to next minimum value */
      i = heap_del_min(i, xPtr);

      /* Compute the required step to satisfy the lambda constraint */
      alpha  = csb / j;

      /* We are done as soon as the constraint can be satisfied */
      /* without exceeding the current minimum value of b       */
      if (alpha >= b)
          break;
   }

   /* Set the solution by applying soft-thresholding with tau */
   for (i = 0; i < n; i++)
   {  b = bPtr[i];
      if (fabs(b) <= tau)
           xPtr[i] = 0;
      else xPtr[i] = b - tau * (b < 0 ? -1 : 1); 
   }

   return j;
}


/* ----------------------------------------------------------------------- */
int projectD(double xPtr[], double bPtr[], double dPtr[], double lambda, int n)
/* ----------------------------------------------------------------------- */
{  int
       i, j;
   double
       csdb,        /* Cumulative sum of d.*b          */
       csd2,        /* Cumulative sum of d.^2          */
       b,           /* Current element of vector b     */
       d,           /* Current element of vector d     */
       bd,          /* Current element of vector b / d */
       alpha  = 0,
       tau    = 0;

   /* Check if lambda is essentially zero.  Exit with x = 0. */
   if (lambda < DBL_EPSILON)
   {   for (i = 0; i < n; i++) xPtr[i] = 0;
       return 0;
   }

   /* Preliminary check on trivial solution x = Db (and scale b) */
   for (csdb = 0, i = 0; i < n; i++)
   {  d = dPtr[i];
      b = xPtr[i];
      csdb += (d * fabs(b));
      xPtr[i] = b / d;
   }
   if (csdb <= lambda)
   {  
       /* Reset the entries of x to b */
      memcpy((void *)xPtr, (void *)bPtr, n * sizeof(double));
      return 0;
   }

   /* Set up the heap (we have to sort on Db) */
   heap_build_2(n, xPtr, dPtr);

   /* Initialise csbd with -lambda so we don't have to subtract this at every iteration */
   csdb = -lambda;
   csd2 =  0;

   /* Determine the threshold level tau */
   for (i = n, j = 0; j < n; tau = alpha)
   {
      bd    = fabs(xPtr[0]);  /* Get current minimum b / d                */
      j    ++;                /* Give compiler some room for optimization */
      d     = dPtr[0];        /* Get current value of d                   */
      d    *= d;              /* Compute d squared                        */
      csd2 += d;              /* Update the cumulative sum of d.*d        */
      csdb += bd * d;         /* Update the cumulative sum of d.*b        */

      /* Move heap to next minimum value */
      i = heap_del_min_2(i, bPtr, dPtr);

      /* Compute the required step to satisfy the lambda constraint */
      alpha  = csdb / csd2;

      /* We are done as soon as the constraint can be satisfied */
      /* without exceeding the current minimum value of b / d   */
      if (alpha >= bd) break;
   }

   /* Set the solution */
   for (i = 0; i < n; i++)
   {  b     = bPtr[i];
      alpha = dPtr[i] * tau;
      if (fabs(b) <= alpha)
           xPtr[i] = 0;
      else xPtr[i] = b - alpha * (b < 0 ? -1 : 1);
      
   }

   return j;
}
