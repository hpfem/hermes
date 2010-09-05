/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*             ********   ***                                 SparseLib++    */
/*          *******  **  ***       ***      ***                              */
/*           *****      ***     ******** ********                            */
/*            *****    ***     ******** ********              R. Pozo        */
/*       **  *******  ***   **   ***      ***                 K. Remington   */
/*        ********   ********                                 A. Lumsdaine   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                     SparseLib++ : Sparse Matrix Library                   */
/*                                                                           */
/*               National Institute of Standards and Technology              */
/*                        University of Notre Dame                           */
/*              Authors: R. Pozo, K. Remington, A. Lumsdaine                 */
/*                                                                           */
/*                                 NOTICE                                    */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that the above notice appear in all copies and supporting        */
/* documentation.                                                            */
/*                                                                           */
/* Neither the Institutions (National Institute of Standards and Technology, */
/* University of Notre Dame) nor the Authors make any representations about  */
/* the suitability of this software for any purpose.  This software is       */
/* provided ``as is'' without expressed or implied warranty.                 */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// An adaptation of Schmidt's new quicksort (taken from GNU emacs)
// for MV++
// by Andrew Lumsdaine

#include "qsort_int.h"

#ifdef sparc
#include <alloca.h>
#endif

// Invoke the < comparison function
#define CMP(A,B) ((A) < (B))

// swap two items
//
static inline void SWAP_int(int &A, int &B)
{
    int tmp = A; A = B; B = tmp;
}
 

static inline void swap_int(int &a, int &b)
{
    int tmp=a; a=b; b=tmp;
}

// This should be replaced by a standard ANSI macro.
#define BYTES_PER_WORD 8

/* The next 4 #defines implement a very fast in-line stack abstraction. */
#define STACK_SIZE (BYTES_PER_WORD * sizeof (long))
#define PUSH(LOW,HIGH) do {top->lo = LOW;top++->hi = HIGH;} while (0)
#define POP(LOW,HIGH)  do {LOW = (--top)->lo;HIGH = top->hi;} while (0)
#define STACK_NOT_EMPTY (stack < top)                

/* Discontinue quicksort algorithm when partition gets below this size.
   This particular magic number was chosen to work best on a Sun 4/260. */
#define MAX_THRESH 4

/* Stack node declarations used to store unfulfilled partition obligations. */
typedef struct {
  int lo;
  int hi;
} stack_node;

/* Order size using quicksort.  This implementation incorporates
   four optimizations discussed in Sedgewick:
   
   1. Non-recursive, using an explicit stack of pointer that store the 
      next array partition to sort.  To save time, this maximum amount 
      of space required to store an array of MAX_INT is allocated on the 
      stack.  Assuming a 32-bit integer, this needs only 32 * 
      sizeof (stack_node) == 136 bits.  Pretty cheap, actually.

   2. Chose the pivot element using a median-of-three decision tree.
      This reduces the probability of selecting a bad pivot value and 
      eliminates certain extraneous comparisons.

   3. Only quicksorts TOTAL_ELEMS / MAX_THRESH partitions, leaving
      insertion sort to order the MAX_THRESH items within each partition.  
      This is a big win, since insertion sort is faster for small, mostly
      sorted array segments.
   
   4. The larger of the two sub-partitions is always pushed onto the
      stack first, with the algorithm then concentrating on the
      smaller partition.  This *guarantees* no more than log (n)
      stack size is needed (actually O(1) in this case)! */
      

int QSort(VECTOR_int& v, int base_ptr, int total_elems)
{
  int pivot_buffer;
  
  if (total_elems > MAX_THRESH) {
    
    int lo = base_ptr;
    int hi = lo + total_elems - 1;
    
    stack_node stack[STACK_SIZE]; /* Largest size needed for 32-bit int!!! */
    stack_node *top = stack + 1;
    
    while (STACK_NOT_EMPTY) {
      int left_ptr;
      int right_ptr;
      {

    {

      /* Select median value from among LO, MID, and HI. Rearrange
         LO and HI so the three values are sorted. This lowers the 
         probability of picking a pathological pivot value and 
         skips a comparison for both the LEFT_PTR and RIGHT_PTR. */
      
      int mid = lo + (hi - lo) / 2;
      
      if (CMP (v[mid], v[lo]))
        SWAP_int(v[mid], v[lo]);
      if (CMP (v[hi], v[mid]))
        SWAP_int(v[hi], v[mid]);
      else 
        goto jump_over;
      
      if (CMP (v[mid], v[lo]))
        SWAP_int(v[mid], v[lo]);
      
    jump_over:

      pivot_buffer = v[mid];
    }
    
    left_ptr  = lo + 1;
    right_ptr = hi - 1;
    
    /* Here's the famous ``collapse the walls'' section of quicksort.  
       Gotta like those tight inner loops!  They are the main reason 
       that this algorithm runs much faster than others. */
    do {
      while (CMP (v[left_ptr], pivot_buffer))
        left_ptr++;
      
      while (CMP (pivot_buffer, v[right_ptr]))
        right_ptr--;
      
      if (left_ptr < right_ptr) {
        SWAP_int(v[left_ptr], v[right_ptr]);
        left_ptr++;
        right_ptr--;
      } else if (left_ptr == right_ptr) {
        left_ptr ++;
        right_ptr --;
        break;
      }
    } while (left_ptr <= right_ptr);
      }
      
      /* Set up pointers for next iteration.  First determine whether
     left and right partitions are below the threshold size. If so, 
     ignore one or both.  Otherwise, push the larger partition's
     bounds on the stack and continue sorting the smaller one. */
      
      if ((right_ptr - lo) <= MAX_THRESH) {
    if ((hi - left_ptr) <= MAX_THRESH)
      POP (lo, hi); 
    else
      lo = left_ptr;
      } else if ((hi - left_ptr) <= MAX_THRESH)
    hi = right_ptr;
      else if ((right_ptr - lo) > (hi - left_ptr)) {                   
    PUSH (lo, right_ptr);
    lo = left_ptr;
      } else {                   
    PUSH (left_ptr, hi);
    hi = right_ptr;
      }
    }
  }
  
  /* Once the BASE_PTR array is partially sorted by quicksort the rest
     is completely sorted using insertion sort, since this is efficient 
     for partitions below MAX_THRESH size. BASE_PTR points to the beginning 
     of the array to sort, and END_PTR points at the very last element in
     the array (*not* one beyond it!). */
  
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

  {
    int end_ptr = base_ptr + total_elems - 1;
    int run_ptr;
    int tmp_ptr = base_ptr;
    int thresh = MIN (end_ptr, base_ptr + MAX_THRESH);
    
    for (run_ptr = tmp_ptr + 1; run_ptr <= thresh; run_ptr++)
      if (CMP (v[run_ptr], v[tmp_ptr]))
    tmp_ptr = run_ptr;
    
    if (tmp_ptr != base_ptr)
      SWAP_int(v[tmp_ptr], v[base_ptr]);
    
    for (run_ptr = base_ptr + 1; (tmp_ptr = run_ptr += 1) <= end_ptr;) {
      
      while (CMP (v[run_ptr], v[tmp_ptr -= 1]))
    ;
      
      if ((tmp_ptr += 1) != run_ptr) {
    int trav;
    
    for (trav = run_ptr + 1; --trav >= run_ptr;) {
      int c;
      c = v[trav];
      int hi, lo;
      
      for (hi = lo = trav; (lo -= 1) >= tmp_ptr; hi = lo)
        v[hi] = v[lo];
      v[hi] = c;
    }
      }
    }
  }

  return 1;
}


int QSort(VECTOR_int & v, VECTOR_int& x, int base_ptr, int total_elems)
{
  int pivot_buffer;
  int pixot_buffer;
  
  if (total_elems > MAX_THRESH) {
    
    int lo = base_ptr;
    int hi = lo + total_elems - 1;
    
    stack_node stack[STACK_SIZE]; /* Largest size needed for 32-bit int!!! */
    stack_node *top = stack + 1;
    
    while (STACK_NOT_EMPTY) {
      int left_ptr;
      int right_ptr;
      {

    {

      /* Select median value from among LO, MID, and HI. Rearrange
         LO and HI so the three values are sorted. This lowers the 
         probability of picking a pathological pivot value and 
         skips a comparison for both the LEFT_PTR and RIGHT_PTR. */
      
      int mid = lo + (hi - lo) / 2;
      
      if (CMP (v[mid], v[lo])) {
        swap_int (v[mid], v[lo]);
        SWAP_int (x[mid], x[lo]);
      }
      if (CMP (v[hi], v[mid])) {
        swap_int (v[hi], v[mid]);
        SWAP_int (x[hi], x[mid]);
      } else 
        goto jump_over;
      
      if (CMP (v[mid], v[lo])) {
        swap_int(v[mid], v[lo]);
        SWAP_int (x[mid], x[lo]);
      }
      
    jump_over:

      pivot_buffer = v[mid];
      pixot_buffer = x[mid];
    }
    
    left_ptr  = lo + 1;
    right_ptr = hi - 1;
    
    /* Here's the famous ``collapse the walls'' section of quicksort.  
       Gotta like those tight inner loops!  They are the main reason 
       that this algorithm runs much faster than others. */
    do {
      while (CMP (v[left_ptr], pivot_buffer))
        left_ptr++;
      
      while (CMP (pivot_buffer, v[right_ptr]))
        right_ptr--;
      
      if (left_ptr < right_ptr) {
        swap_int (v[left_ptr], v[right_ptr]);
        SWAP_int (x[left_ptr], x[right_ptr]);
        left_ptr++;
        right_ptr--;
      } else if (left_ptr == right_ptr) {
        left_ptr ++;
        right_ptr --;
        break;
      }
    } while (left_ptr <= right_ptr);
      }
      
      /* Set up pointers for next iteration.  First determine whether
     left and right partitions are below the threshold size. If so, 
     ignore one or both.  Otherwise, push the larger partition's
     bounds on the stack and continue sorting the smaller one. */
      
      if ((right_ptr - lo) <= MAX_THRESH) {
    if ((hi - left_ptr) <= MAX_THRESH)
      POP (lo, hi); 
    else
      lo = left_ptr;
      } else if ((hi - left_ptr) <= MAX_THRESH)
    hi = right_ptr;
      else if ((right_ptr - lo) > (hi - left_ptr)) {                   
    PUSH (lo, right_ptr);
    lo = left_ptr;
      } else {                   
    PUSH (left_ptr, hi);
    hi = right_ptr;
      }
    }
  }
  
  /* Once the BASE_PTR array is partially sorted by quicksort the rest
     is completely sorted using insertion sort, since this is efficient 
     for partitions below MAX_THRESH size. BASE_PTR points to the beginning 
     of the array to sort, and END_PTR points at the very last element in
     the array (*not* one beyond it!). */
  
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

  {
    int end_ptr = base_ptr + total_elems - 1;
    int run_ptr;
    int tmp_ptr = base_ptr;
    int thresh = MIN (end_ptr, base_ptr + MAX_THRESH);
    
    for (run_ptr = tmp_ptr + 1; run_ptr <= thresh; run_ptr++)
      if (CMP (v[run_ptr], v[tmp_ptr]))
    tmp_ptr = run_ptr;
    
    if (tmp_ptr != base_ptr) {
      swap_int(v[tmp_ptr], v[base_ptr]);
      SWAP_int (x[tmp_ptr], x[base_ptr]);
    }
    
    for (run_ptr = base_ptr + 1; (tmp_ptr = run_ptr += 1) <= end_ptr;) {
      
      while (CMP (v[run_ptr], v[tmp_ptr -= 1]))
    ;
      
      if ((tmp_ptr += 1) != run_ptr) {
    int trav;
    
    for (trav = run_ptr + 1; --trav >= run_ptr;) {
      int  c;
      int d;
      c = v[trav];
      d = x[trav];
      int hi, lo;
      
      for (hi = lo = trav; (lo -= 1) >= tmp_ptr; hi = lo) {
        v[hi] = v[lo];
        x[hi] = x[lo];
      }
      v[hi] = c;
      x[hi] = d;
    }
      }
    }
  }

  return 1;
}
