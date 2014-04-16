#ifndef DYNAMITEprevelvetHEADERFILE
#define DYNAMITEprevelvetHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "kmer_.h"

#define MAXKMERTRACE 4
#define MAXKMERREADPAIR 4

typedef struct KmerTrace {
  kmer_t forward;
  kmer_t backward;
  long int count;
  boolean is_dead;
  struct KmerTrace right[MAXKMERTRACE];
  struct KmerTrace left[MAXKMERTRACE];
  struct KmerTrace pair[MAXREADPAIR];
} KmerTrace;


typedef struct KmerTraceSet {
  KmerIndexInterface * kii;
} KmerTraceSet;




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
