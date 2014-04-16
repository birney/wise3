#ifndef DYNAMITEprobalHEADERFILE
#define DYNAMITEprobalHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna2.h"
#include "dynafunc.h"



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  write_probabilistic_models(dfp,gm,dpi)
 *
 * Descrip:    Makes all the probabilistic routines
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 *
 */
void write_probabilistic_models(DYNFILE * dfp,GenericMatrix * gm,DPImplementation * dpi);


/* Function:  write_forwardbackward_struct(dfp,dpi,gm)
 *
 * Descrip:    Makes the forwardbackward struct
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 */
void write_forwardbackward_struct(DYNFILE * dfp,DPImplementation * dpi,GenericMatrix * gm);


/* Function:  write_backward_mat_sum(dfp,gm,dpi)
 *
 * Descrip:    Makes the backward matrix over all paths
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 *
 */
void write_backward_mat_sum(DYNFILE * dfp,GenericMatrix * gm,DPImplementation * dpi);


/* Function:  write_expl_mat_sum(dfp,gm,dpi)
 *
 * Descrip:    Makes the explicit matrix over all paths
 *             construction method
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 *
 */
void write_expl_mat_sum(DYNFILE * dfp,GenericMatrix * gm,DPImplementation * dpi);


/* Function:  write_score_only_sum(dfp,gm,dpi)
 *
 * Descrip:    Makes the single one-on-one over all paths
 *             searching method
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 *
 */
void write_score_only_sum(DYNFILE * dfp,GenericMatrix * gm,DPImplementation * dpi);


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
