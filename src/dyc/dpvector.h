#ifndef DYNAMITEdpvectorHEADERFILE
#define DYNAMITEdpvectorHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna2.h"
#include "dpimpl.h"




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  write_naive_vector(dfp,*dpi,gm)
 *
 * Descrip:    writes out naive vectorised loop
 *
 *
 * Arg:         dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        *dpi [UNKN ] Undocumented argument [DPImplementation]
 * Arg:          gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean write_naive_vector(DYNFILE * dfp,DPImplementation *dpi,GenericMatrix * gm);


/* Function:  write_naive_vector_calc(dfp,*dpi,gm)
 *
 * Descrip:    writes out an alloc+calc structure for naive vectorisation calc lines
 *
 *
 * Arg:         dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        *dpi [UNKN ] Undocumented argument [DPImplementation]
 * Arg:          gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean write_naive_vector_calc(DYNFILE * dfp,DPImplementation *dpi,GenericMatrix * gm);


/* Function:  write_naive_vector_alloc(dfp,*dpi,gm)
 *
 * Descrip:    writes out an allocation structure for naive vectorisation
 *
 *
 * Arg:         dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        *dpi [UNKN ] Undocumented argument [DPImplementation]
 * Arg:          gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean write_naive_vector_alloc(DYNFILE * dfp,DPImplementation *dpi,GenericMatrix * gm);


/* Function:  write_naive_vector_calc_struct(dfp,dpi,gm)
 *
 * Descrip:    writes out a calc structure for calc layouts
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean write_naive_vector_calc_struct(DYNFILE * dfp,DPImplementation * dpi,GenericMatrix * gm);


/* Function:  write_naive_vector_struct(dfp,*dpi,gm)
 *
 * Descrip:    writes out vector layout structure for naive vectorisation
 *
 *
 * Arg:         dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        *dpi [UNKN ] Undocumented argument [DPImplementation]
 * Arg:          gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean write_naive_vector_struct(DYNFILE * dfp,DPImplementation *dpi,GenericMatrix * gm);


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
