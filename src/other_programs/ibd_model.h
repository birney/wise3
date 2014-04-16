#ifndef DYNAMITEibd_modelHEADERFILE
#define DYNAMITEibd_modelHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "snplocus.h"
#include "ancestral.h"

#include "mousethreestate.h"

#define TargetMatchIBD(j) (mat->target->match[j] == 1 ? mat->query->ibd_match : mat->query->ibd_mismatch)

#define TargetMatchDIFF(j) (mat->target->match[j] == 1 ? mat->query->diff_match : mat->query->diff_mismatch)

struct Wise2_PairMatch {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    long len;    
    char * match;    
    } ;  
/* PairMatch defined */ 
#ifndef DYNAMITE_DEFINED_PairMatch
typedef struct Wise2_PairMatch Wise2_PairMatch;
#define PairMatch Wise2_PairMatch
#define DYNAMITE_DEFINED_PairMatch
#endif


struct Wise2_GenomePairPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score ibd_match;     
    Score ibd_mismatch;  
    Score diff_match;    
    Score diff_mismatch;     
    Score extend_ibd;    
    Score bswitch_ibd;   
    Score extend_diff;   
    Score bswitch_diff;  
    int len;     
    } ;  
/* GenomePairPara defined */ 
#ifndef DYNAMITE_DEFINED_GenomePairPara
typedef struct Wise2_GenomePairPara Wise2_GenomePairPara;
#define GenomePairPara Wise2_GenomePairPara
#define DYNAMITE_DEFINED_GenomePairPara
#endif


struct Wise2_MousePairMatch {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    GenomePairPara* query;   
    PairMatch* target;   
    } ;  
/* MousePairMatch defined */ 
#ifndef DYNAMITE_DEFINED_MousePairMatch
typedef struct Wise2_MousePairMatch Wise2_MousePairMatch;
#define MousePairMatch Wise2_MousePairMatch
#define DYNAMITE_DEFINED_MousePairMatch
#endif


struct Wise2_MousePairMatch_Posterior {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    MousePairMatch * forward;    
    MousePairMatch * backward;   
    } ;  
/* MousePairMatch_Posterior defined */ 
#ifndef DYNAMITE_DEFINED_MousePairMatch_Posterior
typedef struct Wise2_MousePairMatch_Posterior Wise2_MousePairMatch_Posterior;
#define MousePairMatch_Posterior Wise2_MousePairMatch_Posterior
#define DYNAMITE_DEFINED_MousePairMatch_Posterior
#endif


struct Wise2_MousePairMatch_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(MousePairMatch*,int,int,int); 
    int (*access_special)(MousePairMatch*,int,int,int);  
    } ;  
/* MousePairMatch_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_MousePairMatch_access_func_holder
typedef struct Wise2_MousePairMatch_access_func_holder Wise2_MousePairMatch_access_func_holder;
#define MousePairMatch_access_func_holder Wise2_MousePairMatch_access_func_holder
#define DYNAMITE_DEFINED_MousePairMatch_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_MousePairMatch_Posterior(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MousePairMatch_Posterior *]
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch_Posterior *]
 *
 */
MousePairMatch_Posterior * Wise2_hard_link_MousePairMatch_Posterior(MousePairMatch_Posterior * obj);
#define hard_link_MousePairMatch_Posterior Wise2_hard_link_MousePairMatch_Posterior


/* Function:  MousePairMatch_Posterior_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch_Posterior *]
 *
 */
MousePairMatch_Posterior * Wise2_MousePairMatch_Posterior_alloc(void);
#define MousePairMatch_Posterior_alloc Wise2_MousePairMatch_Posterior_alloc


/* Function:  free_MousePairMatch_Posterior(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MousePairMatch_Posterior *]
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch_Posterior *]
 *
 */
MousePairMatch_Posterior * Wise2_free_MousePairMatch_Posterior(MousePairMatch_Posterior * obj);
#define free_MousePairMatch_Posterior Wise2_free_MousePairMatch_Posterior


/* Function:  hard_link_PairMatch(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [PairMatch *]
 *
 */
PairMatch * Wise2_hard_link_PairMatch(PairMatch * obj);
#define hard_link_PairMatch Wise2_hard_link_PairMatch


/* Function:  PairMatch_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairMatch *]
 *
 */
PairMatch * Wise2_PairMatch_alloc(void);
#define PairMatch_alloc Wise2_PairMatch_alloc


/* Function:  free_PairMatch(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [PairMatch *]
 *
 */
PairMatch * Wise2_free_PairMatch(PairMatch * obj);
#define free_PairMatch Wise2_free_PairMatch


/* Function:  hard_link_GenomePairPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenomePairPara *]
 *
 * Return [UNKN ]  Undocumented return value [GenomePairPara *]
 *
 */
GenomePairPara * Wise2_hard_link_GenomePairPara(GenomePairPara * obj);
#define hard_link_GenomePairPara Wise2_hard_link_GenomePairPara


/* Function:  GenomePairPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomePairPara *]
 *
 */
GenomePairPara * Wise2_GenomePairPara_alloc(void);
#define GenomePairPara_alloc Wise2_GenomePairPara_alloc


/* Function:  free_GenomePairPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenomePairPara *]
 *
 * Return [UNKN ]  Undocumented return value [GenomePairPara *]
 *
 */
GenomePairPara * Wise2_free_GenomePairPara(GenomePairPara * obj);
#define free_GenomePairPara Wise2_free_GenomePairPara


/* Function:  PackAln_read_Shatter_MousePairMatch(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MousePairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_MousePairMatch(MousePairMatch * mat);
#define PackAln_read_Shatter_MousePairMatch Wise2_PackAln_read_Shatter_MousePairMatch


/* Function:  calculate_shatter_MousePairMatch(mat,dpenv)
 *
 * Descrip:    This function calculates the MousePairMatch matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [MousePairMatch *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_MousePairMatch(MousePairMatch * mat,DPEnvelope * dpenv);
#define calculate_shatter_MousePairMatch Wise2_calculate_shatter_MousePairMatch


/* Function:  search_MousePairMatch(dbsi,out,query,target)
 *
 * Descrip:    This function makes a database search of MousePairMatch
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:          dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [GenomePairPara*]
 * Arg:        target [UNKN ] Undocumented argument [PairMatch*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_MousePairMatch(DBSearchImpl * dbsi,Hscore * out,GenomePairPara* query,PairMatch* target );
#define search_MousePairMatch Wise2_search_MousePairMatch


/* Function:  serial_search_MousePairMatch(out,query,target)
 *
 * Descrip:    This function makes a database search of MousePairMatch
 *             It is a single processor implementation
 *
 *
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [GenomePairPara*]
 * Arg:        target [UNKN ] Undocumented argument [PairMatch*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_MousePairMatch(Hscore * out,GenomePairPara* query,PairMatch* target );
#define serial_search_MousePairMatch Wise2_serial_search_MousePairMatch


/* Function:  PackAln_bestmemory_MousePairMatch(query,target,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_MousePairMatch
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePairPara*]
 * Arg:        target [UNKN ] target data structure [PairMatch*]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_MousePairMatch(GenomePairPara* query,PairMatch* target ,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_MousePairMatch Wise2_PackAln_bestmemory_MousePairMatch


/* Function:  allocate_Expl_MousePairMatch(query,target,dpri)
 *
 * Descrip:    This function allocates the MousePairMatch structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_MousePairMatch_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePairPara*]
 * Arg:        target [UNKN ] target data structure [PairMatch*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch *]
 *
 */
MousePairMatch * Wise2_allocate_Expl_MousePairMatch(GenomePairPara* query,PairMatch* target ,DPRunImpl * dpri);
#define allocate_Expl_MousePairMatch Wise2_allocate_Expl_MousePairMatch


/* Function:  recalculate_PackAln_MousePairMatch(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by MousePairMatch
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [MousePairMatch *]
 *
 */
void Wise2_recalculate_PackAln_MousePairMatch(PackAln * pal,MousePairMatch * mat);
#define recalculate_PackAln_MousePairMatch Wise2_recalculate_PackAln_MousePairMatch


/* Function:  allocate_Small_MousePairMatch(query,target)
 *
 * Descrip:    This function allocates the MousePairMatch structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_MousePairMatch_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePairPara*]
 * Arg:        target [UNKN ] target data structure [PairMatch*]
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch *]
 *
 */
MousePairMatch * Wise2_allocate_Small_MousePairMatch(GenomePairPara* query,PairMatch* target );
#define allocate_Small_MousePairMatch Wise2_allocate_Small_MousePairMatch


/* Function:  PackAln_calculate_Small_MousePairMatch(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for MousePairMatch structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_MousePairMatch 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_MousePairMatch 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_MousePairMatch(MousePairMatch * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_MousePairMatch Wise2_PackAln_calculate_Small_MousePairMatch


/* Function:  AlnRangeSet_calculate_Small_MousePairMatch(mat)
 *
 * Descrip:    This function calculates an alignment for MousePairMatch structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_MousePairMatch 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_MousePairMatch
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_MousePairMatch 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MousePairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_MousePairMatch(MousePairMatch * mat);
#define AlnRangeSet_calculate_Small_MousePairMatch Wise2_AlnRangeSet_calculate_Small_MousePairMatch


/* Function:  AlnRangeSet_from_MousePairMatch(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for MousePairMatch structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_MousePairMatch 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_MousePairMatch
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MousePairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_MousePairMatch(MousePairMatch * mat);
#define AlnRangeSet_from_MousePairMatch Wise2_AlnRangeSet_from_MousePairMatch


/* Function:  convert_PackAln_to_AlnBlock_MousePairMatch(pal)
 *
 * Descrip:    Converts a path alignment to a label alignment
 *             The label alignment is probably much more useful than the path
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_MousePairMatch(PackAln * pal);
#define convert_PackAln_to_AlnBlock_MousePairMatch Wise2_convert_PackAln_to_AlnBlock_MousePairMatch


/* Function:  PackAln_read_Expl_MousePairMatch(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MousePairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_MousePairMatch(MousePairMatch * mat);
#define PackAln_read_Expl_MousePairMatch Wise2_PackAln_read_Expl_MousePairMatch


/* Function:  PackAln_read_generic_MousePairMatch(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:          h [UNKN ] Undocumented argument [MousePairMatch_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_MousePairMatch(MousePairMatch * mat,MousePairMatch_access_func_holder h);
#define PackAln_read_generic_MousePairMatch Wise2_PackAln_read_generic_MousePairMatch


/* Function:  calculate_MousePairMatch(mat)
 *
 * Descrip:    This function calculates the MousePairMatch matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_MousePairMatch
 *
 *
 * Arg:        mat [UNKN ] MousePairMatch which contains explicit basematrix memory [MousePairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_MousePairMatch(MousePairMatch * mat);
#define calculate_MousePairMatch Wise2_calculate_MousePairMatch


/* Function:  calculate_dpenv_MousePairMatch(mat,dpenv)
 *
 * Descrip:    This function calculates the MousePairMatch matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] MousePairMatch which contains explicit basematrix memory [MousePairMatch *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_MousePairMatch(MousePairMatch * mat,DPEnvelope * dpenv);
#define calculate_dpenv_MousePairMatch Wise2_calculate_dpenv_MousePairMatch


/* Function:  MousePairMatch_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch *]
 *
 */
MousePairMatch * Wise2_MousePairMatch_alloc(void);
#define MousePairMatch_alloc Wise2_MousePairMatch_alloc


/* Function:  free_MousePairMatch(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MousePairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch *]
 *
 */
MousePairMatch * Wise2_free_MousePairMatch(MousePairMatch * obj);
#define free_MousePairMatch Wise2_free_MousePairMatch


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_dump_naive_AlnBlock_IBD(char * strain_a,char * strain_b,GenoVarChr * chr,AlnBlock * alb,FILE * ofp);
#define dump_naive_AlnBlock_IBD Wise2_dump_naive_AlnBlock_IBD
PairMatch * Wise2_new_PairMatch(GenoVarSet * set,GenoVarChr * chr,char * strain_a,char * strain_b);
#define new_PairMatch Wise2_new_PairMatch
GenomePairPara * Wise2_new_GenomePairPara_null(Probability ibd_mismatch,Probability diff_mismatch,Probability ibd_switch,Probability diff_switch);
#define new_GenomePairPara_null Wise2_new_GenomePairPara_null
GenomePairPara * Wise2_new_GenomePairPara(Probability ibd_mismatch,Probability diff_mismatch,Probability ibd_switch,Probability diff_switch);
#define new_GenomePairPara Wise2_new_GenomePairPara


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_MousePairMatch_shatter_access_main(MousePairMatch * mat,int i,int j,int state);
#define MousePairMatch_shatter_access_main Wise2_MousePairMatch_shatter_access_main
int Wise2_MousePairMatch_shatter_access_special(MousePairMatch * mat,int i,int j,int state);
#define MousePairMatch_shatter_access_special Wise2_MousePairMatch_shatter_access_special
Score Wise2_score_only_logsum_MousePairMatch(GenomePairPara* query,PairMatch* target );
#define score_only_logsum_MousePairMatch Wise2_score_only_logsum_MousePairMatch
MousePairMatch * Wise2_forward_logsum_MousePairMatch(GenomePairPara* query,PairMatch* target ,DPRunImpl * dpri);
#define forward_logsum_MousePairMatch Wise2_forward_logsum_MousePairMatch
MousePairMatch * Wise2_backward_logsum_MousePairMatch(GenomePairPara* query,PairMatch* target ,DPRunImpl * dpri);
#define backward_logsum_MousePairMatch Wise2_backward_logsum_MousePairMatch
int Wise2_score_only_MousePairMatch(GenomePairPara* query,PairMatch* target );
#define score_only_MousePairMatch Wise2_score_only_MousePairMatch
MousePairMatch * Wise2_allocate_MousePairMatch_only(GenomePairPara* query,PairMatch* target );
#define allocate_MousePairMatch_only Wise2_allocate_MousePairMatch_only
void Wise2_init_MousePairMatch(MousePairMatch * mat);
#define init_MousePairMatch Wise2_init_MousePairMatch
AlnRange * Wise2_AlnRange_build_MousePairMatch(MousePairMatch * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_MousePairMatch Wise2_AlnRange_build_MousePairMatch
boolean Wise2_read_hidden_MousePairMatch(MousePairMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_MousePairMatch Wise2_read_hidden_MousePairMatch
int Wise2_max_hidden_MousePairMatch(MousePairMatch * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_MousePairMatch Wise2_max_hidden_MousePairMatch
boolean Wise2_read_special_strip_MousePairMatch(MousePairMatch * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_MousePairMatch Wise2_read_special_strip_MousePairMatch
int Wise2_max_special_strip_MousePairMatch(MousePairMatch * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_MousePairMatch Wise2_max_special_strip_MousePairMatch
int Wise2_max_matrix_to_special_MousePairMatch(MousePairMatch * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_MousePairMatch Wise2_max_matrix_to_special_MousePairMatch
void Wise2_calculate_hidden_MousePairMatch(MousePairMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_MousePairMatch Wise2_calculate_hidden_MousePairMatch
void Wise2_init_hidden_MousePairMatch(MousePairMatch * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_MousePairMatch Wise2_init_hidden_MousePairMatch
boolean Wise2_full_dc_MousePairMatch(MousePairMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_MousePairMatch Wise2_full_dc_MousePairMatch
boolean Wise2_do_dc_single_pass_MousePairMatch(MousePairMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_MousePairMatch Wise2_do_dc_single_pass_MousePairMatch
void Wise2_push_dc_at_merge_MousePairMatch(MousePairMatch * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_MousePairMatch Wise2_push_dc_at_merge_MousePairMatch
void Wise2_follow_on_dc_MousePairMatch(MousePairMatch * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_MousePairMatch Wise2_follow_on_dc_MousePairMatch
void Wise2_run_up_dc_MousePairMatch(MousePairMatch * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_MousePairMatch Wise2_run_up_dc_MousePairMatch
void Wise2_init_dc_MousePairMatch(MousePairMatch * mat);
#define init_dc_MousePairMatch Wise2_init_dc_MousePairMatch
int Wise2_start_end_find_end_MousePairMatch(MousePairMatch * mat,int * endj);
#define start_end_find_end_MousePairMatch Wise2_start_end_find_end_MousePairMatch
boolean Wise2_dc_optimised_start_end_calc_MousePairMatch(MousePairMatch *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_MousePairMatch Wise2_dc_optimised_start_end_calc_MousePairMatch
void Wise2_init_start_end_linear_MousePairMatch(MousePairMatch * mat);
#define init_start_end_linear_MousePairMatch Wise2_init_start_end_linear_MousePairMatch
AlnConvertSet * Wise2_AlnConvertSet_MousePairMatch(void);
#define AlnConvertSet_MousePairMatch Wise2_AlnConvertSet_MousePairMatch
int Wise2_MousePairMatch_explicit_access_main(MousePairMatch * mat,int i,int j,int state);
#define MousePairMatch_explicit_access_main Wise2_MousePairMatch_explicit_access_main
int Wise2_MousePairMatch_explicit_access_special(MousePairMatch * mat,int i,int j,int state);
#define MousePairMatch_explicit_access_special Wise2_MousePairMatch_explicit_access_special
int Wise2_find_end_MousePairMatch(MousePairMatch * mat,int * ri,int * rj,int * state,boolean * isspecial,MousePairMatch_access_func_holder h);
#define find_end_MousePairMatch Wise2_find_end_MousePairMatch
void Wise2_MousePairMatch_debug_show_matrix(MousePairMatch * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define MousePairMatch_debug_show_matrix Wise2_MousePairMatch_debug_show_matrix
int Wise2_max_calc_MousePairMatch(MousePairMatch * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,MousePairMatch_access_func_holder h);
#define max_calc_MousePairMatch Wise2_max_calc_MousePairMatch
int Wise2_max_calc_special_MousePairMatch(MousePairMatch * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,MousePairMatch_access_func_holder h);
#define max_calc_special_MousePairMatch Wise2_max_calc_special_MousePairMatch

#ifdef _cplusplus
}
#endif

#endif
