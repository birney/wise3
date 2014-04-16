#ifndef DYNAMITEmousethreestateHEADERFILE
#define DYNAMITEmousethreestateHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "snplocus.h"
#include "ancestral.h"

#define TargetMatchMusculus(j) (mat->target->musculus[j] == 1 ? mat->query->match : mat->query->mismatch)
#define TargetMatchCast(j) (mat->target->cast[j] == 1 ? mat->query->match : mat->query->mismatch)
#define TargetMatchDomesticus(j) (mat->target->domesticus[j] == 1 ? mat->query->match : mat->query->mismatch)

struct Wise2_SnpMatch {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    long len;    
    char * musculus;     
    char * cast;     
    char * domesticus;   
    } ;  
/* SnpMatch defined */ 
#ifndef DYNAMITE_DEFINED_SnpMatch
typedef struct Wise2_SnpMatch Wise2_SnpMatch;
#define SnpMatch Wise2_SnpMatch
#define DYNAMITE_DEFINED_SnpMatch
#endif


struct Wise2_GenomePara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score match;     
    Score mismatch;  
    Score extend;    
    Score bswitch;   
    int len;     
    } ;  
/* GenomePara defined */ 
#ifndef DYNAMITE_DEFINED_GenomePara
typedef struct Wise2_GenomePara Wise2_GenomePara;
#define GenomePara Wise2_GenomePara
#define DYNAMITE_DEFINED_GenomePara
#endif


struct Wise2_MouseSNPMatch {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    GenomePara* query;   
    SnpMatch* target;    
    } ;  
/* MouseSNPMatch defined */ 
#ifndef DYNAMITE_DEFINED_MouseSNPMatch
typedef struct Wise2_MouseSNPMatch Wise2_MouseSNPMatch;
#define MouseSNPMatch Wise2_MouseSNPMatch
#define DYNAMITE_DEFINED_MouseSNPMatch
#endif


struct Wise2_MouseSNPMatch_Posterior {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    MouseSNPMatch * forward;     
    MouseSNPMatch * backward;    
    } ;  
/* MouseSNPMatch_Posterior defined */ 
#ifndef DYNAMITE_DEFINED_MouseSNPMatch_Posterior
typedef struct Wise2_MouseSNPMatch_Posterior Wise2_MouseSNPMatch_Posterior;
#define MouseSNPMatch_Posterior Wise2_MouseSNPMatch_Posterior
#define DYNAMITE_DEFINED_MouseSNPMatch_Posterior
#endif


struct Wise2_MouseSNPMatch_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(MouseSNPMatch*,int,int,int);  
    int (*access_special)(MouseSNPMatch*,int,int,int);   
    } ;  
/* MouseSNPMatch_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_MouseSNPMatch_access_func_holder
typedef struct Wise2_MouseSNPMatch_access_func_holder Wise2_MouseSNPMatch_access_func_holder;
#define MouseSNPMatch_access_func_holder Wise2_MouseSNPMatch_access_func_holder
#define DYNAMITE_DEFINED_MouseSNPMatch_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_MouseSNPMatch_Posterior(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MouseSNPMatch_Posterior *]
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch_Posterior *]
 *
 */
MouseSNPMatch_Posterior * Wise2_hard_link_MouseSNPMatch_Posterior(MouseSNPMatch_Posterior * obj);
#define hard_link_MouseSNPMatch_Posterior Wise2_hard_link_MouseSNPMatch_Posterior


/* Function:  MouseSNPMatch_Posterior_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch_Posterior *]
 *
 */
MouseSNPMatch_Posterior * Wise2_MouseSNPMatch_Posterior_alloc(void);
#define MouseSNPMatch_Posterior_alloc Wise2_MouseSNPMatch_Posterior_alloc


/* Function:  free_MouseSNPMatch_Posterior(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MouseSNPMatch_Posterior *]
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch_Posterior *]
 *
 */
MouseSNPMatch_Posterior * Wise2_free_MouseSNPMatch_Posterior(MouseSNPMatch_Posterior * obj);
#define free_MouseSNPMatch_Posterior Wise2_free_MouseSNPMatch_Posterior


/* Function:  hard_link_SnpMatch(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SnpMatch *]
 *
 * Return [UNKN ]  Undocumented return value [SnpMatch *]
 *
 */
SnpMatch * Wise2_hard_link_SnpMatch(SnpMatch * obj);
#define hard_link_SnpMatch Wise2_hard_link_SnpMatch


/* Function:  SnpMatch_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SnpMatch *]
 *
 */
SnpMatch * Wise2_SnpMatch_alloc(void);
#define SnpMatch_alloc Wise2_SnpMatch_alloc


/* Function:  free_SnpMatch(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SnpMatch *]
 *
 * Return [UNKN ]  Undocumented return value [SnpMatch *]
 *
 */
SnpMatch * Wise2_free_SnpMatch(SnpMatch * obj);
#define free_SnpMatch Wise2_free_SnpMatch


/* Function:  hard_link_GenomePara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenomePara *]
 *
 * Return [UNKN ]  Undocumented return value [GenomePara *]
 *
 */
GenomePara * Wise2_hard_link_GenomePara(GenomePara * obj);
#define hard_link_GenomePara Wise2_hard_link_GenomePara


/* Function:  GenomePara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomePara *]
 *
 */
GenomePara * Wise2_GenomePara_alloc(void);
#define GenomePara_alloc Wise2_GenomePara_alloc


/* Function:  free_GenomePara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenomePara *]
 *
 * Return [UNKN ]  Undocumented return value [GenomePara *]
 *
 */
GenomePara * Wise2_free_GenomePara(GenomePara * obj);
#define free_GenomePara Wise2_free_GenomePara


/* Function:  PackAln_read_Shatter_MouseSNPMatch(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_MouseSNPMatch(MouseSNPMatch * mat);
#define PackAln_read_Shatter_MouseSNPMatch Wise2_PackAln_read_Shatter_MouseSNPMatch


/* Function:  calculate_shatter_MouseSNPMatch(mat,dpenv)
 *
 * Descrip:    This function calculates the MouseSNPMatch matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [MouseSNPMatch *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_MouseSNPMatch(MouseSNPMatch * mat,DPEnvelope * dpenv);
#define calculate_shatter_MouseSNPMatch Wise2_calculate_shatter_MouseSNPMatch


/* Function:  search_MouseSNPMatch(dbsi,out,query,target)
 *
 * Descrip:    This function makes a database search of MouseSNPMatch
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:          dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [GenomePara*]
 * Arg:        target [UNKN ] Undocumented argument [SnpMatch*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_MouseSNPMatch(DBSearchImpl * dbsi,Hscore * out,GenomePara* query,SnpMatch* target );
#define search_MouseSNPMatch Wise2_search_MouseSNPMatch


/* Function:  serial_search_MouseSNPMatch(out,query,target)
 *
 * Descrip:    This function makes a database search of MouseSNPMatch
 *             It is a single processor implementation
 *
 *
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [GenomePara*]
 * Arg:        target [UNKN ] Undocumented argument [SnpMatch*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_MouseSNPMatch(Hscore * out,GenomePara* query,SnpMatch* target );
#define serial_search_MouseSNPMatch Wise2_serial_search_MouseSNPMatch


/* Function:  PackAln_bestmemory_MouseSNPMatch(query,target,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_MouseSNPMatch
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePara*]
 * Arg:        target [UNKN ] target data structure [SnpMatch*]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_MouseSNPMatch(GenomePara* query,SnpMatch* target ,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_MouseSNPMatch Wise2_PackAln_bestmemory_MouseSNPMatch


/* Function:  allocate_Expl_MouseSNPMatch(query,target,dpri)
 *
 * Descrip:    This function allocates the MouseSNPMatch structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_MouseSNPMatch_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePara*]
 * Arg:        target [UNKN ] target data structure [SnpMatch*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch *]
 *
 */
MouseSNPMatch * Wise2_allocate_Expl_MouseSNPMatch(GenomePara* query,SnpMatch* target ,DPRunImpl * dpri);
#define allocate_Expl_MouseSNPMatch Wise2_allocate_Expl_MouseSNPMatch


/* Function:  recalculate_PackAln_MouseSNPMatch(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by MouseSNPMatch
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 *
 */
void Wise2_recalculate_PackAln_MouseSNPMatch(PackAln * pal,MouseSNPMatch * mat);
#define recalculate_PackAln_MouseSNPMatch Wise2_recalculate_PackAln_MouseSNPMatch


/* Function:  allocate_Small_MouseSNPMatch(query,target)
 *
 * Descrip:    This function allocates the MouseSNPMatch structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_MouseSNPMatch_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePara*]
 * Arg:        target [UNKN ] target data structure [SnpMatch*]
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch *]
 *
 */
MouseSNPMatch * Wise2_allocate_Small_MouseSNPMatch(GenomePara* query,SnpMatch* target );
#define allocate_Small_MouseSNPMatch Wise2_allocate_Small_MouseSNPMatch


/* Function:  PackAln_calculate_Small_MouseSNPMatch(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for MouseSNPMatch structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_MouseSNPMatch 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_MouseSNPMatch 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_MouseSNPMatch(MouseSNPMatch * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_MouseSNPMatch Wise2_PackAln_calculate_Small_MouseSNPMatch


/* Function:  AlnRangeSet_calculate_Small_MouseSNPMatch(mat)
 *
 * Descrip:    This function calculates an alignment for MouseSNPMatch structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_MouseSNPMatch 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_MouseSNPMatch
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_MouseSNPMatch 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_MouseSNPMatch(MouseSNPMatch * mat);
#define AlnRangeSet_calculate_Small_MouseSNPMatch Wise2_AlnRangeSet_calculate_Small_MouseSNPMatch


/* Function:  AlnRangeSet_from_MouseSNPMatch(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for MouseSNPMatch structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_MouseSNPMatch 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_MouseSNPMatch
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_MouseSNPMatch(MouseSNPMatch * mat);
#define AlnRangeSet_from_MouseSNPMatch Wise2_AlnRangeSet_from_MouseSNPMatch


/* Function:  convert_PackAln_to_AlnBlock_MouseSNPMatch(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_MouseSNPMatch(PackAln * pal);
#define convert_PackAln_to_AlnBlock_MouseSNPMatch Wise2_convert_PackAln_to_AlnBlock_MouseSNPMatch


/* Function:  PackAln_read_Expl_MouseSNPMatch(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_MouseSNPMatch(MouseSNPMatch * mat);
#define PackAln_read_Expl_MouseSNPMatch Wise2_PackAln_read_Expl_MouseSNPMatch


/* Function:  PackAln_read_generic_MouseSNPMatch(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:          h [UNKN ] Undocumented argument [MouseSNPMatch_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_MouseSNPMatch(MouseSNPMatch * mat,MouseSNPMatch_access_func_holder h);
#define PackAln_read_generic_MouseSNPMatch Wise2_PackAln_read_generic_MouseSNPMatch


/* Function:  calculate_MouseSNPMatch(mat)
 *
 * Descrip:    This function calculates the MouseSNPMatch matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_MouseSNPMatch
 *
 *
 * Arg:        mat [UNKN ] MouseSNPMatch which contains explicit basematrix memory [MouseSNPMatch *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_MouseSNPMatch(MouseSNPMatch * mat);
#define calculate_MouseSNPMatch Wise2_calculate_MouseSNPMatch


/* Function:  calculate_dpenv_MouseSNPMatch(mat,dpenv)
 *
 * Descrip:    This function calculates the MouseSNPMatch matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] MouseSNPMatch which contains explicit basematrix memory [MouseSNPMatch *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_MouseSNPMatch(MouseSNPMatch * mat,DPEnvelope * dpenv);
#define calculate_dpenv_MouseSNPMatch Wise2_calculate_dpenv_MouseSNPMatch


/* Function:  MouseSNPMatch_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch *]
 *
 */
MouseSNPMatch * Wise2_MouseSNPMatch_alloc(void);
#define MouseSNPMatch_alloc Wise2_MouseSNPMatch_alloc


/* Function:  free_MouseSNPMatch(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MouseSNPMatch *]
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch *]
 *
 */
MouseSNPMatch * Wise2_free_MouseSNPMatch(MouseSNPMatch * obj);
#define free_MouseSNPMatch Wise2_free_MouseSNPMatch


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
AncestralVarSet * Wise2_ancestral_prediction_MouseThreeState(GenoVarSet * gvs,char * mus,char * cast,char * dom,GenomePara * gp,DPRunImpl * dpri);
#define ancestral_prediction_MouseThreeState Wise2_ancestral_prediction_MouseThreeState
void Wise2_populate_individual_Ancestral_decoding_ThreeState(AlnBlock * alb,AncestralVarChr * avc,int ind_number,short int mus_index,short int cast_index,short int dom_index);
#define populate_individual_Ancestral_decoding_ThreeState Wise2_populate_individual_Ancestral_decoding_ThreeState
void Wise2_create_blank_AncestralChromosomes(AncestralVarSet * avs,GenoVarSet * gvs);
#define create_blank_AncestralChromosomes Wise2_create_blank_AncestralChromosomes
GenomePara * Wise2_new_GenomePara(Probability match,Probability bswitch);
#define new_GenomePara Wise2_new_GenomePara
void Wise2_show_SnpMatchStats(SnpMatch * snpm,FILE * ofp);
#define show_SnpMatchStats Wise2_show_SnpMatchStats
SnpMatch * Wise2_new_SnpMatch(GenoVarChr * chr,GenoVarSet * set,char * test_strain,char * musculus,char * cast,char * domesticus);
#define new_SnpMatch Wise2_new_SnpMatch


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_MouseSNPMatch_shatter_access_main(MouseSNPMatch * mat,int i,int j,int state);
#define MouseSNPMatch_shatter_access_main Wise2_MouseSNPMatch_shatter_access_main
int Wise2_MouseSNPMatch_shatter_access_special(MouseSNPMatch * mat,int i,int j,int state);
#define MouseSNPMatch_shatter_access_special Wise2_MouseSNPMatch_shatter_access_special
Score Wise2_score_only_logsum_MouseSNPMatch(GenomePara* query,SnpMatch* target );
#define score_only_logsum_MouseSNPMatch Wise2_score_only_logsum_MouseSNPMatch
MouseSNPMatch * Wise2_forward_logsum_MouseSNPMatch(GenomePara* query,SnpMatch* target ,DPRunImpl * dpri);
#define forward_logsum_MouseSNPMatch Wise2_forward_logsum_MouseSNPMatch
MouseSNPMatch * Wise2_backward_logsum_MouseSNPMatch(GenomePara* query,SnpMatch* target ,DPRunImpl * dpri);
#define backward_logsum_MouseSNPMatch Wise2_backward_logsum_MouseSNPMatch
int Wise2_score_only_MouseSNPMatch(GenomePara* query,SnpMatch* target );
#define score_only_MouseSNPMatch Wise2_score_only_MouseSNPMatch
MouseSNPMatch * Wise2_allocate_MouseSNPMatch_only(GenomePara* query,SnpMatch* target );
#define allocate_MouseSNPMatch_only Wise2_allocate_MouseSNPMatch_only
void Wise2_init_MouseSNPMatch(MouseSNPMatch * mat);
#define init_MouseSNPMatch Wise2_init_MouseSNPMatch
AlnRange * Wise2_AlnRange_build_MouseSNPMatch(MouseSNPMatch * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_MouseSNPMatch Wise2_AlnRange_build_MouseSNPMatch
boolean Wise2_read_hidden_MouseSNPMatch(MouseSNPMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_MouseSNPMatch Wise2_read_hidden_MouseSNPMatch
int Wise2_max_hidden_MouseSNPMatch(MouseSNPMatch * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_MouseSNPMatch Wise2_max_hidden_MouseSNPMatch
boolean Wise2_read_special_strip_MouseSNPMatch(MouseSNPMatch * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_MouseSNPMatch Wise2_read_special_strip_MouseSNPMatch
int Wise2_max_special_strip_MouseSNPMatch(MouseSNPMatch * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_MouseSNPMatch Wise2_max_special_strip_MouseSNPMatch
int Wise2_max_matrix_to_special_MouseSNPMatch(MouseSNPMatch * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_MouseSNPMatch Wise2_max_matrix_to_special_MouseSNPMatch
void Wise2_calculate_hidden_MouseSNPMatch(MouseSNPMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_MouseSNPMatch Wise2_calculate_hidden_MouseSNPMatch
void Wise2_init_hidden_MouseSNPMatch(MouseSNPMatch * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_MouseSNPMatch Wise2_init_hidden_MouseSNPMatch
boolean Wise2_full_dc_MouseSNPMatch(MouseSNPMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_MouseSNPMatch Wise2_full_dc_MouseSNPMatch
boolean Wise2_do_dc_single_pass_MouseSNPMatch(MouseSNPMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_MouseSNPMatch Wise2_do_dc_single_pass_MouseSNPMatch
void Wise2_push_dc_at_merge_MouseSNPMatch(MouseSNPMatch * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_MouseSNPMatch Wise2_push_dc_at_merge_MouseSNPMatch
void Wise2_follow_on_dc_MouseSNPMatch(MouseSNPMatch * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_MouseSNPMatch Wise2_follow_on_dc_MouseSNPMatch
void Wise2_run_up_dc_MouseSNPMatch(MouseSNPMatch * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_MouseSNPMatch Wise2_run_up_dc_MouseSNPMatch
void Wise2_init_dc_MouseSNPMatch(MouseSNPMatch * mat);
#define init_dc_MouseSNPMatch Wise2_init_dc_MouseSNPMatch
int Wise2_start_end_find_end_MouseSNPMatch(MouseSNPMatch * mat,int * endj);
#define start_end_find_end_MouseSNPMatch Wise2_start_end_find_end_MouseSNPMatch
boolean Wise2_dc_optimised_start_end_calc_MouseSNPMatch(MouseSNPMatch *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_MouseSNPMatch Wise2_dc_optimised_start_end_calc_MouseSNPMatch
void Wise2_init_start_end_linear_MouseSNPMatch(MouseSNPMatch * mat);
#define init_start_end_linear_MouseSNPMatch Wise2_init_start_end_linear_MouseSNPMatch
AlnConvertSet * Wise2_AlnConvertSet_MouseSNPMatch(void);
#define AlnConvertSet_MouseSNPMatch Wise2_AlnConvertSet_MouseSNPMatch
int Wise2_MouseSNPMatch_explicit_access_main(MouseSNPMatch * mat,int i,int j,int state);
#define MouseSNPMatch_explicit_access_main Wise2_MouseSNPMatch_explicit_access_main
int Wise2_MouseSNPMatch_explicit_access_special(MouseSNPMatch * mat,int i,int j,int state);
#define MouseSNPMatch_explicit_access_special Wise2_MouseSNPMatch_explicit_access_special
int Wise2_find_end_MouseSNPMatch(MouseSNPMatch * mat,int * ri,int * rj,int * state,boolean * isspecial,MouseSNPMatch_access_func_holder h);
#define find_end_MouseSNPMatch Wise2_find_end_MouseSNPMatch
void Wise2_MouseSNPMatch_debug_show_matrix(MouseSNPMatch * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define MouseSNPMatch_debug_show_matrix Wise2_MouseSNPMatch_debug_show_matrix
int Wise2_max_calc_MouseSNPMatch(MouseSNPMatch * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,MouseSNPMatch_access_func_holder h);
#define max_calc_MouseSNPMatch Wise2_max_calc_MouseSNPMatch
int Wise2_max_calc_special_MouseSNPMatch(MouseSNPMatch * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,MouseSNPMatch_access_func_holder h);
#define max_calc_special_MouseSNPMatch Wise2_max_calc_special_MouseSNPMatch

#ifdef _cplusplus
}
#endif

#endif
