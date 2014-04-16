#ifndef DYNAMITEtwo_track_refHEADERFILE
#define DYNAMITEtwo_track_refHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "two_track.h"


struct Wise2_TwoTrackRefSingleState {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    TwoTrack* query;     
    ComplexSequence* target;     
    Score insert;    
    Score delete;    
    } ;  
/* TwoTrackRefSingleState defined */ 
#ifndef DYNAMITE_DEFINED_TwoTrackRefSingleState
typedef struct Wise2_TwoTrackRefSingleState Wise2_TwoTrackRefSingleState;
#define TwoTrackRefSingleState Wise2_TwoTrackRefSingleState
#define DYNAMITE_DEFINED_TwoTrackRefSingleState
#endif


struct Wise2_TwoTrackRefSingleState_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(TwoTrackRefSingleState*,int,int,int); 
    int (*access_special)(TwoTrackRefSingleState*,int,int,int);  
    } ;  
/* TwoTrackRefSingleState_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_TwoTrackRefSingleState_access_func_holder
typedef struct Wise2_TwoTrackRefSingleState_access_func_holder Wise2_TwoTrackRefSingleState_access_func_holder;
#define TwoTrackRefSingleState_access_func_holder Wise2_TwoTrackRefSingleState_access_func_holder
#define DYNAMITE_DEFINED_TwoTrackRefSingleState_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_TwoTrackRefSingleState(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [TwoTrackRefSingleState *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_TwoTrackRefSingleState(TwoTrackRefSingleState * mat);
#define PackAln_read_Shatter_TwoTrackRefSingleState Wise2_PackAln_read_Shatter_TwoTrackRefSingleState


/* Function:  calculate_shatter_TwoTrackRefSingleState(mat,dpenv)
 *
 * Descrip:    This function calculates the TwoTrackRefSingleState matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [TwoTrackRefSingleState *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,DPEnvelope * dpenv);
#define calculate_shatter_TwoTrackRefSingleState Wise2_calculate_shatter_TwoTrackRefSingleState


/* Function:  search_TwoTrackRefSingleState(dbsi,out,query,target,insert,delete)
 *
 * Descrip:    This function makes a database search of TwoTrackRefSingleState
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:          dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [TwoTrack*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        insert [UNKN ] Undocumented argument [Score]
 * Arg:        delete [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_TwoTrackRefSingleState(DBSearchImpl * dbsi,Hscore * out,TwoTrack* query,ComplexSequence* target ,Score insert,Score delete);
#define search_TwoTrackRefSingleState Wise2_search_TwoTrackRefSingleState


/* Function:  serial_search_TwoTrackRefSingleState(out,query,target,insert,delete)
 *
 * Descrip:    This function makes a database search of TwoTrackRefSingleState
 *             It is a single processor implementation
 *
 *
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [TwoTrack*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        insert [UNKN ] Undocumented argument [Score]
 * Arg:        delete [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_TwoTrackRefSingleState(Hscore * out,TwoTrack* query,ComplexSequence* target ,Score insert,Score delete);
#define serial_search_TwoTrackRefSingleState Wise2_serial_search_TwoTrackRefSingleState


/* Function:  PackAln_bestmemory_TwoTrackRefSingleState(query,target,insert,delete,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_TwoTrackRefSingleState
 *
 *
 * Arg:         query [UNKN ] query data structure [TwoTrack*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:        insert [UNKN ] Resource [Score]
 * Arg:        delete [UNKN ] Resource [Score]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_TwoTrackRefSingleState(TwoTrack* query,ComplexSequence* target ,Score insert,Score delete,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_TwoTrackRefSingleState Wise2_PackAln_bestmemory_TwoTrackRefSingleState


/* Function:  allocate_Expl_TwoTrackRefSingleState(query,target,insert,delete,dpri)
 *
 * Descrip:    This function allocates the TwoTrackRefSingleState structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_TwoTrackRefSingleState_only
 *
 *
 * Arg:         query [UNKN ] query data structure [TwoTrack*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:        insert [UNKN ] Resource [Score]
 * Arg:        delete [UNKN ] Resource [Score]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackRefSingleState *]
 *
 */
TwoTrackRefSingleState * Wise2_allocate_Expl_TwoTrackRefSingleState(TwoTrack* query,ComplexSequence* target ,Score insert,Score delete,DPRunImpl * dpri);
#define allocate_Expl_TwoTrackRefSingleState Wise2_allocate_Expl_TwoTrackRefSingleState


/* Function:  recalculate_PackAln_TwoTrackRefSingleState(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by TwoTrackRefSingleState
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [TwoTrackRefSingleState *]
 *
 */
void Wise2_recalculate_PackAln_TwoTrackRefSingleState(PackAln * pal,TwoTrackRefSingleState * mat);
#define recalculate_PackAln_TwoTrackRefSingleState Wise2_recalculate_PackAln_TwoTrackRefSingleState


/* Function:  allocate_Small_TwoTrackRefSingleState(query,target,insert,delete)
 *
 * Descrip:    This function allocates the TwoTrackRefSingleState structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_TwoTrackRefSingleState_only
 *
 *
 * Arg:         query [UNKN ] query data structure [TwoTrack*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:        insert [UNKN ] Resource [Score]
 * Arg:        delete [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackRefSingleState *]
 *
 */
TwoTrackRefSingleState * Wise2_allocate_Small_TwoTrackRefSingleState(TwoTrack* query,ComplexSequence* target ,Score insert,Score delete);
#define allocate_Small_TwoTrackRefSingleState Wise2_allocate_Small_TwoTrackRefSingleState


/* Function:  PackAln_calculate_Small_TwoTrackRefSingleState(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for TwoTrackRefSingleState structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_TwoTrackRefSingleState 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_TwoTrackRefSingleState 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [TwoTrackRefSingleState *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_TwoTrackRefSingleState Wise2_PackAln_calculate_Small_TwoTrackRefSingleState


/* Function:  AlnRangeSet_calculate_Small_TwoTrackRefSingleState(mat)
 *
 * Descrip:    This function calculates an alignment for TwoTrackRefSingleState structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_TwoTrackRefSingleState 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_TwoTrackRefSingleState
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_TwoTrackRefSingleState 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [TwoTrackRefSingleState *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_TwoTrackRefSingleState(TwoTrackRefSingleState * mat);
#define AlnRangeSet_calculate_Small_TwoTrackRefSingleState Wise2_AlnRangeSet_calculate_Small_TwoTrackRefSingleState


/* Function:  AlnRangeSet_from_TwoTrackRefSingleState(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for TwoTrackRefSingleState structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_TwoTrackRefSingleState 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_TwoTrackRefSingleState
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [TwoTrackRefSingleState *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_TwoTrackRefSingleState(TwoTrackRefSingleState * mat);
#define AlnRangeSet_from_TwoTrackRefSingleState Wise2_AlnRangeSet_from_TwoTrackRefSingleState


/* Function:  convert_PackAln_to_AlnBlock_TwoTrackRefSingleState(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_TwoTrackRefSingleState(PackAln * pal);
#define convert_PackAln_to_AlnBlock_TwoTrackRefSingleState Wise2_convert_PackAln_to_AlnBlock_TwoTrackRefSingleState


/* Function:  PackAln_read_Expl_TwoTrackRefSingleState(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [TwoTrackRefSingleState *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_TwoTrackRefSingleState(TwoTrackRefSingleState * mat);
#define PackAln_read_Expl_TwoTrackRefSingleState Wise2_PackAln_read_Expl_TwoTrackRefSingleState


/* Function:  PackAln_read_generic_TwoTrackRefSingleState(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [TwoTrackRefSingleState *]
 * Arg:          h [UNKN ] Undocumented argument [TwoTrackRefSingleState_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,TwoTrackRefSingleState_access_func_holder h);
#define PackAln_read_generic_TwoTrackRefSingleState Wise2_PackAln_read_generic_TwoTrackRefSingleState


/* Function:  calculate_TwoTrackRefSingleState(mat)
 *
 * Descrip:    This function calculates the TwoTrackRefSingleState matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_TwoTrackRefSingleState
 *
 *
 * Arg:        mat [UNKN ] TwoTrackRefSingleState which contains explicit basematrix memory [TwoTrackRefSingleState *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_TwoTrackRefSingleState(TwoTrackRefSingleState * mat);
#define calculate_TwoTrackRefSingleState Wise2_calculate_TwoTrackRefSingleState


/* Function:  calculate_dpenv_TwoTrackRefSingleState(mat,dpenv)
 *
 * Descrip:    This function calculates the TwoTrackRefSingleState matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] TwoTrackRefSingleState which contains explicit basematrix memory [TwoTrackRefSingleState *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,DPEnvelope * dpenv);
#define calculate_dpenv_TwoTrackRefSingleState Wise2_calculate_dpenv_TwoTrackRefSingleState


/* Function:  TwoTrackRefSingleState_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackRefSingleState *]
 *
 */
TwoTrackRefSingleState * Wise2_TwoTrackRefSingleState_alloc(void);
#define TwoTrackRefSingleState_alloc Wise2_TwoTrackRefSingleState_alloc


/* Function:  free_TwoTrackRefSingleState(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TwoTrackRefSingleState *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackRefSingleState *]
 *
 */
TwoTrackRefSingleState * Wise2_free_TwoTrackRefSingleState(TwoTrackRefSingleState * obj);
#define free_TwoTrackRefSingleState Wise2_free_TwoTrackRefSingleState


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_TwoTrackRefSingleState_shatter_access_main(TwoTrackRefSingleState * mat,int i,int j,int state);
#define TwoTrackRefSingleState_shatter_access_main Wise2_TwoTrackRefSingleState_shatter_access_main
int Wise2_TwoTrackRefSingleState_shatter_access_special(TwoTrackRefSingleState * mat,int i,int j,int state);
#define TwoTrackRefSingleState_shatter_access_special Wise2_TwoTrackRefSingleState_shatter_access_special
int Wise2_score_only_TwoTrackRefSingleState(TwoTrack* query,ComplexSequence* target ,Score insert,Score delete);
#define score_only_TwoTrackRefSingleState Wise2_score_only_TwoTrackRefSingleState
TwoTrackRefSingleState * Wise2_allocate_TwoTrackRefSingleState_only(TwoTrack* query,ComplexSequence* target ,Score insert,Score delete);
#define allocate_TwoTrackRefSingleState_only Wise2_allocate_TwoTrackRefSingleState_only
void Wise2_init_TwoTrackRefSingleState(TwoTrackRefSingleState * mat);
#define init_TwoTrackRefSingleState Wise2_init_TwoTrackRefSingleState
AlnRange * Wise2_AlnRange_build_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_TwoTrackRefSingleState Wise2_AlnRange_build_TwoTrackRefSingleState
boolean Wise2_read_hidden_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_TwoTrackRefSingleState Wise2_read_hidden_TwoTrackRefSingleState
int Wise2_max_hidden_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_TwoTrackRefSingleState Wise2_max_hidden_TwoTrackRefSingleState
boolean Wise2_read_special_strip_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_TwoTrackRefSingleState Wise2_read_special_strip_TwoTrackRefSingleState
int Wise2_max_special_strip_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_TwoTrackRefSingleState Wise2_max_special_strip_TwoTrackRefSingleState
int Wise2_max_matrix_to_special_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_TwoTrackRefSingleState Wise2_max_matrix_to_special_TwoTrackRefSingleState
void Wise2_calculate_hidden_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_TwoTrackRefSingleState Wise2_calculate_hidden_TwoTrackRefSingleState
void Wise2_init_hidden_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_TwoTrackRefSingleState Wise2_init_hidden_TwoTrackRefSingleState
boolean Wise2_full_dc_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_TwoTrackRefSingleState Wise2_full_dc_TwoTrackRefSingleState
boolean Wise2_do_dc_single_pass_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_TwoTrackRefSingleState Wise2_do_dc_single_pass_TwoTrackRefSingleState
void Wise2_push_dc_at_merge_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_TwoTrackRefSingleState Wise2_push_dc_at_merge_TwoTrackRefSingleState
void Wise2_follow_on_dc_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_TwoTrackRefSingleState Wise2_follow_on_dc_TwoTrackRefSingleState
void Wise2_run_up_dc_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_TwoTrackRefSingleState Wise2_run_up_dc_TwoTrackRefSingleState
void Wise2_init_dc_TwoTrackRefSingleState(TwoTrackRefSingleState * mat);
#define init_dc_TwoTrackRefSingleState Wise2_init_dc_TwoTrackRefSingleState
int Wise2_start_end_find_end_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int * endj);
#define start_end_find_end_TwoTrackRefSingleState Wise2_start_end_find_end_TwoTrackRefSingleState
boolean Wise2_dc_optimised_start_end_calc_TwoTrackRefSingleState(TwoTrackRefSingleState *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_TwoTrackRefSingleState Wise2_dc_optimised_start_end_calc_TwoTrackRefSingleState
void Wise2_init_start_end_linear_TwoTrackRefSingleState(TwoTrackRefSingleState * mat);
#define init_start_end_linear_TwoTrackRefSingleState Wise2_init_start_end_linear_TwoTrackRefSingleState
AlnConvertSet * Wise2_AlnConvertSet_TwoTrackRefSingleState(void);
#define AlnConvertSet_TwoTrackRefSingleState Wise2_AlnConvertSet_TwoTrackRefSingleState
int Wise2_TwoTrackRefSingleState_explicit_access_main(TwoTrackRefSingleState * mat,int i,int j,int state);
#define TwoTrackRefSingleState_explicit_access_main Wise2_TwoTrackRefSingleState_explicit_access_main
int Wise2_TwoTrackRefSingleState_explicit_access_special(TwoTrackRefSingleState * mat,int i,int j,int state);
#define TwoTrackRefSingleState_explicit_access_special Wise2_TwoTrackRefSingleState_explicit_access_special
int Wise2_find_end_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int * ri,int * rj,int * state,boolean * isspecial,TwoTrackRefSingleState_access_func_holder h);
#define find_end_TwoTrackRefSingleState Wise2_find_end_TwoTrackRefSingleState
void Wise2_TwoTrackRefSingleState_debug_show_matrix(TwoTrackRefSingleState * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define TwoTrackRefSingleState_debug_show_matrix Wise2_TwoTrackRefSingleState_debug_show_matrix
int Wise2_max_calc_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,TwoTrackRefSingleState_access_func_holder h);
#define max_calc_TwoTrackRefSingleState Wise2_max_calc_TwoTrackRefSingleState
int Wise2_max_calc_special_TwoTrackRefSingleState(TwoTrackRefSingleState * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,TwoTrackRefSingleState_access_func_holder h);
#define max_calc_special_TwoTrackRefSingleState Wise2_max_calc_special_TwoTrackRefSingleState

#ifdef _cplusplus
}
#endif

#endif
