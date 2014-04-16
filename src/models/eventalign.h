#ifndef DYNAMITEeventalignHEADERFILE
#define DYNAMITEeventalignHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dualsignal.h"



struct Wise2_EventSignalMat {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    SignalSeq* query;    
    Sequence* target;    
    SignalMap* sm;   
    Score gap;   
    Score gapext;    
    Score seqdiff_open;  
    Score seqdiff_ext;   
    } ;  
/* EventSignalMat defined */ 
#ifndef DYNAMITE_DEFINED_EventSignalMat
typedef struct Wise2_EventSignalMat Wise2_EventSignalMat;
#define EventSignalMat Wise2_EventSignalMat
#define DYNAMITE_DEFINED_EventSignalMat
#endif


#ifdef PTHREAD
struct thread_pool_holder_EventSignalMat {  
    SignalSeq* query;   /* Static query data: never free'd */ 
    Sequence* target;   /* Static target data: never free'd */ 
    SignalMap* sm;   
    Score gap;   
    Score gapext;    
    Score seqdiff_open;  
    Score seqdiff_ext;   
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_EventSignalMat_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(EventSignalMat*,int,int,int); 
    int (*access_special)(EventSignalMat*,int,int,int);  
    } ;  
/* EventSignalMat_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_EventSignalMat_access_func_holder
typedef struct Wise2_EventSignalMat_access_func_holder Wise2_EventSignalMat_access_func_holder;
#define EventSignalMat_access_func_holder Wise2_EventSignalMat_access_func_holder
#define DYNAMITE_DEFINED_EventSignalMat_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_EventSignalMat(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EventSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_EventSignalMat(EventSignalMat * mat);
#define PackAln_read_Shatter_EventSignalMat Wise2_PackAln_read_Shatter_EventSignalMat


/* Function:  calculate_shatter_EventSignalMat(mat,dpenv)
 *
 * Descrip:    This function calculates the EventSignalMat matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [EventSignalMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_EventSignalMat(EventSignalMat * mat,DPEnvelope * dpenv);
#define calculate_shatter_EventSignalMat Wise2_calculate_shatter_EventSignalMat


/* Function:  search_EventSignalMat(dbsi,out,query,target,sm,gap,gapext,seqdiff_open,seqdiff_ext)
 *
 * Descrip:    This function makes a database search of EventSignalMat
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:                dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:                 out [UNKN ] Undocumented argument [Hscore *]
 * Arg:               query [UNKN ] Undocumented argument [SignalSeq*]
 * Arg:              target [UNKN ] Undocumented argument [Sequence*]
 * Arg:                  sm [UNKN ] Undocumented argument [SignalMap*]
 * Arg:                 gap [UNKN ] Undocumented argument [Score]
 * Arg:              gapext [UNKN ] Undocumented argument [Score]
 * Arg:        seqdiff_open [UNKN ] Undocumented argument [Score]
 * Arg:         seqdiff_ext [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_EventSignalMat(DBSearchImpl * dbsi,Hscore * out,SignalSeq* query,Sequence* target ,SignalMap* sm,Score gap,Score gapext,Score seqdiff_open,Score seqdiff_ext);
#define search_EventSignalMat Wise2_search_EventSignalMat


/* Function:  serial_search_EventSignalMat(out,query,target,sm,gap,gapext,seqdiff_open,seqdiff_ext)
 *
 * Descrip:    This function makes a database search of EventSignalMat
 *             It is a single processor implementation
 *
 *
 * Arg:                 out [UNKN ] Undocumented argument [Hscore *]
 * Arg:               query [UNKN ] Undocumented argument [SignalSeq*]
 * Arg:              target [UNKN ] Undocumented argument [Sequence*]
 * Arg:                  sm [UNKN ] Undocumented argument [SignalMap*]
 * Arg:                 gap [UNKN ] Undocumented argument [Score]
 * Arg:              gapext [UNKN ] Undocumented argument [Score]
 * Arg:        seqdiff_open [UNKN ] Undocumented argument [Score]
 * Arg:         seqdiff_ext [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_EventSignalMat(Hscore * out,SignalSeq* query,Sequence* target ,SignalMap* sm,Score gap,Score gapext,Score seqdiff_open,Score seqdiff_ext);
#define serial_search_EventSignalMat Wise2_serial_search_EventSignalMat


/* Function:  PackAln_bestmemory_EventSignalMat(query,target,sm,gap,gapext,seqdiff_open,seqdiff_ext,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_EventSignalMat
 *
 *
 * Arg:               query [UNKN ] query data structure [SignalSeq*]
 * Arg:              target [UNKN ] target data structure [Sequence*]
 * Arg:                  sm [UNKN ] Resource [SignalMap*]
 * Arg:                 gap [UNKN ] Resource [Score]
 * Arg:              gapext [UNKN ] Resource [Score]
 * Arg:        seqdiff_open [UNKN ] Resource [Score]
 * Arg:         seqdiff_ext [UNKN ] Resource [Score]
 * Arg:               dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:                dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_EventSignalMat(SignalSeq* query,Sequence* target ,SignalMap* sm,Score gap,Score gapext,Score seqdiff_open,Score seqdiff_ext,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_EventSignalMat Wise2_PackAln_bestmemory_EventSignalMat


/* Function:  allocate_Expl_EventSignalMat(query,target,sm,gap,gapext,seqdiff_open,seqdiff_ext,dpri)
 *
 * Descrip:    This function allocates the EventSignalMat structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_EventSignalMat_only
 *
 *
 * Arg:               query [UNKN ] query data structure [SignalSeq*]
 * Arg:              target [UNKN ] target data structure [Sequence*]
 * Arg:                  sm [UNKN ] Resource [SignalMap*]
 * Arg:                 gap [UNKN ] Resource [Score]
 * Arg:              gapext [UNKN ] Resource [Score]
 * Arg:        seqdiff_open [UNKN ] Resource [Score]
 * Arg:         seqdiff_ext [UNKN ] Resource [Score]
 * Arg:                dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [EventSignalMat *]
 *
 */
EventSignalMat * Wise2_allocate_Expl_EventSignalMat(SignalSeq* query,Sequence* target ,SignalMap* sm,Score gap,Score gapext,Score seqdiff_open,Score seqdiff_ext,DPRunImpl * dpri);
#define allocate_Expl_EventSignalMat Wise2_allocate_Expl_EventSignalMat


/* Function:  recalculate_PackAln_EventSignalMat(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by EventSignalMat
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [EventSignalMat *]
 *
 */
void Wise2_recalculate_PackAln_EventSignalMat(PackAln * pal,EventSignalMat * mat);
#define recalculate_PackAln_EventSignalMat Wise2_recalculate_PackAln_EventSignalMat


/* Function:  allocate_Small_EventSignalMat(query,target,sm,gap,gapext,seqdiff_open,seqdiff_ext)
 *
 * Descrip:    This function allocates the EventSignalMat structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_EventSignalMat_only
 *
 *
 * Arg:               query [UNKN ] query data structure [SignalSeq*]
 * Arg:              target [UNKN ] target data structure [Sequence*]
 * Arg:                  sm [UNKN ] Resource [SignalMap*]
 * Arg:                 gap [UNKN ] Resource [Score]
 * Arg:              gapext [UNKN ] Resource [Score]
 * Arg:        seqdiff_open [UNKN ] Resource [Score]
 * Arg:         seqdiff_ext [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [EventSignalMat *]
 *
 */
EventSignalMat * Wise2_allocate_Small_EventSignalMat(SignalSeq* query,Sequence* target ,SignalMap* sm,Score gap,Score gapext,Score seqdiff_open,Score seqdiff_ext);
#define allocate_Small_EventSignalMat Wise2_allocate_Small_EventSignalMat


/* Function:  PackAln_calculate_Small_EventSignalMat(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for EventSignalMat structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_EventSignalMat 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_EventSignalMat 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [EventSignalMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_EventSignalMat(EventSignalMat * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_EventSignalMat Wise2_PackAln_calculate_Small_EventSignalMat


/* Function:  AlnRangeSet_calculate_Small_EventSignalMat(mat)
 *
 * Descrip:    This function calculates an alignment for EventSignalMat structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_EventSignalMat 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_EventSignalMat
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_EventSignalMat 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EventSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_EventSignalMat(EventSignalMat * mat);
#define AlnRangeSet_calculate_Small_EventSignalMat Wise2_AlnRangeSet_calculate_Small_EventSignalMat


/* Function:  AlnRangeSet_from_EventSignalMat(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for EventSignalMat structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_EventSignalMat 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_EventSignalMat
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EventSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_EventSignalMat(EventSignalMat * mat);
#define AlnRangeSet_from_EventSignalMat Wise2_AlnRangeSet_from_EventSignalMat


/* Function:  convert_PackAln_to_AlnBlock_EventSignalMat(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_EventSignalMat(PackAln * pal);
#define convert_PackAln_to_AlnBlock_EventSignalMat Wise2_convert_PackAln_to_AlnBlock_EventSignalMat


/* Function:  PackAln_read_Expl_EventSignalMat(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EventSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_EventSignalMat(EventSignalMat * mat);
#define PackAln_read_Expl_EventSignalMat Wise2_PackAln_read_Expl_EventSignalMat


/* Function:  PackAln_read_generic_EventSignalMat(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [EventSignalMat *]
 * Arg:          h [UNKN ] Undocumented argument [EventSignalMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_EventSignalMat(EventSignalMat * mat,EventSignalMat_access_func_holder h);
#define PackAln_read_generic_EventSignalMat Wise2_PackAln_read_generic_EventSignalMat


/* Function:  calculate_EventSignalMat(mat)
 *
 * Descrip:    This function calculates the EventSignalMat matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_EventSignalMat
 *
 *
 * Arg:        mat [UNKN ] EventSignalMat which contains explicit basematrix memory [EventSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_EventSignalMat(EventSignalMat * mat);
#define calculate_EventSignalMat Wise2_calculate_EventSignalMat


/* Function:  calculate_dpenv_EventSignalMat(mat,dpenv)
 *
 * Descrip:    This function calculates the EventSignalMat matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] EventSignalMat which contains explicit basematrix memory [EventSignalMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_EventSignalMat(EventSignalMat * mat,DPEnvelope * dpenv);
#define calculate_dpenv_EventSignalMat Wise2_calculate_dpenv_EventSignalMat


/* Function:  EventSignalMat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [EventSignalMat *]
 *
 */
EventSignalMat * Wise2_EventSignalMat_alloc(void);
#define EventSignalMat_alloc Wise2_EventSignalMat_alloc


/* Function:  free_EventSignalMat(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [EventSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [EventSignalMat *]
 *
 */
EventSignalMat * Wise2_free_EventSignalMat(EventSignalMat * obj);
#define free_EventSignalMat Wise2_free_EventSignalMat


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_EventSignalMat_shatter_access_main(EventSignalMat * mat,int i,int j,int state);
#define EventSignalMat_shatter_access_main Wise2_EventSignalMat_shatter_access_main
int Wise2_EventSignalMat_shatter_access_special(EventSignalMat * mat,int i,int j,int state);
#define EventSignalMat_shatter_access_special Wise2_EventSignalMat_shatter_access_special
void * Wise2_thread_loop_EventSignalMat(void * ptr);
#define thread_loop_EventSignalMat Wise2_thread_loop_EventSignalMat
int Wise2_score_only_EventSignalMat(SignalSeq* query,Sequence* target ,SignalMap* sm,Score gap,Score gapext,Score seqdiff_open,Score seqdiff_ext);
#define score_only_EventSignalMat Wise2_score_only_EventSignalMat
EventSignalMat * Wise2_allocate_EventSignalMat_only(SignalSeq* query,Sequence* target ,SignalMap* sm,Score gap,Score gapext,Score seqdiff_open,Score seqdiff_ext);
#define allocate_EventSignalMat_only Wise2_allocate_EventSignalMat_only
void Wise2_init_EventSignalMat(EventSignalMat * mat);
#define init_EventSignalMat Wise2_init_EventSignalMat
AlnRange * Wise2_AlnRange_build_EventSignalMat(EventSignalMat * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_EventSignalMat Wise2_AlnRange_build_EventSignalMat
boolean Wise2_read_hidden_EventSignalMat(EventSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_EventSignalMat Wise2_read_hidden_EventSignalMat
int Wise2_max_hidden_EventSignalMat(EventSignalMat * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_EventSignalMat Wise2_max_hidden_EventSignalMat
boolean Wise2_read_special_strip_EventSignalMat(EventSignalMat * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_EventSignalMat Wise2_read_special_strip_EventSignalMat
int Wise2_max_special_strip_EventSignalMat(EventSignalMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_EventSignalMat Wise2_max_special_strip_EventSignalMat
int Wise2_max_matrix_to_special_EventSignalMat(EventSignalMat * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_EventSignalMat Wise2_max_matrix_to_special_EventSignalMat
void Wise2_calculate_hidden_EventSignalMat(EventSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_EventSignalMat Wise2_calculate_hidden_EventSignalMat
void Wise2_init_hidden_EventSignalMat(EventSignalMat * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_EventSignalMat Wise2_init_hidden_EventSignalMat
boolean Wise2_full_dc_EventSignalMat(EventSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_EventSignalMat Wise2_full_dc_EventSignalMat
boolean Wise2_do_dc_single_pass_EventSignalMat(EventSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_EventSignalMat Wise2_do_dc_single_pass_EventSignalMat
void Wise2_push_dc_at_merge_EventSignalMat(EventSignalMat * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_EventSignalMat Wise2_push_dc_at_merge_EventSignalMat
void Wise2_follow_on_dc_EventSignalMat(EventSignalMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_EventSignalMat Wise2_follow_on_dc_EventSignalMat
void Wise2_run_up_dc_EventSignalMat(EventSignalMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_EventSignalMat Wise2_run_up_dc_EventSignalMat
void Wise2_init_dc_EventSignalMat(EventSignalMat * mat);
#define init_dc_EventSignalMat Wise2_init_dc_EventSignalMat
int Wise2_start_end_find_end_EventSignalMat(EventSignalMat * mat,int * endj);
#define start_end_find_end_EventSignalMat Wise2_start_end_find_end_EventSignalMat
boolean Wise2_dc_optimised_start_end_calc_EventSignalMat(EventSignalMat *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_EventSignalMat Wise2_dc_optimised_start_end_calc_EventSignalMat
void Wise2_init_start_end_linear_EventSignalMat(EventSignalMat * mat);
#define init_start_end_linear_EventSignalMat Wise2_init_start_end_linear_EventSignalMat
AlnConvertSet * Wise2_AlnConvertSet_EventSignalMat(void);
#define AlnConvertSet_EventSignalMat Wise2_AlnConvertSet_EventSignalMat
int Wise2_EventSignalMat_explicit_access_main(EventSignalMat * mat,int i,int j,int state);
#define EventSignalMat_explicit_access_main Wise2_EventSignalMat_explicit_access_main
int Wise2_EventSignalMat_explicit_access_special(EventSignalMat * mat,int i,int j,int state);
#define EventSignalMat_explicit_access_special Wise2_EventSignalMat_explicit_access_special
int Wise2_find_end_EventSignalMat(EventSignalMat * mat,int * ri,int * rj,int * state,boolean * isspecial,EventSignalMat_access_func_holder h);
#define find_end_EventSignalMat Wise2_find_end_EventSignalMat
void Wise2_EventSignalMat_debug_show_matrix(EventSignalMat * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define EventSignalMat_debug_show_matrix Wise2_EventSignalMat_debug_show_matrix
int Wise2_max_calc_EventSignalMat(EventSignalMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,EventSignalMat_access_func_holder h);
#define max_calc_EventSignalMat Wise2_max_calc_EventSignalMat
int Wise2_max_calc_special_EventSignalMat(EventSignalMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,EventSignalMat_access_func_holder h);
#define max_calc_special_EventSignalMat Wise2_max_calc_special_EventSignalMat

#ifdef _cplusplus
}
#endif

#endif
