#ifndef DYNAMITEsignalalignHEADERFILE
#define DYNAMITEsignalalignHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dualsignal.h"



struct Wise2_SimpleSignalMat {  
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
/* SimpleSignalMat defined */ 
#ifndef DYNAMITE_DEFINED_SimpleSignalMat
typedef struct Wise2_SimpleSignalMat Wise2_SimpleSignalMat;
#define SimpleSignalMat Wise2_SimpleSignalMat
#define DYNAMITE_DEFINED_SimpleSignalMat
#endif


#ifdef PTHREAD
struct thread_pool_holder_SimpleSignalMat {  
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
struct Wise2_SimpleSignalMat_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(SimpleSignalMat*,int,int,int);    
    int (*access_special)(SimpleSignalMat*,int,int,int); 
    } ;  
/* SimpleSignalMat_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_SimpleSignalMat_access_func_holder
typedef struct Wise2_SimpleSignalMat_access_func_holder Wise2_SimpleSignalMat_access_func_holder;
#define SimpleSignalMat_access_func_holder Wise2_SimpleSignalMat_access_func_holder
#define DYNAMITE_DEFINED_SimpleSignalMat_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_SimpleSignalMat(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SimpleSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_SimpleSignalMat(SimpleSignalMat * mat);
#define PackAln_read_Shatter_SimpleSignalMat Wise2_PackAln_read_Shatter_SimpleSignalMat


/* Function:  calculate_shatter_SimpleSignalMat(mat,dpenv)
 *
 * Descrip:    This function calculates the SimpleSignalMat matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [SimpleSignalMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_SimpleSignalMat(SimpleSignalMat * mat,DPEnvelope * dpenv);
#define calculate_shatter_SimpleSignalMat Wise2_calculate_shatter_SimpleSignalMat


/* Function:  search_SimpleSignalMat(dbsi,out,query,target,sm,gap,gapext,seqdiff_open,seqdiff_ext)
 *
 * Descrip:    This function makes a database search of SimpleSignalMat
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
Search_Return_Type Wise2_search_SimpleSignalMat(DBSearchImpl * dbsi,Hscore * out,SignalSeq* query,Sequence* target ,SignalMap* sm,Score gap,Score gapext,Score seqdiff_open,Score seqdiff_ext);
#define search_SimpleSignalMat Wise2_search_SimpleSignalMat


/* Function:  serial_search_SimpleSignalMat(out,query,target,sm,gap,gapext,seqdiff_open,seqdiff_ext)
 *
 * Descrip:    This function makes a database search of SimpleSignalMat
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
Search_Return_Type Wise2_serial_search_SimpleSignalMat(Hscore * out,SignalSeq* query,Sequence* target ,SignalMap* sm,Score gap,Score gapext,Score seqdiff_open,Score seqdiff_ext);
#define serial_search_SimpleSignalMat Wise2_serial_search_SimpleSignalMat


/* Function:  PackAln_bestmemory_SimpleSignalMat(query,target,sm,gap,gapext,seqdiff_open,seqdiff_ext,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_SimpleSignalMat
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
PackAln * Wise2_PackAln_bestmemory_SimpleSignalMat(SignalSeq* query,Sequence* target ,SignalMap* sm,Score gap,Score gapext,Score seqdiff_open,Score seqdiff_ext,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_SimpleSignalMat Wise2_PackAln_bestmemory_SimpleSignalMat


/* Function:  allocate_Expl_SimpleSignalMat(query,target,sm,gap,gapext,seqdiff_open,seqdiff_ext,dpri)
 *
 * Descrip:    This function allocates the SimpleSignalMat structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_SimpleSignalMat_only
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
 * Return [UNKN ]  Undocumented return value [SimpleSignalMat *]
 *
 */
SimpleSignalMat * Wise2_allocate_Expl_SimpleSignalMat(SignalSeq* query,Sequence* target ,SignalMap* sm,Score gap,Score gapext,Score seqdiff_open,Score seqdiff_ext,DPRunImpl * dpri);
#define allocate_Expl_SimpleSignalMat Wise2_allocate_Expl_SimpleSignalMat


/* Function:  recalculate_PackAln_SimpleSignalMat(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by SimpleSignalMat
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [SimpleSignalMat *]
 *
 */
void Wise2_recalculate_PackAln_SimpleSignalMat(PackAln * pal,SimpleSignalMat * mat);
#define recalculate_PackAln_SimpleSignalMat Wise2_recalculate_PackAln_SimpleSignalMat


/* Function:  allocate_Small_SimpleSignalMat(query,target,sm,gap,gapext,seqdiff_open,seqdiff_ext)
 *
 * Descrip:    This function allocates the SimpleSignalMat structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_SimpleSignalMat_only
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
 * Return [UNKN ]  Undocumented return value [SimpleSignalMat *]
 *
 */
SimpleSignalMat * Wise2_allocate_Small_SimpleSignalMat(SignalSeq* query,Sequence* target ,SignalMap* sm,Score gap,Score gapext,Score seqdiff_open,Score seqdiff_ext);
#define allocate_Small_SimpleSignalMat Wise2_allocate_Small_SimpleSignalMat


/* Function:  PackAln_calculate_Small_SimpleSignalMat(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for SimpleSignalMat structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_SimpleSignalMat 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_SimpleSignalMat 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [SimpleSignalMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_SimpleSignalMat(SimpleSignalMat * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_SimpleSignalMat Wise2_PackAln_calculate_Small_SimpleSignalMat


/* Function:  AlnRangeSet_calculate_Small_SimpleSignalMat(mat)
 *
 * Descrip:    This function calculates an alignment for SimpleSignalMat structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_SimpleSignalMat 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_SimpleSignalMat
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_SimpleSignalMat 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SimpleSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_SimpleSignalMat(SimpleSignalMat * mat);
#define AlnRangeSet_calculate_Small_SimpleSignalMat Wise2_AlnRangeSet_calculate_Small_SimpleSignalMat


/* Function:  AlnRangeSet_from_SimpleSignalMat(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for SimpleSignalMat structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_SimpleSignalMat 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_SimpleSignalMat
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SimpleSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_SimpleSignalMat(SimpleSignalMat * mat);
#define AlnRangeSet_from_SimpleSignalMat Wise2_AlnRangeSet_from_SimpleSignalMat


/* Function:  convert_PackAln_to_AlnBlock_SimpleSignalMat(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_SimpleSignalMat(PackAln * pal);
#define convert_PackAln_to_AlnBlock_SimpleSignalMat Wise2_convert_PackAln_to_AlnBlock_SimpleSignalMat


/* Function:  PackAln_read_Expl_SimpleSignalMat(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SimpleSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_SimpleSignalMat(SimpleSignalMat * mat);
#define PackAln_read_Expl_SimpleSignalMat Wise2_PackAln_read_Expl_SimpleSignalMat


/* Function:  PackAln_read_generic_SimpleSignalMat(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [SimpleSignalMat *]
 * Arg:          h [UNKN ] Undocumented argument [SimpleSignalMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_SimpleSignalMat(SimpleSignalMat * mat,SimpleSignalMat_access_func_holder h);
#define PackAln_read_generic_SimpleSignalMat Wise2_PackAln_read_generic_SimpleSignalMat


/* Function:  calculate_SimpleSignalMat(mat)
 *
 * Descrip:    This function calculates the SimpleSignalMat matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_SimpleSignalMat
 *
 *
 * Arg:        mat [UNKN ] SimpleSignalMat which contains explicit basematrix memory [SimpleSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_SimpleSignalMat(SimpleSignalMat * mat);
#define calculate_SimpleSignalMat Wise2_calculate_SimpleSignalMat


/* Function:  calculate_dpenv_SimpleSignalMat(mat,dpenv)
 *
 * Descrip:    This function calculates the SimpleSignalMat matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] SimpleSignalMat which contains explicit basematrix memory [SimpleSignalMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_SimpleSignalMat(SimpleSignalMat * mat,DPEnvelope * dpenv);
#define calculate_dpenv_SimpleSignalMat Wise2_calculate_dpenv_SimpleSignalMat


/* Function:  SimpleSignalMat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SimpleSignalMat *]
 *
 */
SimpleSignalMat * Wise2_SimpleSignalMat_alloc(void);
#define SimpleSignalMat_alloc Wise2_SimpleSignalMat_alloc


/* Function:  free_SimpleSignalMat(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SimpleSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [SimpleSignalMat *]
 *
 */
SimpleSignalMat * Wise2_free_SimpleSignalMat(SimpleSignalMat * obj);
#define free_SimpleSignalMat Wise2_free_SimpleSignalMat


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_show_alignment_with_fit_SimpleSignalMat(AlnBlock * alb,SignalEventList * sel,Sequence * comp,SignalMap * sm,FILE * ofp);
#define show_alignment_with_fit_SimpleSignalMat Wise2_show_alignment_with_fit_SimpleSignalMat


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_SimpleSignalMat_shatter_access_main(SimpleSignalMat * mat,int i,int j,int state);
#define SimpleSignalMat_shatter_access_main Wise2_SimpleSignalMat_shatter_access_main
int Wise2_SimpleSignalMat_shatter_access_special(SimpleSignalMat * mat,int i,int j,int state);
#define SimpleSignalMat_shatter_access_special Wise2_SimpleSignalMat_shatter_access_special
void * Wise2_thread_loop_SimpleSignalMat(void * ptr);
#define thread_loop_SimpleSignalMat Wise2_thread_loop_SimpleSignalMat
int Wise2_score_only_SimpleSignalMat(SignalSeq* query,Sequence* target ,SignalMap* sm,Score gap,Score gapext,Score seqdiff_open,Score seqdiff_ext);
#define score_only_SimpleSignalMat Wise2_score_only_SimpleSignalMat
SimpleSignalMat * Wise2_allocate_SimpleSignalMat_only(SignalSeq* query,Sequence* target ,SignalMap* sm,Score gap,Score gapext,Score seqdiff_open,Score seqdiff_ext);
#define allocate_SimpleSignalMat_only Wise2_allocate_SimpleSignalMat_only
void Wise2_init_SimpleSignalMat(SimpleSignalMat * mat);
#define init_SimpleSignalMat Wise2_init_SimpleSignalMat
AlnRange * Wise2_AlnRange_build_SimpleSignalMat(SimpleSignalMat * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_SimpleSignalMat Wise2_AlnRange_build_SimpleSignalMat
boolean Wise2_read_hidden_SimpleSignalMat(SimpleSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_SimpleSignalMat Wise2_read_hidden_SimpleSignalMat
int Wise2_max_hidden_SimpleSignalMat(SimpleSignalMat * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_SimpleSignalMat Wise2_max_hidden_SimpleSignalMat
boolean Wise2_read_special_strip_SimpleSignalMat(SimpleSignalMat * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_SimpleSignalMat Wise2_read_special_strip_SimpleSignalMat
int Wise2_max_special_strip_SimpleSignalMat(SimpleSignalMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_SimpleSignalMat Wise2_max_special_strip_SimpleSignalMat
int Wise2_max_matrix_to_special_SimpleSignalMat(SimpleSignalMat * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_SimpleSignalMat Wise2_max_matrix_to_special_SimpleSignalMat
void Wise2_calculate_hidden_SimpleSignalMat(SimpleSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_SimpleSignalMat Wise2_calculate_hidden_SimpleSignalMat
void Wise2_init_hidden_SimpleSignalMat(SimpleSignalMat * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_SimpleSignalMat Wise2_init_hidden_SimpleSignalMat
boolean Wise2_full_dc_SimpleSignalMat(SimpleSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_SimpleSignalMat Wise2_full_dc_SimpleSignalMat
boolean Wise2_do_dc_single_pass_SimpleSignalMat(SimpleSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_SimpleSignalMat Wise2_do_dc_single_pass_SimpleSignalMat
void Wise2_push_dc_at_merge_SimpleSignalMat(SimpleSignalMat * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_SimpleSignalMat Wise2_push_dc_at_merge_SimpleSignalMat
void Wise2_follow_on_dc_SimpleSignalMat(SimpleSignalMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_SimpleSignalMat Wise2_follow_on_dc_SimpleSignalMat
void Wise2_run_up_dc_SimpleSignalMat(SimpleSignalMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_SimpleSignalMat Wise2_run_up_dc_SimpleSignalMat
void Wise2_init_dc_SimpleSignalMat(SimpleSignalMat * mat);
#define init_dc_SimpleSignalMat Wise2_init_dc_SimpleSignalMat
int Wise2_start_end_find_end_SimpleSignalMat(SimpleSignalMat * mat,int * endj);
#define start_end_find_end_SimpleSignalMat Wise2_start_end_find_end_SimpleSignalMat
boolean Wise2_dc_optimised_start_end_calc_SimpleSignalMat(SimpleSignalMat *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_SimpleSignalMat Wise2_dc_optimised_start_end_calc_SimpleSignalMat
void Wise2_init_start_end_linear_SimpleSignalMat(SimpleSignalMat * mat);
#define init_start_end_linear_SimpleSignalMat Wise2_init_start_end_linear_SimpleSignalMat
AlnConvertSet * Wise2_AlnConvertSet_SimpleSignalMat(void);
#define AlnConvertSet_SimpleSignalMat Wise2_AlnConvertSet_SimpleSignalMat
int Wise2_SimpleSignalMat_explicit_access_main(SimpleSignalMat * mat,int i,int j,int state);
#define SimpleSignalMat_explicit_access_main Wise2_SimpleSignalMat_explicit_access_main
int Wise2_SimpleSignalMat_explicit_access_special(SimpleSignalMat * mat,int i,int j,int state);
#define SimpleSignalMat_explicit_access_special Wise2_SimpleSignalMat_explicit_access_special
int Wise2_find_end_SimpleSignalMat(SimpleSignalMat * mat,int * ri,int * rj,int * state,boolean * isspecial,SimpleSignalMat_access_func_holder h);
#define find_end_SimpleSignalMat Wise2_find_end_SimpleSignalMat
void Wise2_SimpleSignalMat_debug_show_matrix(SimpleSignalMat * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define SimpleSignalMat_debug_show_matrix Wise2_SimpleSignalMat_debug_show_matrix
int Wise2_max_calc_SimpleSignalMat(SimpleSignalMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,SimpleSignalMat_access_func_holder h);
#define max_calc_SimpleSignalMat Wise2_max_calc_SimpleSignalMat
int Wise2_max_calc_special_SimpleSignalMat(SimpleSignalMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,SimpleSignalMat_access_func_holder h);
#define max_calc_special_SimpleSignalMat Wise2_max_calc_special_SimpleSignalMat

#ifdef _cplusplus
}
#endif

#endif
