#ifndef DYNAMITEfoldbackmatHEADERFILE
#define DYNAMITEfoldbackmatHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"


struct Wise2_FoldBackMat {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    ComplexSequence* query;  
    ComplexSequence* target;     
    DnaMatrix* comp;     
    DnaMatrix* stutter;  
    Score qgap;  
    Score qext;  
    Score tgap;  
    Score text;  
    Score lgap;  
    Score lext;  
    Score stgap;     
    } ;  
/* FoldBackMat defined */ 
#ifndef DYNAMITE_DEFINED_FoldBackMat
typedef struct Wise2_FoldBackMat Wise2_FoldBackMat;
#define FoldBackMat Wise2_FoldBackMat
#define DYNAMITE_DEFINED_FoldBackMat
#endif


struct Wise2_FoldBackMat_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(FoldBackMat*,int,int,int);    
    int (*access_special)(FoldBackMat*,int,int,int); 
    } ;  
/* FoldBackMat_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_FoldBackMat_access_func_holder
typedef struct Wise2_FoldBackMat_access_func_holder Wise2_FoldBackMat_access_func_holder;
#define FoldBackMat_access_func_holder Wise2_FoldBackMat_access_func_holder
#define DYNAMITE_DEFINED_FoldBackMat_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_FoldBackMat(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FoldBackMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_FoldBackMat(FoldBackMat * mat);
#define PackAln_read_Shatter_FoldBackMat Wise2_PackAln_read_Shatter_FoldBackMat


/* Function:  calculate_shatter_FoldBackMat(mat,dpenv)
 *
 * Descrip:    This function calculates the FoldBackMat matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [FoldBackMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_FoldBackMat(FoldBackMat * mat,DPEnvelope * dpenv);
#define calculate_shatter_FoldBackMat Wise2_calculate_shatter_FoldBackMat


/* Function:  search_FoldBackMat(dbsi,out,query,target,comp,stutter,qgap,qext,tgap,text,lgap,lext,stgap)
 *
 * Descrip:    This function makes a database search of FoldBackMat
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:           dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:            out [UNKN ] Undocumented argument [Hscore *]
 * Arg:          query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:         target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:           comp [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        stutter [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:           qgap [UNKN ] Undocumented argument [Score]
 * Arg:           qext [UNKN ] Undocumented argument [Score]
 * Arg:           tgap [UNKN ] Undocumented argument [Score]
 * Arg:           text [UNKN ] Undocumented argument [Score]
 * Arg:           lgap [UNKN ] Undocumented argument [Score]
 * Arg:           lext [UNKN ] Undocumented argument [Score]
 * Arg:          stgap [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_FoldBackMat(DBSearchImpl * dbsi,Hscore * out,ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,DnaMatrix* stutter,Score qgap,Score qext,Score tgap,Score text,Score lgap,Score lext,Score stgap);
#define search_FoldBackMat Wise2_search_FoldBackMat


/* Function:  serial_search_FoldBackMat(out,query,target,comp,stutter,qgap,qext,tgap,text,lgap,lext,stgap)
 *
 * Descrip:    This function makes a database search of FoldBackMat
 *             It is a single processor implementation
 *
 *
 * Arg:            out [UNKN ] Undocumented argument [Hscore *]
 * Arg:          query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:         target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:           comp [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        stutter [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:           qgap [UNKN ] Undocumented argument [Score]
 * Arg:           qext [UNKN ] Undocumented argument [Score]
 * Arg:           tgap [UNKN ] Undocumented argument [Score]
 * Arg:           text [UNKN ] Undocumented argument [Score]
 * Arg:           lgap [UNKN ] Undocumented argument [Score]
 * Arg:           lext [UNKN ] Undocumented argument [Score]
 * Arg:          stgap [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_FoldBackMat(Hscore * out,ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,DnaMatrix* stutter,Score qgap,Score qext,Score tgap,Score text,Score lgap,Score lext,Score stgap);
#define serial_search_FoldBackMat Wise2_serial_search_FoldBackMat


/* Function:  PackAln_bestmemory_FoldBackMat(query,target,comp,stutter,qgap,qext,tgap,text,lgap,lext,stgap,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_FoldBackMat
 *
 *
 * Arg:          query [UNKN ] query data structure [ComplexSequence*]
 * Arg:         target [UNKN ] target data structure [ComplexSequence*]
 * Arg:           comp [UNKN ] Resource [DnaMatrix*]
 * Arg:        stutter [UNKN ] Resource [DnaMatrix*]
 * Arg:           qgap [UNKN ] Resource [Score]
 * Arg:           qext [UNKN ] Resource [Score]
 * Arg:           tgap [UNKN ] Resource [Score]
 * Arg:           text [UNKN ] Resource [Score]
 * Arg:           lgap [UNKN ] Resource [Score]
 * Arg:           lext [UNKN ] Resource [Score]
 * Arg:          stgap [UNKN ] Resource [Score]
 * Arg:          dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:           dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_FoldBackMat(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,DnaMatrix* stutter,Score qgap,Score qext,Score tgap,Score text,Score lgap,Score lext,Score stgap,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_FoldBackMat Wise2_PackAln_bestmemory_FoldBackMat


/* Function:  allocate_Expl_FoldBackMat(query,target,comp,stutter,qgap,qext,tgap,text,lgap,lext,stgap,dpri)
 *
 * Descrip:    This function allocates the FoldBackMat structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_FoldBackMat_only
 *
 *
 * Arg:          query [UNKN ] query data structure [ComplexSequence*]
 * Arg:         target [UNKN ] target data structure [ComplexSequence*]
 * Arg:           comp [UNKN ] Resource [DnaMatrix*]
 * Arg:        stutter [UNKN ] Resource [DnaMatrix*]
 * Arg:           qgap [UNKN ] Resource [Score]
 * Arg:           qext [UNKN ] Resource [Score]
 * Arg:           tgap [UNKN ] Resource [Score]
 * Arg:           text [UNKN ] Resource [Score]
 * Arg:           lgap [UNKN ] Resource [Score]
 * Arg:           lext [UNKN ] Resource [Score]
 * Arg:          stgap [UNKN ] Resource [Score]
 * Arg:           dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [FoldBackMat *]
 *
 */
FoldBackMat * Wise2_allocate_Expl_FoldBackMat(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,DnaMatrix* stutter,Score qgap,Score qext,Score tgap,Score text,Score lgap,Score lext,Score stgap,DPRunImpl * dpri);
#define allocate_Expl_FoldBackMat Wise2_allocate_Expl_FoldBackMat


/* Function:  recalculate_PackAln_FoldBackMat(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by FoldBackMat
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [FoldBackMat *]
 *
 */
void Wise2_recalculate_PackAln_FoldBackMat(PackAln * pal,FoldBackMat * mat);
#define recalculate_PackAln_FoldBackMat Wise2_recalculate_PackAln_FoldBackMat


/* Function:  allocate_Small_FoldBackMat(query,target,comp,stutter,qgap,qext,tgap,text,lgap,lext,stgap)
 *
 * Descrip:    This function allocates the FoldBackMat structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_FoldBackMat_only
 *
 *
 * Arg:          query [UNKN ] query data structure [ComplexSequence*]
 * Arg:         target [UNKN ] target data structure [ComplexSequence*]
 * Arg:           comp [UNKN ] Resource [DnaMatrix*]
 * Arg:        stutter [UNKN ] Resource [DnaMatrix*]
 * Arg:           qgap [UNKN ] Resource [Score]
 * Arg:           qext [UNKN ] Resource [Score]
 * Arg:           tgap [UNKN ] Resource [Score]
 * Arg:           text [UNKN ] Resource [Score]
 * Arg:           lgap [UNKN ] Resource [Score]
 * Arg:           lext [UNKN ] Resource [Score]
 * Arg:          stgap [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [FoldBackMat *]
 *
 */
FoldBackMat * Wise2_allocate_Small_FoldBackMat(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,DnaMatrix* stutter,Score qgap,Score qext,Score tgap,Score text,Score lgap,Score lext,Score stgap);
#define allocate_Small_FoldBackMat Wise2_allocate_Small_FoldBackMat


/* Function:  PackAln_calculate_Small_FoldBackMat(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for FoldBackMat structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_FoldBackMat 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_FoldBackMat 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_FoldBackMat(FoldBackMat * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_FoldBackMat Wise2_PackAln_calculate_Small_FoldBackMat


/* Function:  AlnRangeSet_calculate_Small_FoldBackMat(mat)
 *
 * Descrip:    This function calculates an alignment for FoldBackMat structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_FoldBackMat 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_FoldBackMat
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_FoldBackMat 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FoldBackMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_FoldBackMat(FoldBackMat * mat);
#define AlnRangeSet_calculate_Small_FoldBackMat Wise2_AlnRangeSet_calculate_Small_FoldBackMat


/* Function:  AlnRangeSet_from_FoldBackMat(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for FoldBackMat structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_FoldBackMat 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_FoldBackMat
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FoldBackMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_FoldBackMat(FoldBackMat * mat);
#define AlnRangeSet_from_FoldBackMat Wise2_AlnRangeSet_from_FoldBackMat


/* Function:  convert_PackAln_to_AlnBlock_FoldBackMat(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_FoldBackMat(PackAln * pal);
#define convert_PackAln_to_AlnBlock_FoldBackMat Wise2_convert_PackAln_to_AlnBlock_FoldBackMat


/* Function:  PackAln_read_Expl_FoldBackMat(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FoldBackMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_FoldBackMat(FoldBackMat * mat);
#define PackAln_read_Expl_FoldBackMat Wise2_PackAln_read_Expl_FoldBackMat


/* Function:  PackAln_read_generic_FoldBackMat(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:          h [UNKN ] Undocumented argument [FoldBackMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_FoldBackMat(FoldBackMat * mat,FoldBackMat_access_func_holder h);
#define PackAln_read_generic_FoldBackMat Wise2_PackAln_read_generic_FoldBackMat


/* Function:  calculate_FoldBackMat(mat)
 *
 * Descrip:    This function calculates the FoldBackMat matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_FoldBackMat
 *
 *
 * Arg:        mat [UNKN ] FoldBackMat which contains explicit basematrix memory [FoldBackMat *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_FoldBackMat(FoldBackMat * mat);
#define calculate_FoldBackMat Wise2_calculate_FoldBackMat


/* Function:  calculate_dpenv_FoldBackMat(mat,dpenv)
 *
 * Descrip:    This function calculates the FoldBackMat matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] FoldBackMat which contains explicit basematrix memory [FoldBackMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_FoldBackMat(FoldBackMat * mat,DPEnvelope * dpenv);
#define calculate_dpenv_FoldBackMat Wise2_calculate_dpenv_FoldBackMat


/* Function:  FoldBackMat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FoldBackMat *]
 *
 */
FoldBackMat * Wise2_FoldBackMat_alloc(void);
#define FoldBackMat_alloc Wise2_FoldBackMat_alloc


/* Function:  free_FoldBackMat(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FoldBackMat *]
 *
 * Return [UNKN ]  Undocumented return value [FoldBackMat *]
 *
 */
FoldBackMat * Wise2_free_FoldBackMat(FoldBackMat * obj);
#define free_FoldBackMat Wise2_free_FoldBackMat


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_FoldBackMat_shatter_access_main(FoldBackMat * mat,int i,int j,int state);
#define FoldBackMat_shatter_access_main Wise2_FoldBackMat_shatter_access_main
int Wise2_FoldBackMat_shatter_access_special(FoldBackMat * mat,int i,int j,int state);
#define FoldBackMat_shatter_access_special Wise2_FoldBackMat_shatter_access_special
int Wise2_score_only_FoldBackMat(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,DnaMatrix* stutter,Score qgap,Score qext,Score tgap,Score text,Score lgap,Score lext,Score stgap);
#define score_only_FoldBackMat Wise2_score_only_FoldBackMat
FoldBackMat * Wise2_allocate_FoldBackMat_only(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,DnaMatrix* stutter,Score qgap,Score qext,Score tgap,Score text,Score lgap,Score lext,Score stgap);
#define allocate_FoldBackMat_only Wise2_allocate_FoldBackMat_only
void Wise2_init_FoldBackMat(FoldBackMat * mat);
#define init_FoldBackMat Wise2_init_FoldBackMat
AlnRange * Wise2_AlnRange_build_FoldBackMat(FoldBackMat * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_FoldBackMat Wise2_AlnRange_build_FoldBackMat
boolean Wise2_read_hidden_FoldBackMat(FoldBackMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_FoldBackMat Wise2_read_hidden_FoldBackMat
int Wise2_max_hidden_FoldBackMat(FoldBackMat * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_FoldBackMat Wise2_max_hidden_FoldBackMat
boolean Wise2_read_special_strip_FoldBackMat(FoldBackMat * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_FoldBackMat Wise2_read_special_strip_FoldBackMat
int Wise2_max_special_strip_FoldBackMat(FoldBackMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_FoldBackMat Wise2_max_special_strip_FoldBackMat
int Wise2_max_matrix_to_special_FoldBackMat(FoldBackMat * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_FoldBackMat Wise2_max_matrix_to_special_FoldBackMat
void Wise2_calculate_hidden_FoldBackMat(FoldBackMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_FoldBackMat Wise2_calculate_hidden_FoldBackMat
void Wise2_init_hidden_FoldBackMat(FoldBackMat * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_FoldBackMat Wise2_init_hidden_FoldBackMat
boolean Wise2_full_dc_FoldBackMat(FoldBackMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_FoldBackMat Wise2_full_dc_FoldBackMat
boolean Wise2_do_dc_single_pass_FoldBackMat(FoldBackMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_FoldBackMat Wise2_do_dc_single_pass_FoldBackMat
void Wise2_push_dc_at_merge_FoldBackMat(FoldBackMat * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_FoldBackMat Wise2_push_dc_at_merge_FoldBackMat
void Wise2_follow_on_dc_FoldBackMat(FoldBackMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_FoldBackMat Wise2_follow_on_dc_FoldBackMat
void Wise2_run_up_dc_FoldBackMat(FoldBackMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_FoldBackMat Wise2_run_up_dc_FoldBackMat
void Wise2_init_dc_FoldBackMat(FoldBackMat * mat);
#define init_dc_FoldBackMat Wise2_init_dc_FoldBackMat
int Wise2_start_end_find_end_FoldBackMat(FoldBackMat * mat,int * endj);
#define start_end_find_end_FoldBackMat Wise2_start_end_find_end_FoldBackMat
boolean Wise2_dc_optimised_start_end_calc_FoldBackMat(FoldBackMat *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_FoldBackMat Wise2_dc_optimised_start_end_calc_FoldBackMat
void Wise2_init_start_end_linear_FoldBackMat(FoldBackMat * mat);
#define init_start_end_linear_FoldBackMat Wise2_init_start_end_linear_FoldBackMat
AlnConvertSet * Wise2_AlnConvertSet_FoldBackMat(void);
#define AlnConvertSet_FoldBackMat Wise2_AlnConvertSet_FoldBackMat
int Wise2_FoldBackMat_explicit_access_main(FoldBackMat * mat,int i,int j,int state);
#define FoldBackMat_explicit_access_main Wise2_FoldBackMat_explicit_access_main
int Wise2_FoldBackMat_explicit_access_special(FoldBackMat * mat,int i,int j,int state);
#define FoldBackMat_explicit_access_special Wise2_FoldBackMat_explicit_access_special
int Wise2_find_end_FoldBackMat(FoldBackMat * mat,int * ri,int * rj,int * state,boolean * isspecial,FoldBackMat_access_func_holder h);
#define find_end_FoldBackMat Wise2_find_end_FoldBackMat
void Wise2_FoldBackMat_debug_show_matrix(FoldBackMat * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define FoldBackMat_debug_show_matrix Wise2_FoldBackMat_debug_show_matrix
int Wise2_max_calc_FoldBackMat(FoldBackMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,FoldBackMat_access_func_holder h);
#define max_calc_FoldBackMat Wise2_max_calc_FoldBackMat
int Wise2_max_calc_special_FoldBackMat(FoldBackMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,FoldBackMat_access_func_holder h);
#define max_calc_special_FoldBackMat Wise2_max_calc_special_FoldBackMat

#ifdef _cplusplus
}
#endif

#endif
