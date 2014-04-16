#ifndef DYNAMITEproteinswHEADERFILE
#define DYNAMITEproteinswHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

struct ProteinSW {  
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    ComplexSequence* query;  
    ComplexSequence* target;     
    CompMat* comp;   
    int gap;     
    int ext;     
    } ;  
/* ProteinSW defined */ 
#ifndef DYNAMITE_DEFINED_ProteinSW
typedef struct ProteinSW ProteinSW;
#define DYNAMITE_DEFINED_ProteinSW
#endif


struct naive_vector_ProteinSW {  
    ComplexSequence* query;  
    ComplexSequence* target; 
    CompMat* comp;   
    int gap;     
    int ext;     
    int ** non_shift;    
    int ** shift;    
    } ;  
struct naive_vector_calc_ProteinSW {  
    int ** state_MATCH_source_MATCH; 
    int ** state_MATCH_source_INSERT;    
    int ** state_MATCH_source_DELETE;    
    int ** state_MATCH_source_START; 
    int ** state_INSERT_source_MATCH;    
    int ** state_INSERT_source_INSERT;   
    int ** state_DELETE_source_MATCH;    
    int ** state_DELETE_source_DELETE;   
    } ;  
struct ProteinSW_access_func_holder {  
    int (*access_main)(ProteinSW*,int,int,int);  
    int (*access_special)(ProteinSW*,int,int,int);   
    } ;  
/* ProteinSW_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_ProteinSW_access_func_holder
typedef struct ProteinSW_access_func_holder ProteinSW_access_func_holder;
#define DYNAMITE_DEFINED_ProteinSW_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  PackAln_read_Shatter_ProteinSW(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinSW *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_ProteinSW(ProteinSW * mat);


/* Function:  calculate_shatter_ProteinSW(mat,dpenv)
 *
 * Descrip:    This function calculates the ProteinSW matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [ProteinSW *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_ProteinSW(ProteinSW * mat,DPEnvelope * dpenv);


/* Function:  search_ProteinSW(dbsi,out,querydb,targetdb,comp,gap,ext)
 *
 * Descrip:    This function makes a database search of ProteinSW
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:            dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         querydb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:        targetdb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:            comp [UNKN ] Undocumented argument [CompMat*]
 * Arg:             gap [UNKN ] Undocumented argument [int]
 * Arg:             ext [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_ProteinSW(DBSearchImpl * dbsi,Hscore * out,ProteinDB* querydb,ProteinDB* targetdb ,CompMat* comp,int gap,int ext);


/* Function:  serial_search_ProteinSW(out,querydb,targetdb,comp,gap,ext)
 *
 * Descrip:    This function makes a database search of ProteinSW
 *             It is a single processor implementation
 *
 *
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         querydb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:        targetdb [UNKN ] Undocumented argument [ProteinDB*]
 * Arg:            comp [UNKN ] Undocumented argument [CompMat*]
 * Arg:             gap [UNKN ] Undocumented argument [int]
 * Arg:             ext [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type serial_search_ProteinSW(Hscore * out,ProteinDB* querydb,ProteinDB* targetdb ,CompMat* comp,int gap,int ext);


/* Function:  PackAln_bestmemory_ProteinSW(query,target,comp,gap,ext,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_ProteinSW
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          comp [UNKN ] Resource [CompMat*]
 * Arg:           gap [UNKN ] Resource [int]
 * Arg:           ext [UNKN ] Resource [int]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_ProteinSW(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int gap,int ext,DPEnvelope * dpenv,DPRunImpl * dpri);


/* Function:  allocate_Expl_ProteinSW(query,target,comp,gap,ext,dpri)
 *
 * Descrip:    This function allocates the ProteinSW structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_ProteinSW_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          comp [UNKN ] Resource [CompMat*]
 * Arg:           gap [UNKN ] Resource [int]
 * Arg:           ext [UNKN ] Resource [int]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinSW *]
 *
 */
ProteinSW * allocate_Expl_ProteinSW(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int gap,int ext,DPRunImpl * dpri);


/* Function:  recalculate_PackAln_ProteinSW(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by ProteinSW
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [ProteinSW *]
 *
 */
void recalculate_PackAln_ProteinSW(PackAln * pal,ProteinSW * mat);


/* Function:  allocate_Small_ProteinSW(query,target,comp,gap,ext)
 *
 * Descrip:    This function allocates the ProteinSW structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_ProteinSW_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:          comp [UNKN ] Resource [CompMat*]
 * Arg:           gap [UNKN ] Resource [int]
 * Arg:           ext [UNKN ] Resource [int]
 *
 * Return [UNKN ]  Undocumented return value [ProteinSW *]
 *
 */
ProteinSW * allocate_Small_ProteinSW(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int gap,int ext);


/* Function:  PackAln_calculate_Small_ProteinSW(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for ProteinSW structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_ProteinSW 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_ProteinSW 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [ProteinSW *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_ProteinSW(ProteinSW * mat,DPEnvelope * dpenv);


/* Function:  AlnRangeSet_calculate_Small_ProteinSW(mat)
 *
 * Descrip:    This function calculates an alignment for ProteinSW structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_ProteinSW 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_ProteinSW
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_ProteinSW 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinSW *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_ProteinSW(ProteinSW * mat);


/* Function:  AlnRangeSet_from_ProteinSW(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for ProteinSW structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_ProteinSW 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_ProteinSW
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinSW *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_ProteinSW(ProteinSW * mat);


/* Function:  convert_PackAln_to_AlnBlock_ProteinSW(pal)
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
AlnBlock * convert_PackAln_to_AlnBlock_ProteinSW(PackAln * pal);


/* Function:  PackAln_read_Expl_ProteinSW(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinSW *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_ProteinSW(ProteinSW * mat);


/* Function:  PackAln_read_generic_ProteinSW(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [ProteinSW *]
 * Arg:          h [UNKN ] Undocumented argument [ProteinSW_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_ProteinSW(ProteinSW * mat,ProteinSW_access_func_holder h);


/* Function:  calculate_ProteinSW(mat)
 *
 * Descrip:    This function calculates the ProteinSW matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_ProteinSW
 *
 *
 * Arg:        mat [UNKN ] ProteinSW which contains explicit basematrix memory [ProteinSW *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_ProteinSW(ProteinSW * mat);


/* Function:  calculate_dpenv_ProteinSW(mat,dpenv)
 *
 * Descrip:    This function calculates the ProteinSW matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] ProteinSW which contains explicit basematrix memory [ProteinSW *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_ProteinSW(ProteinSW * mat,DPEnvelope * dpenv);


/* Function:  ProteinSW_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ProteinSW *]
 *
 */
ProteinSW * ProteinSW_alloc(void);


/* Function:  free_ProteinSW(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ProteinSW *]
 *
 * Return [UNKN ]  Undocumented return value [ProteinSW *]
 *
 */
ProteinSW * free_ProteinSW(ProteinSW * obj);


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int ProteinSW_shatter_access_main(ProteinSW * mat,int i,int j,int state);
int ProteinSW_shatter_access_special(ProteinSW * mat,int i,int j,int state);
struct naive_vector_ProteinSW * alloc_naive_vector_ProteinSW(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int gap,int ext);
int score_only_ProteinSW(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int gap,int ext);
ProteinSW * allocate_ProteinSW_only(ComplexSequence* query,ComplexSequence* target ,CompMat* comp,int gap,int ext);
void init_ProteinSW(ProteinSW * mat);
AlnRange * AlnRange_build_ProteinSW(ProteinSW * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
boolean read_hidden_ProteinSW(ProteinSW * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
int max_hidden_ProteinSW(ProteinSW * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
boolean read_special_strip_ProteinSW(ProteinSW * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
int max_special_strip_ProteinSW(ProteinSW * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
int max_matrix_to_special_ProteinSW(ProteinSW * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
void calculate_hidden_ProteinSW(ProteinSW * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
void init_hidden_ProteinSW(ProteinSW * mat,int starti,int startj,int stopi,int stopj);
boolean full_dc_ProteinSW(ProteinSW * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
boolean do_dc_single_pass_ProteinSW(ProteinSW * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
void push_dc_at_merge_ProteinSW(ProteinSW * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
void follow_on_dc_ProteinSW(ProteinSW * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
void run_up_dc_ProteinSW(ProteinSW * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
void init_dc_ProteinSW(ProteinSW * mat);
int start_end_find_end_ProteinSW(ProteinSW * mat,int * endj);
boolean dc_optimised_start_end_calc_ProteinSW(ProteinSW *mat,DPEnvelope * dpenv);
void init_start_end_linear_ProteinSW(ProteinSW * mat);
AlnConvertSet * AlnConvertSet_ProteinSW(void);
int ProteinSW_explicit_access_main(ProteinSW * mat,int i,int j,int state);
int ProteinSW_explicit_access_special(ProteinSW * mat,int i,int j,int state);
int find_end_ProteinSW(ProteinSW * mat,int * ri,int * rj,int * state,boolean * isspecial,ProteinSW_access_func_holder h);
void ProteinSW_debug_show_matrix(ProteinSW * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
int max_calc_ProteinSW(ProteinSW * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,ProteinSW_access_func_holder h);
int max_calc_special_ProteinSW(ProteinSW * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,ProteinSW_access_func_holder h);

#ifdef _cplusplus
}
#endif

#endif
