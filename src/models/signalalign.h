#ifndef DYNAMITEsignalalignHEADERFILE
#define DYNAMITEsignalalignHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dualsignal.h"



struct Wise2_RawSignalMatParaProb {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability single_event_signal;     
    Probability double_event_signal;     
    Probability triple_event_signal;     
    Probability quad_event_signal;   
    Probability event_ext;   
    Probability between_event_signal;    
    Probability between_event_ext;   
    Probability rogue_signal;    
    Probability rogue_ext;   
    Probability dark_kmer;   
    } ;  
/* RawSignalMatParaProb defined */ 
#ifndef DYNAMITE_DEFINED_RawSignalMatParaProb
typedef struct Wise2_RawSignalMatParaProb Wise2_RawSignalMatParaProb;
#define RawSignalMatParaProb Wise2_RawSignalMatParaProb
#define DYNAMITE_DEFINED_RawSignalMatParaProb
#endif


struct Wise2_RawSignalMatParaScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Score single_event_signal;   
    Score double_event_signal;   
    Score triple_event_signal;   
    Score quad_event_signal;     
    Score event_ext;     
    Score between_event_signal;  
    Score between_event_ext;     
    Score rogue_signal;  
    Score rogue_ext;     
    Score dark_kmer;     
    } ;  
/* RawSignalMatParaScore defined */ 
#ifndef DYNAMITE_DEFINED_RawSignalMatParaScore
typedef struct Wise2_RawSignalMatParaScore Wise2_RawSignalMatParaScore;
#define RawSignalMatParaScore Wise2_RawSignalMatParaScore
#define DYNAMITE_DEFINED_RawSignalMatParaScore
#endif


struct Wise2_RawSignalMat {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BaseMatrix * basematrix;     
    ShatterMatrix * shatter;     
    int leni;    
    int lenj;    
    Sequence* seq;   
    RawSignalSeq* signal;    
    SignalMap* sm;   
    RawSignalMatParaScore* para;     
    } ;  
/* RawSignalMat defined */ 
#ifndef DYNAMITE_DEFINED_RawSignalMat
typedef struct Wise2_RawSignalMat Wise2_RawSignalMat;
#define RawSignalMat Wise2_RawSignalMat
#define DYNAMITE_DEFINED_RawSignalMat
#endif


#ifdef PTHREAD
struct thread_pool_holder_RawSignalMat {  
    Sequence* seq;  /* Static query data: never free'd */ 
    RawSignalSeq* signal;   /* Static target data: never free'd */ 
    SignalMap* sm;   
    RawSignalMatParaScore* para;     
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */
struct Wise2_RawSignalMat_access_func_holder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int (*access_main)(RawSignalMat*,int,int,int);   
    int (*access_special)(RawSignalMat*,int,int,int);    
    } ;  
/* RawSignalMat_access_func_holder defined */ 
#ifndef DYNAMITE_DEFINED_RawSignalMat_access_func_holder
typedef struct Wise2_RawSignalMat_access_func_holder Wise2_RawSignalMat_access_func_holder;
#define RawSignalMat_access_func_holder Wise2_RawSignalMat_access_func_holder
#define DYNAMITE_DEFINED_RawSignalMat_access_func_holder
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_RawSignalMatParaProb(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RawSignalMatParaProb *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMatParaProb *]
 *
 */
RawSignalMatParaProb * Wise2_hard_link_RawSignalMatParaProb(RawSignalMatParaProb * obj);
#define hard_link_RawSignalMatParaProb Wise2_hard_link_RawSignalMatParaProb


/* Function:  RawSignalMatParaProb_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMatParaProb *]
 *
 */
RawSignalMatParaProb * Wise2_RawSignalMatParaProb_alloc(void);
#define RawSignalMatParaProb_alloc Wise2_RawSignalMatParaProb_alloc


/* Function:  free_RawSignalMatParaProb(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RawSignalMatParaProb *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMatParaProb *]
 *
 */
RawSignalMatParaProb * Wise2_free_RawSignalMatParaProb(RawSignalMatParaProb * obj);
#define free_RawSignalMatParaProb Wise2_free_RawSignalMatParaProb


/* Function:  hard_link_RawSignalMatParaScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RawSignalMatParaScore *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMatParaScore *]
 *
 */
RawSignalMatParaScore * Wise2_hard_link_RawSignalMatParaScore(RawSignalMatParaScore * obj);
#define hard_link_RawSignalMatParaScore Wise2_hard_link_RawSignalMatParaScore


/* Function:  RawSignalMatParaScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMatParaScore *]
 *
 */
RawSignalMatParaScore * Wise2_RawSignalMatParaScore_alloc(void);
#define RawSignalMatParaScore_alloc Wise2_RawSignalMatParaScore_alloc


/* Function:  free_RawSignalMatParaScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RawSignalMatParaScore *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMatParaScore *]
 *
 */
RawSignalMatParaScore * Wise2_free_RawSignalMatParaScore(RawSignalMatParaScore * obj);
#define free_RawSignalMatParaScore Wise2_free_RawSignalMatParaScore


/* Function:  PackAln_read_Shatter_RawSignalMat(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [RawSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Shatter_RawSignalMat(RawSignalMat * mat);
#define PackAln_read_Shatter_RawSignalMat Wise2_PackAln_read_Shatter_RawSignalMat


/* Function:  calculate_shatter_RawSignalMat(mat,dpenv)
 *
 * Descrip:    This function calculates the RawSignalMat matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [RawSignalMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_shatter_RawSignalMat(RawSignalMat * mat,DPEnvelope * dpenv);
#define calculate_shatter_RawSignalMat Wise2_calculate_shatter_RawSignalMat


/* Function:  search_RawSignalMat(dbsi,out,seq,signal,sm,para)
 *
 * Descrip:    This function makes a database search of RawSignalMat
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:          dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           seq [UNKN ] Undocumented argument [Sequence*]
 * Arg:        signal [UNKN ] Undocumented argument [RawSignalSeq*]
 * Arg:            sm [UNKN ] Undocumented argument [SignalMap*]
 * Arg:          para [UNKN ] Undocumented argument [RawSignalMatParaScore*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_RawSignalMat(DBSearchImpl * dbsi,Hscore * out,Sequence* seq,RawSignalSeq* signal ,SignalMap* sm,RawSignalMatParaScore* para);
#define search_RawSignalMat Wise2_search_RawSignalMat


/* Function:  serial_search_RawSignalMat(out,seq,signal,sm,para)
 *
 * Descrip:    This function makes a database search of RawSignalMat
 *             It is a single processor implementation
 *
 *
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           seq [UNKN ] Undocumented argument [Sequence*]
 * Arg:        signal [UNKN ] Undocumented argument [RawSignalSeq*]
 * Arg:            sm [UNKN ] Undocumented argument [SignalMap*]
 * Arg:          para [UNKN ] Undocumented argument [RawSignalMatParaScore*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_RawSignalMat(Hscore * out,Sequence* seq,RawSignalSeq* signal ,SignalMap* sm,RawSignalMatParaScore* para);
#define serial_search_RawSignalMat Wise2_serial_search_RawSignalMat


/* Function:  PackAln_bestmemory_RawSignalMat(seq,signal,sm,para,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_RawSignalMat
 *
 *
 * Arg:           seq [UNKN ] query data structure [Sequence*]
 * Arg:        signal [UNKN ] target data structure [RawSignalSeq*]
 * Arg:            sm [UNKN ] Resource [SignalMap*]
 * Arg:          para [UNKN ] Resource [RawSignalMatParaScore*]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_RawSignalMat(Sequence* seq,RawSignalSeq* signal ,SignalMap* sm,RawSignalMatParaScore* para,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_RawSignalMat Wise2_PackAln_bestmemory_RawSignalMat


/* Function:  allocate_Expl_RawSignalMat(seq,signal,sm,para,dpri)
 *
 * Descrip:    This function allocates the RawSignalMat structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_RawSignalMat_only
 *
 *
 * Arg:           seq [UNKN ] query data structure [Sequence*]
 * Arg:        signal [UNKN ] target data structure [RawSignalSeq*]
 * Arg:            sm [UNKN ] Resource [SignalMap*]
 * Arg:          para [UNKN ] Resource [RawSignalMatParaScore*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMat *]
 *
 */
RawSignalMat * Wise2_allocate_Expl_RawSignalMat(Sequence* seq,RawSignalSeq* signal ,SignalMap* sm,RawSignalMatParaScore* para,DPRunImpl * dpri);
#define allocate_Expl_RawSignalMat Wise2_allocate_Expl_RawSignalMat


/* Function:  recalculate_PackAln_RawSignalMat(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by RawSignalMat
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [RawSignalMat *]
 *
 */
void Wise2_recalculate_PackAln_RawSignalMat(PackAln * pal,RawSignalMat * mat);
#define recalculate_PackAln_RawSignalMat Wise2_recalculate_PackAln_RawSignalMat


/* Function:  allocate_Small_RawSignalMat(seq,signal,sm,para)
 *
 * Descrip:    This function allocates the RawSignalMat structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_RawSignalMat_only
 *
 *
 * Arg:           seq [UNKN ] query data structure [Sequence*]
 * Arg:        signal [UNKN ] target data structure [RawSignalSeq*]
 * Arg:            sm [UNKN ] Resource [SignalMap*]
 * Arg:          para [UNKN ] Resource [RawSignalMatParaScore*]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMat *]
 *
 */
RawSignalMat * Wise2_allocate_Small_RawSignalMat(Sequence* seq,RawSignalSeq* signal ,SignalMap* sm,RawSignalMatParaScore* para);
#define allocate_Small_RawSignalMat Wise2_allocate_Small_RawSignalMat


/* Function:  PackAln_calculate_Small_RawSignalMat(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for RawSignalMat structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_RawSignalMat 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_RawSignalMat 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_RawSignalMat(RawSignalMat * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_RawSignalMat Wise2_PackAln_calculate_Small_RawSignalMat


/* Function:  AlnRangeSet_calculate_Small_RawSignalMat(mat)
 *
 * Descrip:    This function calculates an alignment for RawSignalMat structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_RawSignalMat 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_RawSignalMat
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_RawSignalMat 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [RawSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_RawSignalMat(RawSignalMat * mat);
#define AlnRangeSet_calculate_Small_RawSignalMat Wise2_AlnRangeSet_calculate_Small_RawSignalMat


/* Function:  AlnRangeSet_from_RawSignalMat(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for RawSignalMat structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_RawSignalMat 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_RawSignalMat
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [RawSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_RawSignalMat(RawSignalMat * mat);
#define AlnRangeSet_from_RawSignalMat Wise2_AlnRangeSet_from_RawSignalMat


/* Function:  convert_PackAln_to_AlnBlock_RawSignalMat(pal)
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
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_RawSignalMat(PackAln * pal);
#define convert_PackAln_to_AlnBlock_RawSignalMat Wise2_convert_PackAln_to_AlnBlock_RawSignalMat


/* Function:  PackAln_read_Expl_RawSignalMat(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [RawSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_RawSignalMat(RawSignalMat * mat);
#define PackAln_read_Expl_RawSignalMat Wise2_PackAln_read_Expl_RawSignalMat


/* Function:  PackAln_read_generic_RawSignalMat(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:          h [UNKN ] Undocumented argument [RawSignalMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_generic_RawSignalMat(RawSignalMat * mat,RawSignalMat_access_func_holder h);
#define PackAln_read_generic_RawSignalMat Wise2_PackAln_read_generic_RawSignalMat


/* Function:  calculate_RawSignalMat(mat)
 *
 * Descrip:    This function calculates the RawSignalMat matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_RawSignalMat
 *
 *
 * Arg:        mat [UNKN ] RawSignalMat which contains explicit basematrix memory [RawSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_RawSignalMat(RawSignalMat * mat);
#define calculate_RawSignalMat Wise2_calculate_RawSignalMat


/* Function:  calculate_dpenv_RawSignalMat(mat,dpenv)
 *
 * Descrip:    This function calculates the RawSignalMat matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] RawSignalMat which contains explicit basematrix memory [RawSignalMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_dpenv_RawSignalMat(RawSignalMat * mat,DPEnvelope * dpenv);
#define calculate_dpenv_RawSignalMat Wise2_calculate_dpenv_RawSignalMat


/* Function:  RawSignalMat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMat *]
 *
 */
RawSignalMat * Wise2_RawSignalMat_alloc(void);
#define RawSignalMat_alloc Wise2_RawSignalMat_alloc


/* Function:  free_RawSignalMat(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RawSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMat *]
 *
 */
RawSignalMat * Wise2_free_RawSignalMat(RawSignalMat * obj);
#define free_RawSignalMat Wise2_free_RawSignalMat


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
SignalEventList * Wise2_implied_event_from_RawSignalSeq_align(Sequence * seq,RawSignalSeq * rss,SignalMap * sm,AlnBlock * alb);
#define implied_event_from_RawSignalSeq_align Wise2_implied_event_from_RawSignalSeq_align
void Wise2_show_fit_RawSignalMat(Sequence * seq,RawSignalSeq * rss,RawSignalMatParaScore * rsmp,SignalMap * sm,AlnBlock * alb,FILE * ofp);
#define show_fit_RawSignalMat Wise2_show_fit_RawSignalMat
void Wise2_show_help_RawSignalMatParaProb_from_argv(FILE * ofp);
#define show_help_RawSignalMatParaProb_from_argv Wise2_show_help_RawSignalMatParaProb_from_argv
RawSignalMatParaScore * Wise2_make_RawSignalMatParaScore(RawSignalMatParaProb * rsmp);
#define make_RawSignalMatParaScore Wise2_make_RawSignalMatParaScore
RawSignalMatParaProb * Wise2_new_RawSignalMatParaProb_from_argv(int * argc,char ** argv);
#define new_RawSignalMatParaProb_from_argv Wise2_new_RawSignalMatParaProb_from_argv


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_RawSignalMat_shatter_access_main(RawSignalMat * mat,int i,int j,int state);
#define RawSignalMat_shatter_access_main Wise2_RawSignalMat_shatter_access_main
int Wise2_RawSignalMat_shatter_access_special(RawSignalMat * mat,int i,int j,int state);
#define RawSignalMat_shatter_access_special Wise2_RawSignalMat_shatter_access_special
void * Wise2_thread_loop_RawSignalMat(void * ptr);
#define thread_loop_RawSignalMat Wise2_thread_loop_RawSignalMat
int Wise2_score_only_RawSignalMat(Sequence* seq,RawSignalSeq* signal ,SignalMap* sm,RawSignalMatParaScore* para);
#define score_only_RawSignalMat Wise2_score_only_RawSignalMat
RawSignalMat * Wise2_allocate_RawSignalMat_only(Sequence* seq,RawSignalSeq* signal ,SignalMap* sm,RawSignalMatParaScore* para);
#define allocate_RawSignalMat_only Wise2_allocate_RawSignalMat_only
void Wise2_init_RawSignalMat(RawSignalMat * mat);
#define init_RawSignalMat Wise2_init_RawSignalMat
AlnRange * Wise2_AlnRange_build_RawSignalMat(RawSignalMat * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_RawSignalMat Wise2_AlnRange_build_RawSignalMat
boolean Wise2_read_hidden_RawSignalMat(RawSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_RawSignalMat Wise2_read_hidden_RawSignalMat
int Wise2_max_hidden_RawSignalMat(RawSignalMat * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_RawSignalMat Wise2_max_hidden_RawSignalMat
boolean Wise2_read_special_strip_RawSignalMat(RawSignalMat * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_RawSignalMat Wise2_read_special_strip_RawSignalMat
int Wise2_max_special_strip_RawSignalMat(RawSignalMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_RawSignalMat Wise2_max_special_strip_RawSignalMat
int Wise2_max_matrix_to_special_RawSignalMat(RawSignalMat * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_RawSignalMat Wise2_max_matrix_to_special_RawSignalMat
void Wise2_calculate_hidden_RawSignalMat(RawSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_RawSignalMat Wise2_calculate_hidden_RawSignalMat
void Wise2_init_hidden_RawSignalMat(RawSignalMat * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_RawSignalMat Wise2_init_hidden_RawSignalMat
boolean Wise2_full_dc_RawSignalMat(RawSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_RawSignalMat Wise2_full_dc_RawSignalMat
boolean Wise2_do_dc_single_pass_RawSignalMat(RawSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_RawSignalMat Wise2_do_dc_single_pass_RawSignalMat
void Wise2_push_dc_at_merge_RawSignalMat(RawSignalMat * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_RawSignalMat Wise2_push_dc_at_merge_RawSignalMat
void Wise2_follow_on_dc_RawSignalMat(RawSignalMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_RawSignalMat Wise2_follow_on_dc_RawSignalMat
void Wise2_run_up_dc_RawSignalMat(RawSignalMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_RawSignalMat Wise2_run_up_dc_RawSignalMat
void Wise2_init_dc_RawSignalMat(RawSignalMat * mat);
#define init_dc_RawSignalMat Wise2_init_dc_RawSignalMat
int Wise2_start_end_find_end_RawSignalMat(RawSignalMat * mat,int * endj);
#define start_end_find_end_RawSignalMat Wise2_start_end_find_end_RawSignalMat
boolean Wise2_dc_optimised_start_end_calc_RawSignalMat(RawSignalMat *mat,DPEnvelope * dpenv);
#define dc_optimised_start_end_calc_RawSignalMat Wise2_dc_optimised_start_end_calc_RawSignalMat
void Wise2_init_start_end_linear_RawSignalMat(RawSignalMat * mat);
#define init_start_end_linear_RawSignalMat Wise2_init_start_end_linear_RawSignalMat
AlnConvertSet * Wise2_AlnConvertSet_RawSignalMat(void);
#define AlnConvertSet_RawSignalMat Wise2_AlnConvertSet_RawSignalMat
int Wise2_RawSignalMat_explicit_access_main(RawSignalMat * mat,int i,int j,int state);
#define RawSignalMat_explicit_access_main Wise2_RawSignalMat_explicit_access_main
int Wise2_RawSignalMat_explicit_access_special(RawSignalMat * mat,int i,int j,int state);
#define RawSignalMat_explicit_access_special Wise2_RawSignalMat_explicit_access_special
int Wise2_find_end_RawSignalMat(RawSignalMat * mat,int * ri,int * rj,int * state,boolean * isspecial,RawSignalMat_access_func_holder h);
#define find_end_RawSignalMat Wise2_find_end_RawSignalMat
void Wise2_RawSignalMat_debug_show_matrix(RawSignalMat * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define RawSignalMat_debug_show_matrix Wise2_RawSignalMat_debug_show_matrix
int Wise2_max_calc_RawSignalMat(RawSignalMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,RawSignalMat_access_func_holder h);
#define max_calc_RawSignalMat Wise2_max_calc_RawSignalMat
int Wise2_max_calc_special_RawSignalMat(RawSignalMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,RawSignalMat_access_func_holder h);
#define max_calc_special_RawSignalMat Wise2_max_calc_special_RawSignalMat

#ifdef _cplusplus
}
#endif

#endif
