#ifndef DYNAMITEdualsignalHEADERFILE
#define DYNAMITEdualsignalHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "dnamapping.h"

#define SignalMapLISTLENGTH 128

#define MAXRANDOM 2147483647

#define MAX_SIGNAL_EMISSION 512

#define MAX_EVENT_KMER 24

#define DUALSIGNAL_INIT_SIZE 1024

#define SignalEventListLISTLENGTH 1024


#define CSEQ_DUAL_SIGNAL_KMER(cseq,j) (ComplexSequence_data(cseq,j,3))

struct Wise2_SignalEvent {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char base;   
    double mean;     
    double std;  
    char kmer[MAX_EVENT_KMER];   
    double time_pos;     
    double time_length;  
    } ;  
/* SignalEvent defined */ 
#ifndef DYNAMITE_DEFINED_SignalEvent
typedef struct Wise2_SignalEvent Wise2_SignalEvent;
#define SignalEvent Wise2_SignalEvent
#define DYNAMITE_DEFINED_SignalEvent
#endif


struct Wise2_SignalEventList {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    SignalEvent ** event;    
    int len;/* len for above event  */ 
    int maxlen; /* maxlen for above event */ 
    } ;  
/* SignalEventList defined */ 
#ifndef DYNAMITE_DEFINED_SignalEventList
typedef struct Wise2_SignalEventList Wise2_SignalEventList;
#define SignalEventList Wise2_SignalEventList
#define DYNAMITE_DEFINED_SignalEventList
#endif


struct Wise2_RawSignalSeq {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    double * signal;     
    double * time_len;   
    double * signal_var;     
    int len;     
    double start_time;   
    char * name;     
    } ;  
/* RawSignalSeq defined */ 
#ifndef DYNAMITE_DEFINED_RawSignalSeq
typedef struct Wise2_RawSignalSeq Wise2_RawSignalSeq;
#define RawSignalSeq Wise2_RawSignalSeq
#define DYNAMITE_DEFINED_RawSignalSeq
#endif


struct Wise2_SignalSeq {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    double * signal;     
    int len;    /*  adjusted for the kmer length */ 
    int true_len;    
    Sequence * called;   
    } ;  
/* SignalSeq defined */ 
#ifndef DYNAMITE_DEFINED_SignalSeq
typedef struct Wise2_SignalSeq Wise2_SignalSeq;
#define SignalSeq Wise2_SignalSeq
#define DYNAMITE_DEFINED_SignalSeq
#endif


struct Wise2_SignalComp {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int kmer;    
    char kseq[20];   
    double mean;     
    double sd;   
    Probability em_prob[MAX_SIGNAL_EMISSION];    
    Score em_score[MAX_SIGNAL_EMISSION];     
    double em_weight[MAX_SIGNAL_EMISSION];   
    } ;  
/* SignalComp defined */ 
#ifndef DYNAMITE_DEFINED_SignalComp
typedef struct Wise2_SignalComp Wise2_SignalComp;
#define SignalComp Wise2_SignalComp
#define DYNAMITE_DEFINED_SignalComp
#endif


struct Wise2_SignalMap {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int kbasis;  
    SignalComp ** comp;  
    int len;     
    int emission_length;     
    double emission_start;   
    double emission_end;     
    double emission_step;    
    Probability em_null[MAX_SIGNAL_EMISSION];    
    } ;  
/* SignalMap defined */ 
#ifndef DYNAMITE_DEFINED_SignalMap
typedef struct Wise2_SignalMap Wise2_SignalMap;
#define SignalMap Wise2_SignalMap
#define DYNAMITE_DEFINED_SignalMap
#endif


struct Wise2_SignalSim {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Probability insert;  
    Probability stay_insert;     
    Probability skip;    
    } ;  
/* SignalSim defined */ 
#ifndef DYNAMITE_DEFINED_SignalSim
typedef struct Wise2_SignalSim Wise2_SignalSim;
#define SignalSim Wise2_SignalSim
#define DYNAMITE_DEFINED_SignalSim
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_SignalEvent(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SignalEvent *]
 *
 * Return [UNKN ]  Undocumented return value [SignalEvent *]
 *
 */
SignalEvent * Wise2_hard_link_SignalEvent(SignalEvent * obj);
#define hard_link_SignalEvent Wise2_hard_link_SignalEvent


/* Function:  SignalEvent_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SignalEvent *]
 *
 */
SignalEvent * Wise2_SignalEvent_alloc(void);
#define SignalEvent_alloc Wise2_SignalEvent_alloc


/* Function:  free_SignalEvent(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SignalEvent *]
 *
 * Return [UNKN ]  Undocumented return value [SignalEvent *]
 *
 */
SignalEvent * Wise2_free_SignalEvent(SignalEvent * obj);
#define free_SignalEvent Wise2_free_SignalEvent


/* Function:  add_SignalEventList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SignalEventList *]
 * Arg:        add [OWNER] Object to add to the list [SignalEvent *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_SignalEventList(SignalEventList * obj,SignalEvent * add);
#define add_SignalEventList Wise2_add_SignalEventList


/* Function:  flush_SignalEventList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SignalEventList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_SignalEventList(SignalEventList * obj);
#define flush_SignalEventList Wise2_flush_SignalEventList


/* Function:  SignalEventList_alloc_std(void)
 *
 * Descrip:    Equivalent to SignalEventList_alloc_len(SignalEventListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SignalEventList *]
 *
 */
SignalEventList * Wise2_SignalEventList_alloc_std(void);
#define SignalEventList_alloc_std Wise2_SignalEventList_alloc_std


/* Function:  SignalEventList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SignalEventList *]
 *
 */
SignalEventList * Wise2_SignalEventList_alloc_len(int len);
#define SignalEventList_alloc_len Wise2_SignalEventList_alloc_len


/* Function:  hard_link_SignalEventList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SignalEventList *]
 *
 * Return [UNKN ]  Undocumented return value [SignalEventList *]
 *
 */
SignalEventList * Wise2_hard_link_SignalEventList(SignalEventList * obj);
#define hard_link_SignalEventList Wise2_hard_link_SignalEventList


/* Function:  SignalEventList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SignalEventList *]
 *
 */
SignalEventList * Wise2_SignalEventList_alloc(void);
#define SignalEventList_alloc Wise2_SignalEventList_alloc


/* Function:  free_SignalEventList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SignalEventList *]
 *
 * Return [UNKN ]  Undocumented return value [SignalEventList *]
 *
 */
SignalEventList * Wise2_free_SignalEventList(SignalEventList * obj);
#define free_SignalEventList Wise2_free_SignalEventList


/* Function:  hard_link_RawSignalSeq(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RawSignalSeq *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalSeq *]
 *
 */
RawSignalSeq * Wise2_hard_link_RawSignalSeq(RawSignalSeq * obj);
#define hard_link_RawSignalSeq Wise2_hard_link_RawSignalSeq


/* Function:  RawSignalSeq_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RawSignalSeq *]
 *
 */
RawSignalSeq * Wise2_RawSignalSeq_alloc(void);
#define RawSignalSeq_alloc Wise2_RawSignalSeq_alloc


/* Function:  free_RawSignalSeq(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RawSignalSeq *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalSeq *]
 *
 */
RawSignalSeq * Wise2_free_RawSignalSeq(RawSignalSeq * obj);
#define free_RawSignalSeq Wise2_free_RawSignalSeq


/* Function:  hard_link_SignalSeq(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SignalSeq *]
 *
 * Return [UNKN ]  Undocumented return value [SignalSeq *]
 *
 */
SignalSeq * Wise2_hard_link_SignalSeq(SignalSeq * obj);
#define hard_link_SignalSeq Wise2_hard_link_SignalSeq


/* Function:  SignalSeq_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SignalSeq *]
 *
 */
SignalSeq * Wise2_SignalSeq_alloc(void);
#define SignalSeq_alloc Wise2_SignalSeq_alloc


/* Function:  free_SignalSeq(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SignalSeq *]
 *
 * Return [UNKN ]  Undocumented return value [SignalSeq *]
 *
 */
SignalSeq * Wise2_free_SignalSeq(SignalSeq * obj);
#define free_SignalSeq Wise2_free_SignalSeq


/* Function:  hard_link_SignalComp(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SignalComp *]
 *
 * Return [UNKN ]  Undocumented return value [SignalComp *]
 *
 */
SignalComp * Wise2_hard_link_SignalComp(SignalComp * obj);
#define hard_link_SignalComp Wise2_hard_link_SignalComp


/* Function:  SignalComp_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SignalComp *]
 *
 */
SignalComp * Wise2_SignalComp_alloc(void);
#define SignalComp_alloc Wise2_SignalComp_alloc


/* Function:  free_SignalComp(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SignalComp *]
 *
 * Return [UNKN ]  Undocumented return value [SignalComp *]
 *
 */
SignalComp * Wise2_free_SignalComp(SignalComp * obj);
#define free_SignalComp Wise2_free_SignalComp


/* Function:  hard_link_SignalMap(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SignalMap *]
 *
 * Return [UNKN ]  Undocumented return value [SignalMap *]
 *
 */
SignalMap * Wise2_hard_link_SignalMap(SignalMap * obj);
#define hard_link_SignalMap Wise2_hard_link_SignalMap


/* Function:  SignalMap_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SignalMap *]
 *
 */
SignalMap * Wise2_SignalMap_alloc(void);
#define SignalMap_alloc Wise2_SignalMap_alloc


/* Function:  free_SignalMap(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SignalMap *]
 *
 * Return [UNKN ]  Undocumented return value [SignalMap *]
 *
 */
SignalMap * Wise2_free_SignalMap(SignalMap * obj);
#define free_SignalMap Wise2_free_SignalMap


/* Function:  hard_link_SignalSim(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SignalSim *]
 *
 * Return [UNKN ]  Undocumented return value [SignalSim *]
 *
 */
SignalSim * Wise2_hard_link_SignalSim(SignalSim * obj);
#define hard_link_SignalSim Wise2_hard_link_SignalSim


/* Function:  SignalSim_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SignalSim *]
 *
 */
SignalSim * Wise2_SignalSim_alloc(void);
#define SignalSim_alloc Wise2_SignalSim_alloc


/* Function:  free_SignalSim(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SignalSim *]
 *
 * Return [UNKN ]  Undocumented return value [SignalSim *]
 *
 */
SignalSim * Wise2_free_SignalSim(SignalSim * obj);
#define free_SignalSim Wise2_free_SignalSim


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
ComplexSequenceEvalSet * Wise2_kmer_ComplexSequenceEvalSet(int kmer_size);
#define kmer_ComplexSequenceEvalSet Wise2_kmer_ComplexSequenceEvalSet
ComplexSequenceEval * Wise2_kmer_number_ComplexSequenceEval(int kmer_size);
#define kmer_number_ComplexSequenceEval Wise2_kmer_number_ComplexSequenceEval
int Wise2_kmer_complexseq_eval(int type,void * data,char * seq);
#define kmer_complexseq_eval Wise2_kmer_complexseq_eval
void Wise2_show_alignment_with_fit_SimpleSignalMat(AlnBlock * alb,SignalEventList * sel,Sequence * comp,SignalMap * sm,FILE * ofp);
#define show_alignment_with_fit_SimpleSignalMat Wise2_show_alignment_with_fit_SimpleSignalMat
SignalMap * Wise2_reverse_SignalMap(SignalMap * sm);
#define reverse_SignalMap Wise2_reverse_SignalMap
Score Wise2_Score_offset_RawSignalMap(SignalMap * sm,RawSignalSeq * raw,int i,Sequence * comp,int j);
#define Score_offset_RawSignalMap Wise2_Score_offset_RawSignalMap
Score Wise2_Score_offset_SignalMap(SignalMap * sm,SignalSeq * sseq,int i,Sequence * comp,int j);
#define Score_offset_SignalMap Wise2_Score_offset_SignalMap
boolean Wise2_flip_coin(Probability p);
#define flip_coin Wise2_flip_coin
Probability Wise2_generate_uniform(void) ;
#define generate_uniform Wise2_generate_uniform
double Wise2_draw_from_normal(double mean,double sd);
#define draw_from_normal Wise2_draw_from_normal
SignalSim * Wise2_new_SignalSim(double insert,double stay_insert,double skip);
#define new_SignalSim Wise2_new_SignalSim
RawSignalSeq * Wise2_simulated_RawSignalSeq_from_Seq(Sequence * seq,SignalMap * sm,SignalSim * sim,int length_of_event);
#define simulated_RawSignalSeq_from_Seq Wise2_simulated_RawSignalSeq_from_Seq
SignalSeq * Wise2_simulated_SignalSeq_from_Seq(Sequence * seq,SignalMap * sm,SignalSim * sim);
#define simulated_SignalSeq_from_Seq Wise2_simulated_SignalSeq_from_Seq
void Wise2_write_RawSignalSeq(RawSignalSeq * rss,FILE *ofp);
#define write_RawSignalSeq Wise2_write_RawSignalSeq
RawSignalSeq * Wise2_read_RawSignalSeq(FILE * ifp);
#define read_RawSignalSeq Wise2_read_RawSignalSeq
SignalSeq * Wise2_read_SignalSeq(FILE * ifp);
#define read_SignalSeq Wise2_read_SignalSeq
void Wise2_write_SignalSeq(SignalSeq * ss,FILE * ofp);
#define write_SignalSeq Wise2_write_SignalSeq
void Wise2_write_SignalMap_absolute(SignalMap * sm,FILE * ofp);
#define write_SignalMap_absolute Wise2_write_SignalMap_absolute
void Wise2_dump_SignalMap_Scores(SignalMap * sm, FILE * ofp);
#define dump_SignalMap_Scores Wise2_dump_SignalMap_Scores
void Wise2_write_SignalMap_normal(SignalMap * sm,FILE * ofp);
#define write_SignalMap_normal Wise2_write_SignalMap_normal
void Wise2_prepare_SignalMap(SignalMap * sm,double pseudocount);
#define prepare_SignalMap Wise2_prepare_SignalMap
void Wise2_convert_normal_to_absolute_SignalMap(SignalMap * sm,double weight);
#define convert_normal_to_absolute_SignalMap Wise2_convert_normal_to_absolute_SignalMap
void Wise2_convert_weight_to_Probability_SignalMap(SignalMap * sm);
#define convert_weight_to_Probability_SignalMap Wise2_convert_weight_to_Probability_SignalMap
SignalMap * Wise2_read_SignalMap_absolute(FILE * ifp);
#define read_SignalMap_absolute Wise2_read_SignalMap_absolute
SignalMap * Wise2_read_SignalMap_tsv(int ksize,FILE * ifp);
#define read_SignalMap_tsv Wise2_read_SignalMap_tsv
boolean Wise2_validate_SignalMap_with_warnings(SignalMap * sm);
#define validate_SignalMap_with_warnings Wise2_validate_SignalMap_with_warnings
SignalMap * Wise2_read_SignalMap_normal(FILE * ifp);
#define read_SignalMap_normal Wise2_read_SignalMap_normal
SignalMap * Wise2_SignalMap_alloc_len(int len);
#define SignalMap_alloc_len Wise2_SignalMap_alloc_len
SignalSeq * Wise2_SignalSeq_from_SignalEventList(SignalEventList * sel);
#define SignalSeq_from_SignalEventList Wise2_SignalSeq_from_SignalEventList
void Wise2_write_SignalEventList(SignalEventList * sel,FILE * ofp);
#define write_SignalEventList Wise2_write_SignalEventList
SignalEventList * Wise2_read_SignalEventList(FILE * ifp);
#define read_SignalEventList Wise2_read_SignalEventList


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_SignalEventList(SignalEvent ** list,int i,int j) ;
#define swap_SignalEventList Wise2_swap_SignalEventList
void Wise2_qsort_SignalEventList(SignalEvent ** list,int left,int right,int (*comp)(SignalEvent * ,SignalEvent * ));
#define qsort_SignalEventList Wise2_qsort_SignalEventList
void Wise2_sort_SignalEventList(SignalEventList * obj,int (*comp)(SignalEvent *, SignalEvent *));
#define sort_SignalEventList Wise2_sort_SignalEventList
boolean Wise2_expand_SignalEventList(SignalEventList * obj,int len);
#define expand_SignalEventList Wise2_expand_SignalEventList

#ifdef _cplusplus
}
#endif

#endif
