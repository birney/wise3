#ifndef DYNAMITElocus_model_estimatorsHEADERFILE
#define DYNAMITElocus_model_estimatorsHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "frequency_count.h"

#define ModelProbabilitySetLISTLENGTH 500

typedef enum LocusModelType {
	LocusModel_Normal_One,
	LocusModel_Normal_Multi,
	LocusModel_SingleHomozygous_First_Selection,
	LocusModel_SingleHomozygous_Second_Selection,
	LocusModel_Hetrozygous_Selection,
	LocusModel_Hemi
} LocusModelType;

#define MAX_POP 6


struct Wise2_ModelProbabilityPoint {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    double estimate_selection_hemi;  
    double estimate_first_freq[MAX_POP];     
    Score likelihood_score;  
    } ;  
/* ModelProbabilityPoint defined */ 
#ifndef DYNAMITE_DEFINED_ModelProbabilityPoint
typedef struct Wise2_ModelProbabilityPoint Wise2_ModelProbabilityPoint;
#define ModelProbabilityPoint Wise2_ModelProbabilityPoint
#define DYNAMITE_DEFINED_ModelProbabilityPoint
#endif


struct Wise2_ModelProbabilitySet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    LocusGenotypeCount * data;   
    LocusModelType type;     
    boolean keep;    
    ModelProbabilityPoint ** mpp;    
    int len;/* len for above mpp  */ 
    int maxlen; /* maxlen for above mpp */ 
    ModelProbabilityPoint * best;    
    double total_likelihood;     
    } ;  
/* ModelProbabilitySet defined */ 
#ifndef DYNAMITE_DEFINED_ModelProbabilitySet
typedef struct Wise2_ModelProbabilitySet Wise2_ModelProbabilitySet;
#define ModelProbabilitySet Wise2_ModelProbabilitySet
#define DYNAMITE_DEFINED_ModelProbabilitySet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_ModelProbabilityPoint(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ModelProbabilityPoint *]
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilityPoint *]
 *
 */
ModelProbabilityPoint * Wise2_hard_link_ModelProbabilityPoint(ModelProbabilityPoint * obj);
#define hard_link_ModelProbabilityPoint Wise2_hard_link_ModelProbabilityPoint


/* Function:  ModelProbabilityPoint_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilityPoint *]
 *
 */
ModelProbabilityPoint * Wise2_ModelProbabilityPoint_alloc(void);
#define ModelProbabilityPoint_alloc Wise2_ModelProbabilityPoint_alloc


/* Function:  free_ModelProbabilityPoint(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ModelProbabilityPoint *]
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilityPoint *]
 *
 */
ModelProbabilityPoint * Wise2_free_ModelProbabilityPoint(ModelProbabilityPoint * obj);
#define free_ModelProbabilityPoint Wise2_free_ModelProbabilityPoint


/* Function:  add_ModelProbabilitySet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ModelProbabilitySet *]
 * Arg:        add [OWNER] Object to add to the list [ModelProbabilityPoint *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ModelProbabilitySet(ModelProbabilitySet * obj,ModelProbabilityPoint * add);
#define add_ModelProbabilitySet Wise2_add_ModelProbabilitySet


/* Function:  flush_ModelProbabilitySet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ModelProbabilitySet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_ModelProbabilitySet(ModelProbabilitySet * obj);
#define flush_ModelProbabilitySet Wise2_flush_ModelProbabilitySet


/* Function:  ModelProbabilitySet_alloc_std(void)
 *
 * Descrip:    Equivalent to ModelProbabilitySet_alloc_len(ModelProbabilitySetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilitySet *]
 *
 */
ModelProbabilitySet * Wise2_ModelProbabilitySet_alloc_std(void);
#define ModelProbabilitySet_alloc_std Wise2_ModelProbabilitySet_alloc_std


/* Function:  ModelProbabilitySet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilitySet *]
 *
 */
ModelProbabilitySet * Wise2_ModelProbabilitySet_alloc_len(int len);
#define ModelProbabilitySet_alloc_len Wise2_ModelProbabilitySet_alloc_len


/* Function:  hard_link_ModelProbabilitySet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ModelProbabilitySet *]
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilitySet *]
 *
 */
ModelProbabilitySet * Wise2_hard_link_ModelProbabilitySet(ModelProbabilitySet * obj);
#define hard_link_ModelProbabilitySet Wise2_hard_link_ModelProbabilitySet


/* Function:  ModelProbabilitySet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilitySet *]
 *
 */
ModelProbabilitySet * Wise2_ModelProbabilitySet_alloc(void);
#define ModelProbabilitySet_alloc Wise2_ModelProbabilitySet_alloc


/* Function:  free_ModelProbabilitySet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ModelProbabilitySet *]
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilitySet *]
 *
 */
ModelProbabilitySet * Wise2_free_ModelProbabilitySet(ModelProbabilitySet * obj);
#define free_ModelProbabilitySet Wise2_free_ModelProbabilitySet


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
ModelProbabilitySet * Wise2_estimate_model_ModelProbabilitySet(LocusGenotypeCount * lgc,double sweep,int points,LocusModelType type,double sel_sweep,int sweep_points);
#define estimate_model_ModelProbabilitySet Wise2_estimate_model_ModelProbabilitySet
void Wise2_estimate_normal_multi(LocusGenotypeCount * lgc,ModelProbabilityPoint * mpp);
#define estimate_normal_multi Wise2_estimate_normal_multi
void Wise2_estimate_homozygous_first_selection_multi(LocusGenotypeCount * lgc,ModelProbabilityPoint * mpp,int show);
#define estimate_homozygous_first_selection_multi Wise2_estimate_homozygous_first_selection_multi
void Wise2_estimate_hetreozygous_selection_multi(LocusGenotypeCount * lgc,ModelProbabilityPoint * mpp);
#define estimate_hetreozygous_selection_multi Wise2_estimate_hetreozygous_selection_multi
void Wise2_estimate_homozygous_second_selection_multi(LocusGenotypeCount * lgc,ModelProbabilityPoint * mpp);
#define estimate_homozygous_second_selection_multi Wise2_estimate_homozygous_second_selection_multi
void Wise2_add_and_update_ModelProbabilitySet(ModelProbabilitySet * out,ModelProbabilityPoint * mpp);
#define add_and_update_ModelProbabilitySet Wise2_add_and_update_ModelProbabilitySet
ModelProbabilitySet * Wise2_estimate_normal_set(LocusGenotypeCount * lgc,double sweep, double points);
#define estimate_normal_set Wise2_estimate_normal_set
ModelProbabilityPoint * Wise2_estimate_normal(LocusGenotypeCount * lgc,double estimate_first_freq);
#define estimate_normal Wise2_estimate_normal


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_ModelProbabilitySet(ModelProbabilityPoint ** list,int i,int j) ;
#define swap_ModelProbabilitySet Wise2_swap_ModelProbabilitySet
void Wise2_qsort_ModelProbabilitySet(ModelProbabilityPoint ** list,int left,int right,int (*comp)(ModelProbabilityPoint * ,ModelProbabilityPoint * ));
#define qsort_ModelProbabilitySet Wise2_qsort_ModelProbabilitySet
void Wise2_sort_ModelProbabilitySet(ModelProbabilitySet * obj,int (*comp)(ModelProbabilityPoint *, ModelProbabilityPoint *));
#define sort_ModelProbabilitySet Wise2_sort_ModelProbabilitySet
boolean Wise2_expand_ModelProbabilitySet(ModelProbabilitySet * obj,int len);
#define expand_ModelProbabilitySet Wise2_expand_ModelProbabilitySet

#ifdef _cplusplus
}
#endif

#endif
