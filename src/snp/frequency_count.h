#ifndef DYNAMITEfrequency_countHEADERFILE
#define DYNAMITEfrequency_countHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "genotype.h"
#include "probability.h"


#define LocusGenotypeCountSetLISTLENGTH 128
#define LocusGenotypeCountLISTLENGTH 128


struct Wise2_PopulationGenotypeCount {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int count[BiGenotypeLength];     
    int defined_total;   
    Population * pop;    
    } ;  
/* PopulationGenotypeCount defined */ 
#ifndef DYNAMITE_DEFINED_PopulationGenotypeCount
typedef struct Wise2_PopulationGenotypeCount Wise2_PopulationGenotypeCount;
#define PopulationGenotypeCount Wise2_PopulationGenotypeCount
#define DYNAMITE_DEFINED_PopulationGenotypeCount
#endif


struct Wise2_LocusGenotypeCount {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BiLocus * bl;    
    PopulationGenotypeCount ** pgc;  
    int len;/* len for above pgc  */ 
    int maxlen; /* maxlen for above pgc */ 
    } ;  
/* LocusGenotypeCount defined */ 
#ifndef DYNAMITE_DEFINED_LocusGenotypeCount
typedef struct Wise2_LocusGenotypeCount Wise2_LocusGenotypeCount;
#define LocusGenotypeCount Wise2_LocusGenotypeCount
#define DYNAMITE_DEFINED_LocusGenotypeCount
#endif


struct Wise2_LocusGenotypeCountSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    LocusGenotypeCount ** lgc;   
    int len;/* len for above lgc  */ 
    int maxlen; /* maxlen for above lgc */ 
    } ;  
/* LocusGenotypeCountSet defined */ 
#ifndef DYNAMITE_DEFINED_LocusGenotypeCountSet
typedef struct Wise2_LocusGenotypeCountSet Wise2_LocusGenotypeCountSet;
#define LocusGenotypeCountSet Wise2_LocusGenotypeCountSet
#define DYNAMITE_DEFINED_LocusGenotypeCountSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_PopulationGenotypeCount(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PopulationGenotypeCount *]
 *
 * Return [UNKN ]  Undocumented return value [PopulationGenotypeCount *]
 *
 */
PopulationGenotypeCount * Wise2_hard_link_PopulationGenotypeCount(PopulationGenotypeCount * obj);
#define hard_link_PopulationGenotypeCount Wise2_hard_link_PopulationGenotypeCount


/* Function:  PopulationGenotypeCount_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PopulationGenotypeCount *]
 *
 */
PopulationGenotypeCount * Wise2_PopulationGenotypeCount_alloc(void);
#define PopulationGenotypeCount_alloc Wise2_PopulationGenotypeCount_alloc


/* Function:  free_PopulationGenotypeCount(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PopulationGenotypeCount *]
 *
 * Return [UNKN ]  Undocumented return value [PopulationGenotypeCount *]
 *
 */
PopulationGenotypeCount * Wise2_free_PopulationGenotypeCount(PopulationGenotypeCount * obj);
#define free_PopulationGenotypeCount Wise2_free_PopulationGenotypeCount


/* Function:  add_LocusGenotypeCount(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LocusGenotypeCount *]
 * Arg:        add [OWNER] Object to add to the list [PopulationGenotypeCount *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_LocusGenotypeCount(LocusGenotypeCount * obj,PopulationGenotypeCount * add);
#define add_LocusGenotypeCount Wise2_add_LocusGenotypeCount


/* Function:  flush_LocusGenotypeCount(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LocusGenotypeCount *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_LocusGenotypeCount(LocusGenotypeCount * obj);
#define flush_LocusGenotypeCount Wise2_flush_LocusGenotypeCount


/* Function:  LocusGenotypeCount_alloc_std(void)
 *
 * Descrip:    Equivalent to LocusGenotypeCount_alloc_len(LocusGenotypeCountLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCount *]
 *
 */
LocusGenotypeCount * Wise2_LocusGenotypeCount_alloc_std(void);
#define LocusGenotypeCount_alloc_std Wise2_LocusGenotypeCount_alloc_std


/* Function:  LocusGenotypeCount_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCount *]
 *
 */
LocusGenotypeCount * Wise2_LocusGenotypeCount_alloc_len(int len);
#define LocusGenotypeCount_alloc_len Wise2_LocusGenotypeCount_alloc_len


/* Function:  hard_link_LocusGenotypeCount(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LocusGenotypeCount *]
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCount *]
 *
 */
LocusGenotypeCount * Wise2_hard_link_LocusGenotypeCount(LocusGenotypeCount * obj);
#define hard_link_LocusGenotypeCount Wise2_hard_link_LocusGenotypeCount


/* Function:  LocusGenotypeCount_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCount *]
 *
 */
LocusGenotypeCount * Wise2_LocusGenotypeCount_alloc(void);
#define LocusGenotypeCount_alloc Wise2_LocusGenotypeCount_alloc


/* Function:  free_LocusGenotypeCount(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocusGenotypeCount *]
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCount *]
 *
 */
LocusGenotypeCount * Wise2_free_LocusGenotypeCount(LocusGenotypeCount * obj);
#define free_LocusGenotypeCount Wise2_free_LocusGenotypeCount


/* Function:  add_LocusGenotypeCountSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LocusGenotypeCountSet *]
 * Arg:        add [OWNER] Object to add to the list [LocusGenotypeCount *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_LocusGenotypeCountSet(LocusGenotypeCountSet * obj,LocusGenotypeCount * add);
#define add_LocusGenotypeCountSet Wise2_add_LocusGenotypeCountSet


/* Function:  flush_LocusGenotypeCountSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LocusGenotypeCountSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_LocusGenotypeCountSet(LocusGenotypeCountSet * obj);
#define flush_LocusGenotypeCountSet Wise2_flush_LocusGenotypeCountSet


/* Function:  LocusGenotypeCountSet_alloc_std(void)
 *
 * Descrip:    Equivalent to LocusGenotypeCountSet_alloc_len(LocusGenotypeCountSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCountSet *]
 *
 */
LocusGenotypeCountSet * Wise2_LocusGenotypeCountSet_alloc_std(void);
#define LocusGenotypeCountSet_alloc_std Wise2_LocusGenotypeCountSet_alloc_std


/* Function:  LocusGenotypeCountSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCountSet *]
 *
 */
LocusGenotypeCountSet * Wise2_LocusGenotypeCountSet_alloc_len(int len);
#define LocusGenotypeCountSet_alloc_len Wise2_LocusGenotypeCountSet_alloc_len


/* Function:  hard_link_LocusGenotypeCountSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LocusGenotypeCountSet *]
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCountSet *]
 *
 */
LocusGenotypeCountSet * Wise2_hard_link_LocusGenotypeCountSet(LocusGenotypeCountSet * obj);
#define hard_link_LocusGenotypeCountSet Wise2_hard_link_LocusGenotypeCountSet


/* Function:  LocusGenotypeCountSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCountSet *]
 *
 */
LocusGenotypeCountSet * Wise2_LocusGenotypeCountSet_alloc(void);
#define LocusGenotypeCountSet_alloc Wise2_LocusGenotypeCountSet_alloc


/* Function:  free_LocusGenotypeCountSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocusGenotypeCountSet *]
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCountSet *]
 *
 */
LocusGenotypeCountSet * Wise2_free_LocusGenotypeCountSet(LocusGenotypeCountSet * obj);
#define free_LocusGenotypeCountSet Wise2_free_LocusGenotypeCountSet


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
double Wise2_chisquared_stat_LocusGenotypeCount(LocusGenotypeCount * lgc);
#define chisquared_stat_LocusGenotypeCount Wise2_chisquared_stat_LocusGenotypeCount
double Wise2_chisquared_stat_PopulationGenotypeCount(PopulationGenotypeCount * pgc);
#define chisquared_stat_PopulationGenotypeCount Wise2_chisquared_stat_PopulationGenotypeCount
boolean Wise2_seen_each_genotype_in_all_populations_LocusGenotypeCount(LocusGenotypeCount * lgc);
#define seen_each_genotype_in_all_populations_LocusGenotypeCount Wise2_seen_each_genotype_in_all_populations_LocusGenotypeCount
double Wise2_central_first_frequency_PopulationGenotypeCount(PopulationGenotypeCount * pgc);
#define central_first_frequency_PopulationGenotypeCount Wise2_central_first_frequency_PopulationGenotypeCount
double Wise2_smallest_minor_allele_LocusGenotypeCount(LocusGenotypeCount * lgc) ;
#define smallest_minor_allele_LocusGenotypeCount Wise2_smallest_minor_allele_LocusGenotypeCount
LocusGenotypeCount * Wise2_resampled_LocusGenotypeCount(LocusGenotypeCount * lgc);
#define resampled_LocusGenotypeCount Wise2_resampled_LocusGenotypeCount
LocusGenotypeCountSet * Wise2_LocusGenotypeCountSet_from_BiGenotypeSet(BiGenotypeSet * bgs);
#define LocusGenotypeCountSet_from_BiGenotypeSet Wise2_LocusGenotypeCountSet_from_BiGenotypeSet


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_LocusGenotypeCount(PopulationGenotypeCount ** list,int i,int j) ;
#define swap_LocusGenotypeCount Wise2_swap_LocusGenotypeCount
void Wise2_qsort_LocusGenotypeCount(PopulationGenotypeCount ** list,int left,int right,int (*comp)(PopulationGenotypeCount * ,PopulationGenotypeCount * ));
#define qsort_LocusGenotypeCount Wise2_qsort_LocusGenotypeCount
void Wise2_sort_LocusGenotypeCount(LocusGenotypeCount * obj,int (*comp)(PopulationGenotypeCount *, PopulationGenotypeCount *));
#define sort_LocusGenotypeCount Wise2_sort_LocusGenotypeCount
boolean Wise2_expand_LocusGenotypeCount(LocusGenotypeCount * obj,int len);
#define expand_LocusGenotypeCount Wise2_expand_LocusGenotypeCount
void Wise2_swap_LocusGenotypeCountSet(LocusGenotypeCount ** list,int i,int j) ;
#define swap_LocusGenotypeCountSet Wise2_swap_LocusGenotypeCountSet
void Wise2_qsort_LocusGenotypeCountSet(LocusGenotypeCount ** list,int left,int right,int (*comp)(LocusGenotypeCount * ,LocusGenotypeCount * ));
#define qsort_LocusGenotypeCountSet Wise2_qsort_LocusGenotypeCountSet
void Wise2_sort_LocusGenotypeCountSet(LocusGenotypeCountSet * obj,int (*comp)(LocusGenotypeCount *, LocusGenotypeCount *));
#define sort_LocusGenotypeCountSet Wise2_sort_LocusGenotypeCountSet
boolean Wise2_expand_LocusGenotypeCountSet(LocusGenotypeCountSet * obj,int len);
#define expand_LocusGenotypeCountSet Wise2_expand_LocusGenotypeCountSet

#ifdef _cplusplus
}
#endif

#endif
