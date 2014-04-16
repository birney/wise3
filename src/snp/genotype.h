#ifndef DYNAMITEgenotypeHEADERFILE
#define DYNAMITEgenotypeHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "locus_framework.h"
#include "person.h"


#define BiGenotypeSetLISTLENGTH     128
#define BiGenotypeLocusLISTLENGTH  128

#define MAX_HAPMAP_LINE 10000

typedef enum PhaseStatus {
	PhaseStatus_Unknown = 45,
	PhaseStatus_FirstUp,
	PhaseStatus_FirstDown,
} PhaseStatus;


typedef struct BiGenotype {
  BiGenotypeType type;
  PhaseStatus phase;
  Person * person;
} BiGenotype;

#define MAX_GENOTYPE_ALLOC 1000000
#define GENOTYPE_ALLOC_LINE 4096

typedef struct BiGenotypeAllocSlab {
  BiGenotype * alloc_array[MAX_GENOTYPE_ALLOC];
  int array_pos;
  int array_y;
} BiGenotypeAllocSlab;



struct Wise2_BiGenotypeLocus {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BiGenotype ** big;   
    int len;/* len for above big  */ 
    int maxlen; /* maxlen for above big */ 
    BiLocus * locus;     
    } ;  
/* BiGenotypeLocus defined */ 
#ifndef DYNAMITE_DEFINED_BiGenotypeLocus
typedef struct Wise2_BiGenotypeLocus Wise2_BiGenotypeLocus;
#define BiGenotypeLocus Wise2_BiGenotypeLocus
#define DYNAMITE_DEFINED_BiGenotypeLocus
#endif


struct Wise2_BiGenotypeSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BiGenotypeLocus ** locus;    
    int len;/* len for above locus  */ 
    int maxlen; /* maxlen for above locus */ 
    BiLocusFramework * framework;    
    PopulationSet * ps;  
    BiGenotypeAllocSlab * slab;  
    } ;  
/* BiGenotypeSet defined */ 
#ifndef DYNAMITE_DEFINED_BiGenotypeSet
typedef struct Wise2_BiGenotypeSet Wise2_BiGenotypeSet;
#define BiGenotypeSet Wise2_BiGenotypeSet
#define DYNAMITE_DEFINED_BiGenotypeSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  add_BiGenotypeLocus(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [BiGenotypeLocus *]
 * Arg:        add [OWNER] Object to add to the list [BiGenotype *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_BiGenotypeLocus(BiGenotypeLocus * obj,BiGenotype * add);
#define add_BiGenotypeLocus Wise2_add_BiGenotypeLocus


/* Function:  flush_BiGenotypeLocus(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [BiGenotypeLocus *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_BiGenotypeLocus(BiGenotypeLocus * obj);
#define flush_BiGenotypeLocus Wise2_flush_BiGenotypeLocus


/* Function:  BiGenotypeLocus_alloc_std(void)
 *
 * Descrip:    Equivalent to BiGenotypeLocus_alloc_len(BiGenotypeLocusLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeLocus *]
 *
 */
BiGenotypeLocus * Wise2_BiGenotypeLocus_alloc_std(void);
#define BiGenotypeLocus_alloc_std Wise2_BiGenotypeLocus_alloc_std


/* Function:  BiGenotypeLocus_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeLocus *]
 *
 */
BiGenotypeLocus * Wise2_BiGenotypeLocus_alloc_len(int len);
#define BiGenotypeLocus_alloc_len Wise2_BiGenotypeLocus_alloc_len


/* Function:  hard_link_BiGenotypeLocus(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [BiGenotypeLocus *]
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeLocus *]
 *
 */
BiGenotypeLocus * Wise2_hard_link_BiGenotypeLocus(BiGenotypeLocus * obj);
#define hard_link_BiGenotypeLocus Wise2_hard_link_BiGenotypeLocus


/* Function:  BiGenotypeLocus_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeLocus *]
 *
 */
BiGenotypeLocus * Wise2_BiGenotypeLocus_alloc(void);
#define BiGenotypeLocus_alloc Wise2_BiGenotypeLocus_alloc


/* Function:  free_BiGenotypeLocus(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [BiGenotypeLocus *]
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeLocus *]
 *
 */
BiGenotypeLocus * Wise2_free_BiGenotypeLocus(BiGenotypeLocus * obj);
#define free_BiGenotypeLocus Wise2_free_BiGenotypeLocus


/* Function:  add_BiGenotypeSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [BiGenotypeSet *]
 * Arg:        add [OWNER] Object to add to the list [BiGenotypeLocus *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_BiGenotypeSet(BiGenotypeSet * obj,BiGenotypeLocus * add);
#define add_BiGenotypeSet Wise2_add_BiGenotypeSet


/* Function:  flush_BiGenotypeSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [BiGenotypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_BiGenotypeSet(BiGenotypeSet * obj);
#define flush_BiGenotypeSet Wise2_flush_BiGenotypeSet


/* Function:  BiGenotypeSet_alloc_std(void)
 *
 * Descrip:    Equivalent to BiGenotypeSet_alloc_len(BiGenotypeSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeSet *]
 *
 */
BiGenotypeSet * Wise2_BiGenotypeSet_alloc_std(void);
#define BiGenotypeSet_alloc_std Wise2_BiGenotypeSet_alloc_std


/* Function:  BiGenotypeSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeSet *]
 *
 */
BiGenotypeSet * Wise2_BiGenotypeSet_alloc_len(int len);
#define BiGenotypeSet_alloc_len Wise2_BiGenotypeSet_alloc_len


/* Function:  hard_link_BiGenotypeSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [BiGenotypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeSet *]
 *
 */
BiGenotypeSet * Wise2_hard_link_BiGenotypeSet(BiGenotypeSet * obj);
#define hard_link_BiGenotypeSet Wise2_hard_link_BiGenotypeSet


/* Function:  BiGenotypeSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeSet *]
 *
 */
BiGenotypeSet * Wise2_BiGenotypeSet_alloc(void);
#define BiGenotypeSet_alloc Wise2_BiGenotypeSet_alloc


/* Function:  free_BiGenotypeSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [BiGenotypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeSet *]
 *
 */
BiGenotypeSet * Wise2_free_BiGenotypeSet(BiGenotypeSet * obj);
#define free_BiGenotypeSet Wise2_free_BiGenotypeSet


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
BiGenotype * Wise2_BiGenotype_alloc_slab(BiGenotypeAllocSlab * slab);
#define BiGenotype_alloc_slab Wise2_BiGenotype_alloc_slab
BiGenotypeAllocSlab * Wise2_new_BiGenotypeAllocSlab(void);
#define new_BiGenotypeAllocSlab Wise2_new_BiGenotypeAllocSlab
BiGenotypeAllocSlab * Wise2_free_BiGenotypeAllocSlab(BiGenotypeAllocSlab * s);
#define free_BiGenotypeAllocSlab Wise2_free_BiGenotypeAllocSlab
BiGenotype * Wise2_BiGenotype_alloc(void);
#define BiGenotype_alloc Wise2_BiGenotype_alloc
BiGenotype * Wise2_free_BiGenotype(BiGenotype * bi);
#define free_BiGenotype Wise2_free_BiGenotype
void Wise2_write_simple_genotype(BiGenotypeSet * bgs,FILE * ofp);
#define write_simple_genotype Wise2_write_simple_genotype
boolean Wise2_read_hapmap_genotype_file(BiGenotypeSet * bgs,char * pop_name,FILE * ifp);
#define read_hapmap_genotype_file Wise2_read_hapmap_genotype_file
BiGenotypeLocus * Wise2_find_or_new_BiGenotypeLocus_from_BiGenotypeSet(BiGenotypeSet * bgs,BiLocus * bi);
#define find_or_new_BiGenotypeLocus_from_BiGenotypeSet Wise2_find_or_new_BiGenotypeLocus_from_BiGenotypeSet
BiGenotypeSet * Wise2_read_simple_genotype_file(FILE * ifp);
#define read_simple_genotype_file Wise2_read_simple_genotype_file
BiGenotypeSet * Wise2_new_BiGenotypeSet(void);
#define new_BiGenotypeSet Wise2_new_BiGenotypeSet


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_BiGenotypeLocus(BiGenotype ** list,int i,int j) ;
#define swap_BiGenotypeLocus Wise2_swap_BiGenotypeLocus
void Wise2_qsort_BiGenotypeLocus(BiGenotype ** list,int left,int right,int (*comp)(BiGenotype * ,BiGenotype * ));
#define qsort_BiGenotypeLocus Wise2_qsort_BiGenotypeLocus
void Wise2_sort_BiGenotypeLocus(BiGenotypeLocus * obj,int (*comp)(BiGenotype *, BiGenotype *));
#define sort_BiGenotypeLocus Wise2_sort_BiGenotypeLocus
boolean Wise2_expand_BiGenotypeLocus(BiGenotypeLocus * obj,int len);
#define expand_BiGenotypeLocus Wise2_expand_BiGenotypeLocus
void Wise2_swap_BiGenotypeSet(BiGenotypeLocus ** list,int i,int j) ;
#define swap_BiGenotypeSet Wise2_swap_BiGenotypeSet
void Wise2_qsort_BiGenotypeSet(BiGenotypeLocus ** list,int left,int right,int (*comp)(BiGenotypeLocus * ,BiGenotypeLocus * ));
#define qsort_BiGenotypeSet Wise2_qsort_BiGenotypeSet
void Wise2_sort_BiGenotypeSet(BiGenotypeSet * obj,int (*comp)(BiGenotypeLocus *, BiGenotypeLocus *));
#define sort_BiGenotypeSet Wise2_sort_BiGenotypeSet
boolean Wise2_expand_BiGenotypeSet(BiGenotypeSet * obj,int len);
#define expand_BiGenotypeSet Wise2_expand_BiGenotypeSet

#ifdef _cplusplus
}
#endif

#endif
