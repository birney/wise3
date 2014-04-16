#ifndef DYNAMITElocus_frameworkHEADERFILE
#define DYNAMITElocus_frameworkHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"
#include "glib.h"

typedef enum BiGenotypeType {
	BiGenotypeHomozygousFirst = 0,
	BiGenotypeHetrozygous,
	BiGenotypeHomozygousSecond,
	BiGenotypeUnknown,
	BiGenotypeError,
	BiGenotypeLength
} BiGenotypeType;


#define BiLocusFrameworkLISTLENGTH 128

struct Wise2_BiLocus {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * locus_id;     
    char first_allele_char;  
    char second_allele_char;     
    } ;  
/* BiLocus defined */ 
#ifndef DYNAMITE_DEFINED_BiLocus
typedef struct Wise2_BiLocus Wise2_BiLocus;
#define BiLocus Wise2_BiLocus
#define DYNAMITE_DEFINED_BiLocus
#endif


struct Wise2_BiLocusFramework {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BiLocus ** locus;    
    int len;/* len for above locus  */ 
    int maxlen; /* maxlen for above locus */ 
    GHashTable * hash;   
    } ;  
/* BiLocusFramework defined */ 
#ifndef DYNAMITE_DEFINED_BiLocusFramework
typedef struct Wise2_BiLocusFramework Wise2_BiLocusFramework;
#define BiLocusFramework Wise2_BiLocusFramework
#define DYNAMITE_DEFINED_BiLocusFramework
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_BiLocus(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [BiLocus *]
 *
 * Return [UNKN ]  Undocumented return value [BiLocus *]
 *
 */
BiLocus * Wise2_hard_link_BiLocus(BiLocus * obj);
#define hard_link_BiLocus Wise2_hard_link_BiLocus


/* Function:  BiLocus_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BiLocus *]
 *
 */
BiLocus * Wise2_BiLocus_alloc(void);
#define BiLocus_alloc Wise2_BiLocus_alloc


/* Function:  free_BiLocus(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [BiLocus *]
 *
 * Return [UNKN ]  Undocumented return value [BiLocus *]
 *
 */
BiLocus * Wise2_free_BiLocus(BiLocus * obj);
#define free_BiLocus Wise2_free_BiLocus


/* Function:  add_BiLocusFramework(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [BiLocusFramework *]
 * Arg:        add [OWNER] Object to add to the list [BiLocus *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_BiLocusFramework(BiLocusFramework * obj,BiLocus * add);
#define add_BiLocusFramework Wise2_add_BiLocusFramework


/* Function:  flush_BiLocusFramework(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [BiLocusFramework *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_BiLocusFramework(BiLocusFramework * obj);
#define flush_BiLocusFramework Wise2_flush_BiLocusFramework


/* Function:  BiLocusFramework_alloc_std(void)
 *
 * Descrip:    Equivalent to BiLocusFramework_alloc_len(BiLocusFrameworkLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BiLocusFramework *]
 *
 */
BiLocusFramework * Wise2_BiLocusFramework_alloc_std(void);
#define BiLocusFramework_alloc_std Wise2_BiLocusFramework_alloc_std


/* Function:  BiLocusFramework_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [BiLocusFramework *]
 *
 */
BiLocusFramework * Wise2_BiLocusFramework_alloc_len(int len);
#define BiLocusFramework_alloc_len Wise2_BiLocusFramework_alloc_len


/* Function:  hard_link_BiLocusFramework(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [BiLocusFramework *]
 *
 * Return [UNKN ]  Undocumented return value [BiLocusFramework *]
 *
 */
BiLocusFramework * Wise2_hard_link_BiLocusFramework(BiLocusFramework * obj);
#define hard_link_BiLocusFramework Wise2_hard_link_BiLocusFramework


/* Function:  BiLocusFramework_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BiLocusFramework *]
 *
 */
BiLocusFramework * Wise2_BiLocusFramework_alloc(void);
#define BiLocusFramework_alloc Wise2_BiLocusFramework_alloc


/* Function:  free_BiLocusFramework(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [BiLocusFramework *]
 *
 * Return [UNKN ]  Undocumented return value [BiLocusFramework *]
 *
 */
BiLocusFramework * Wise2_free_BiLocusFramework(BiLocusFramework * obj);
#define free_BiLocusFramework Wise2_free_BiLocusFramework


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
BiGenotypeType Wise2_type_from_hapmap_string_BiLocus(BiLocus * bi,char * g);
#define type_from_hapmap_string_BiLocus Wise2_type_from_hapmap_string_BiLocus
BiLocusFramework * Wise2_new_BiLocusFramework(void);
#define new_BiLocusFramework Wise2_new_BiLocusFramework
BiLocus * Wise2_find_or_new_BiLocus_from_BiLocusFramework(BiLocusFramework * bgf,char * ref_id,char first,char second);
#define find_or_new_BiLocus_from_BiLocusFramework Wise2_find_or_new_BiLocus_from_BiLocusFramework


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_BiLocusFramework(BiLocus ** list,int i,int j) ;
#define swap_BiLocusFramework Wise2_swap_BiLocusFramework
void Wise2_qsort_BiLocusFramework(BiLocus ** list,int left,int right,int (*comp)(BiLocus * ,BiLocus * ));
#define qsort_BiLocusFramework Wise2_qsort_BiLocusFramework
void Wise2_sort_BiLocusFramework(BiLocusFramework * obj,int (*comp)(BiLocus *, BiLocus *));
#define sort_BiLocusFramework Wise2_sort_BiLocusFramework
boolean Wise2_expand_BiLocusFramework(BiLocusFramework * obj,int len);
#define expand_BiLocusFramework Wise2_expand_BiLocusFramework

#ifdef _cplusplus
}
#endif

#endif
