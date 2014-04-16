#ifndef DYNAMITEsnplocusHEADERFILE
#define DYNAMITEsnplocusHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"
#include "wisestring.h"

#define GenoVarSetLISTLENGTH 256
#define GenoVarChrLISTLENGTH 256

#define HOMO_REF 0
#define SIMPLE_HET 1
#define HOMO_ALT 2
#define GENO_NOT_DEFINED 100

#define SIMPLE_SNP_LOCUS 56
#define COMPLEX_SNP_LOCUS 58

#ifndef DYNAMITE_DEFINED_GenoVarChr
typedef struct Wise2_GenoVarChr Wise2_GenoVarChr;
#define GenoVarChr Wise2_GenoVarChr
#define DYNAMITE_DEFINED_GenoVarChr
#endif

struct Wise2_VarLocus {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GenoVarChr * chr;    
    long pos;    
    char * ref_allele;   
    char * alt_allele;   
    int locus_type;  
    } ;  
/* VarLocus defined */ 
#ifndef DYNAMITE_DEFINED_VarLocus
typedef struct Wise2_VarLocus Wise2_VarLocus;
#define VarLocus Wise2_VarLocus
#define DYNAMITE_DEFINED_VarLocus
#endif


struct Wise2_GenoVarLocus {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    VarLocus * var;  
    char * ind;  
    } ;  
/* GenoVarLocus defined */ 
#ifndef DYNAMITE_DEFINED_GenoVarLocus
typedef struct Wise2_GenoVarLocus Wise2_GenoVarLocus;
#define GenoVarLocus Wise2_GenoVarLocus
#define DYNAMITE_DEFINED_GenoVarLocus
#endif


struct Wise2_Individual {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * id;   
    } ;  
/* Individual defined */ 
#ifndef DYNAMITE_DEFINED_Individual
typedef struct Wise2_Individual Wise2_Individual;
#define Individual Wise2_Individual
#define DYNAMITE_DEFINED_Individual
#endif


struct Wise2_GenoVarChr {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * chr;  
    GenoVarLocus ** loci;    
    int len;/* len for above loci  */ 
    int maxlen; /* maxlen for above loci */ 
    } ;  
/* GenoVarChr defined */ 
#ifndef DYNAMITE_DEFINED_GenoVarChr
typedef struct Wise2_GenoVarChr Wise2_GenoVarChr;
#define GenoVarChr Wise2_GenoVarChr
#define DYNAMITE_DEFINED_GenoVarChr
#endif


struct Wise2_GenoVarSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GenoVarChr ** chr;   
    int len;/* len for above chr  */ 
    int maxlen; /* maxlen for above chr */ 
    Individual ** ind;   
    int ind_len;/* len for above ind  */ 
    int ind_maxlen; /* maxlen for above ind */ 
    } ;  
/* GenoVarSet defined */ 
#ifndef DYNAMITE_DEFINED_GenoVarSet
typedef struct Wise2_GenoVarSet Wise2_GenoVarSet;
#define GenoVarSet Wise2_GenoVarSet
#define DYNAMITE_DEFINED_GenoVarSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_VarLocus(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [VarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [VarLocus *]
 *
 */
VarLocus * Wise2_hard_link_VarLocus(VarLocus * obj);
#define hard_link_VarLocus Wise2_hard_link_VarLocus


/* Function:  VarLocus_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [VarLocus *]
 *
 */
VarLocus * Wise2_VarLocus_alloc(void);
#define VarLocus_alloc Wise2_VarLocus_alloc


/* Function:  free_VarLocus(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [VarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [VarLocus *]
 *
 */
VarLocus * Wise2_free_VarLocus(VarLocus * obj);
#define free_VarLocus Wise2_free_VarLocus


/* Function:  hard_link_GenoVarLocus(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenoVarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarLocus *]
 *
 */
GenoVarLocus * Wise2_hard_link_GenoVarLocus(GenoVarLocus * obj);
#define hard_link_GenoVarLocus Wise2_hard_link_GenoVarLocus


/* Function:  GenoVarLocus_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenoVarLocus *]
 *
 */
GenoVarLocus * Wise2_GenoVarLocus_alloc(void);
#define GenoVarLocus_alloc Wise2_GenoVarLocus_alloc


/* Function:  free_GenoVarLocus(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenoVarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarLocus *]
 *
 */
GenoVarLocus * Wise2_free_GenoVarLocus(GenoVarLocus * obj);
#define free_GenoVarLocus Wise2_free_GenoVarLocus


/* Function:  hard_link_Individual(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Individual *]
 *
 * Return [UNKN ]  Undocumented return value [Individual *]
 *
 */
Individual * Wise2_hard_link_Individual(Individual * obj);
#define hard_link_Individual Wise2_hard_link_Individual


/* Function:  Individual_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Individual *]
 *
 */
Individual * Wise2_Individual_alloc(void);
#define Individual_alloc Wise2_Individual_alloc


/* Function:  free_Individual(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Individual *]
 *
 * Return [UNKN ]  Undocumented return value [Individual *]
 *
 */
Individual * Wise2_free_Individual(Individual * obj);
#define free_Individual Wise2_free_Individual


/* Function:  add_GenoVarChr(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenoVarChr *]
 * Arg:        add [OWNER] Object to add to the list [GenoVarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GenoVarChr(GenoVarChr * obj,GenoVarLocus * add);
#define add_GenoVarChr Wise2_add_GenoVarChr


/* Function:  flush_GenoVarChr(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenoVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_GenoVarChr(GenoVarChr * obj);
#define flush_GenoVarChr Wise2_flush_GenoVarChr


/* Function:  GenoVarChr_alloc_std(void)
 *
 * Descrip:    Equivalent to GenoVarChr_alloc_len(GenoVarChrLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenoVarChr *]
 *
 */
GenoVarChr * Wise2_GenoVarChr_alloc_std(void);
#define GenoVarChr_alloc_std Wise2_GenoVarChr_alloc_std


/* Function:  GenoVarChr_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarChr *]
 *
 */
GenoVarChr * Wise2_GenoVarChr_alloc_len(int len);
#define GenoVarChr_alloc_len Wise2_GenoVarChr_alloc_len


/* Function:  hard_link_GenoVarChr(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenoVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarChr *]
 *
 */
GenoVarChr * Wise2_hard_link_GenoVarChr(GenoVarChr * obj);
#define hard_link_GenoVarChr Wise2_hard_link_GenoVarChr


/* Function:  GenoVarChr_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenoVarChr *]
 *
 */
GenoVarChr * Wise2_GenoVarChr_alloc(void);
#define GenoVarChr_alloc Wise2_GenoVarChr_alloc


/* Function:  free_GenoVarChr(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenoVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarChr *]
 *
 */
GenoVarChr * Wise2_free_GenoVarChr(GenoVarChr * obj);
#define free_GenoVarChr Wise2_free_GenoVarChr


/* Function:  add_GenoVarSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenoVarSet *]
 * Arg:        add [OWNER] Object to add to the list [GenoVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GenoVarSet(GenoVarSet * obj,GenoVarChr * add);
#define add_GenoVarSet Wise2_add_GenoVarSet


/* Function:  flush_GenoVarSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenoVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_GenoVarSet(GenoVarSet * obj);
#define flush_GenoVarSet Wise2_flush_GenoVarSet


/* Function:  add_ind_GenoVarSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenoVarSet *]
 * Arg:        add [OWNER] Object to add to the list [Individual *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ind_GenoVarSet(GenoVarSet * obj,Individual * add);
#define add_ind_GenoVarSet Wise2_add_ind_GenoVarSet


/* Function:  flush_ind_GenoVarSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenoVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_ind_GenoVarSet(GenoVarSet * obj);
#define flush_ind_GenoVarSet Wise2_flush_ind_GenoVarSet


/* Function:  GenoVarSet_alloc_std(void)
 *
 * Descrip:    Equivalent to GenoVarSet_alloc_len(GenoVarSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenoVarSet *]
 *
 */
GenoVarSet * Wise2_GenoVarSet_alloc_std(void);
#define GenoVarSet_alloc_std Wise2_GenoVarSet_alloc_std


/* Function:  GenoVarSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarSet *]
 *
 */
GenoVarSet * Wise2_GenoVarSet_alloc_len(int len);
#define GenoVarSet_alloc_len Wise2_GenoVarSet_alloc_len


/* Function:  hard_link_GenoVarSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenoVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarSet *]
 *
 */
GenoVarSet * Wise2_hard_link_GenoVarSet(GenoVarSet * obj);
#define hard_link_GenoVarSet Wise2_hard_link_GenoVarSet


/* Function:  GenoVarSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenoVarSet *]
 *
 */
GenoVarSet * Wise2_GenoVarSet_alloc(void);
#define GenoVarSet_alloc Wise2_GenoVarSet_alloc


/* Function:  free_GenoVarSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenoVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarSet *]
 *
 */
GenoVarSet * Wise2_free_GenoVarSet(GenoVarSet * obj);
#define free_GenoVarSet Wise2_free_GenoVarSet


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
GenoVarSet * Wise2_only_simple_snp_loci_GenoVarSet(GenoVarSet * gvs);
#define only_simple_snp_loci_GenoVarSet Wise2_only_simple_snp_loci_GenoVarSet
int Wise2_individual_index_from_string_GenoVarSet(GenoVarSet * s, char * strain);
#define individual_index_from_string_GenoVarSet Wise2_individual_index_from_string_GenoVarSet
int Wise2_number_of_simple_snp_loci_GenoVarChr(GenoVarChr * c);
#define number_of_simple_snp_loci_GenoVarChr Wise2_number_of_simple_snp_loci_GenoVarChr
void Wise2_write_sanger_GenoVarSet(GenoVarSet * s,FILE * ofp);
#define write_sanger_GenoVarSet Wise2_write_sanger_GenoVarSet
GenoVarSet * Wise2_read_sanger_genotype_file(FILE * ifp);
#define read_sanger_genotype_file Wise2_read_sanger_genotype_file


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_GenoVarChr(GenoVarLocus ** list,int i,int j) ;
#define swap_GenoVarChr Wise2_swap_GenoVarChr
void Wise2_qsort_GenoVarChr(GenoVarLocus ** list,int left,int right,int (*comp)(GenoVarLocus * ,GenoVarLocus * ));
#define qsort_GenoVarChr Wise2_qsort_GenoVarChr
void Wise2_sort_GenoVarChr(GenoVarChr * obj,int (*comp)(GenoVarLocus *, GenoVarLocus *));
#define sort_GenoVarChr Wise2_sort_GenoVarChr
boolean Wise2_expand_GenoVarChr(GenoVarChr * obj,int len);
#define expand_GenoVarChr Wise2_expand_GenoVarChr
void Wise2_swap_GenoVarSet(GenoVarChr ** list,int i,int j) ;
#define swap_GenoVarSet Wise2_swap_GenoVarSet
void Wise2_qsort_GenoVarSet(GenoVarChr ** list,int left,int right,int (*comp)(GenoVarChr * ,GenoVarChr * ));
#define qsort_GenoVarSet Wise2_qsort_GenoVarSet
void Wise2_sort_GenoVarSet(GenoVarSet * obj,int (*comp)(GenoVarChr *, GenoVarChr *));
#define sort_GenoVarSet Wise2_sort_GenoVarSet
boolean Wise2_expand_GenoVarSet(GenoVarSet * obj,int len);
#define expand_GenoVarSet Wise2_expand_GenoVarSet
void Wise2_swap_ind_GenoVarSet(Individual ** list,int i,int j) ;
#define swap_ind_GenoVarSet Wise2_swap_ind_GenoVarSet
void Wise2_qsort_ind_GenoVarSet(Individual ** list,int left,int right,int (*comp)(Individual * ,Individual * ));
#define qsort_ind_GenoVarSet Wise2_qsort_ind_GenoVarSet
void Wise2_sort_ind_GenoVarSet(GenoVarSet * obj,int (*comp)(Individual *, Individual *));
#define sort_ind_GenoVarSet Wise2_sort_ind_GenoVarSet
boolean Wise2_expand_ind_GenoVarSet(GenoVarSet * obj,int len);
#define expand_ind_GenoVarSet Wise2_expand_ind_GenoVarSet

#ifdef _cplusplus
}
#endif

#endif
