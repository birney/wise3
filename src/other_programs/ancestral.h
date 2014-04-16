#ifndef DYNAMITEancestralHEADERFILE
#define DYNAMITEancestralHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "snplocus.h"
#include "plot.h"

#define AncestralVarChrLISTLENGTH 1048
#define AncestralVarSetLISTLENGTH 60

#define AncestralBlockStrainLISTLENGTH 60
#define AncestralBlockSetLISTLENGTH 60

struct Wise2_AncestralVarLocus {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    VarLocus * var;  
    char * ind; /*  short int of length individual, number of ancestor number */ 
    char * anc; /*  shor int of length ancestor number, 0 for reference, 1 for alternate */ 
    } ;  
/* AncestralVarLocus defined */ 
#ifndef DYNAMITE_DEFINED_AncestralVarLocus
typedef struct Wise2_AncestralVarLocus Wise2_AncestralVarLocus;
#define AncestralVarLocus Wise2_AncestralVarLocus
#define DYNAMITE_DEFINED_AncestralVarLocus
#endif


struct Wise2_AncestralVarChr {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * chr;  
    AncestralVarLocus ** loci;   
    int len;/* len for above loci  */ 
    int maxlen; /* maxlen for above loci */ 
    } ;  
/* AncestralVarChr defined */ 
#ifndef DYNAMITE_DEFINED_AncestralVarChr
typedef struct Wise2_AncestralVarChr Wise2_AncestralVarChr;
#define AncestralVarChr Wise2_AncestralVarChr
#define DYNAMITE_DEFINED_AncestralVarChr
#endif


struct Wise2_AncestralIndividual {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * name;     
    char single_letter;  
    } ;  
/* AncestralIndividual defined */ 
#ifndef DYNAMITE_DEFINED_AncestralIndividual
typedef struct Wise2_AncestralIndividual Wise2_AncestralIndividual;
#define AncestralIndividual Wise2_AncestralIndividual
#define DYNAMITE_DEFINED_AncestralIndividual
#endif


struct Wise2_AncestralVarSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    AncestralVarChr ** chr;  
    int len;/* len for above chr  */ 
    int maxlen; /* maxlen for above chr */ 
    Individual ** ind;   
    int ind_len;/* len for above ind  */ 
    int ind_maxlen; /* maxlen for above ind */ 
    AncestralIndividual ** anc;  
    int anc_len;/* len for above anc  */ 
    int anc_maxlen; /* maxlen for above anc */ 
    } ;  
/* AncestralVarSet defined */ 
#ifndef DYNAMITE_DEFINED_AncestralVarSet
typedef struct Wise2_AncestralVarSet Wise2_AncestralVarSet;
#define AncestralVarSet Wise2_AncestralVarSet
#define DYNAMITE_DEFINED_AncestralVarSet
#endif


struct Wise2_AncestralBlock {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * chr;  
    long start;  
    long end;    
    int snp_count;   
    int anc_diff;    
    AncestralIndividual * ai;    
    } ;  
/* AncestralBlock defined */ 
#ifndef DYNAMITE_DEFINED_AncestralBlock
typedef struct Wise2_AncestralBlock Wise2_AncestralBlock;
#define AncestralBlock Wise2_AncestralBlock
#define DYNAMITE_DEFINED_AncestralBlock
#endif


struct Wise2_AncestralBlockStrain {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    AncestralBlock ** block;     
    int len;/* len for above block  */ 
    int maxlen; /* maxlen for above block */ 
    Individual * ind;    
    } ;  
/* AncestralBlockStrain defined */ 
#ifndef DYNAMITE_DEFINED_AncestralBlockStrain
typedef struct Wise2_AncestralBlockStrain Wise2_AncestralBlockStrain;
#define AncestralBlockStrain Wise2_AncestralBlockStrain
#define DYNAMITE_DEFINED_AncestralBlockStrain
#endif


struct Wise2_AncestralBlockSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    AncestralBlockStrain ** strain;  
    int len;/* len for above strain  */ 
    int maxlen; /* maxlen for above strain */ 
    } ;  
/* AncestralBlockSet defined */ 
#ifndef DYNAMITE_DEFINED_AncestralBlockSet
typedef struct Wise2_AncestralBlockSet Wise2_AncestralBlockSet;
#define AncestralBlockSet Wise2_AncestralBlockSet
#define DYNAMITE_DEFINED_AncestralBlockSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_AncestralVarLocus(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AncestralVarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarLocus *]
 *
 */
AncestralVarLocus * Wise2_hard_link_AncestralVarLocus(AncestralVarLocus * obj);
#define hard_link_AncestralVarLocus Wise2_hard_link_AncestralVarLocus


/* Function:  AncestralVarLocus_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarLocus *]
 *
 */
AncestralVarLocus * Wise2_AncestralVarLocus_alloc(void);
#define AncestralVarLocus_alloc Wise2_AncestralVarLocus_alloc


/* Function:  free_AncestralVarLocus(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AncestralVarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarLocus *]
 *
 */
AncestralVarLocus * Wise2_free_AncestralVarLocus(AncestralVarLocus * obj);
#define free_AncestralVarLocus Wise2_free_AncestralVarLocus


/* Function:  add_AncestralVarChr(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralVarChr *]
 * Arg:        add [OWNER] Object to add to the list [AncestralVarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_AncestralVarChr(AncestralVarChr * obj,AncestralVarLocus * add);
#define add_AncestralVarChr Wise2_add_AncestralVarChr


/* Function:  flush_AncestralVarChr(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AncestralVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_AncestralVarChr(AncestralVarChr * obj);
#define flush_AncestralVarChr Wise2_flush_AncestralVarChr


/* Function:  AncestralVarChr_alloc_std(void)
 *
 * Descrip:    Equivalent to AncestralVarChr_alloc_len(AncestralVarChrLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarChr *]
 *
 */
AncestralVarChr * Wise2_AncestralVarChr_alloc_std(void);
#define AncestralVarChr_alloc_std Wise2_AncestralVarChr_alloc_std


/* Function:  AncestralVarChr_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarChr *]
 *
 */
AncestralVarChr * Wise2_AncestralVarChr_alloc_len(int len);
#define AncestralVarChr_alloc_len Wise2_AncestralVarChr_alloc_len


/* Function:  hard_link_AncestralVarChr(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AncestralVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarChr *]
 *
 */
AncestralVarChr * Wise2_hard_link_AncestralVarChr(AncestralVarChr * obj);
#define hard_link_AncestralVarChr Wise2_hard_link_AncestralVarChr


/* Function:  AncestralVarChr_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarChr *]
 *
 */
AncestralVarChr * Wise2_AncestralVarChr_alloc(void);
#define AncestralVarChr_alloc Wise2_AncestralVarChr_alloc


/* Function:  free_AncestralVarChr(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AncestralVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarChr *]
 *
 */
AncestralVarChr * Wise2_free_AncestralVarChr(AncestralVarChr * obj);
#define free_AncestralVarChr Wise2_free_AncestralVarChr


/* Function:  hard_link_AncestralIndividual(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AncestralIndividual *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralIndividual *]
 *
 */
AncestralIndividual * Wise2_hard_link_AncestralIndividual(AncestralIndividual * obj);
#define hard_link_AncestralIndividual Wise2_hard_link_AncestralIndividual


/* Function:  AncestralIndividual_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralIndividual *]
 *
 */
AncestralIndividual * Wise2_AncestralIndividual_alloc(void);
#define AncestralIndividual_alloc Wise2_AncestralIndividual_alloc


/* Function:  free_AncestralIndividual(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AncestralIndividual *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralIndividual *]
 *
 */
AncestralIndividual * Wise2_free_AncestralIndividual(AncestralIndividual * obj);
#define free_AncestralIndividual Wise2_free_AncestralIndividual


/* Function:  add_AncestralVarSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralVarSet *]
 * Arg:        add [OWNER] Object to add to the list [AncestralVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_AncestralVarSet(AncestralVarSet * obj,AncestralVarChr * add);
#define add_AncestralVarSet Wise2_add_AncestralVarSet


/* Function:  flush_AncestralVarSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AncestralVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_AncestralVarSet(AncestralVarSet * obj);
#define flush_AncestralVarSet Wise2_flush_AncestralVarSet


/* Function:  add_ind_AncestralVarSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralVarSet *]
 * Arg:        add [OWNER] Object to add to the list [Individual *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ind_AncestralVarSet(AncestralVarSet * obj,Individual * add);
#define add_ind_AncestralVarSet Wise2_add_ind_AncestralVarSet


/* Function:  flush_ind_AncestralVarSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AncestralVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_ind_AncestralVarSet(AncestralVarSet * obj);
#define flush_ind_AncestralVarSet Wise2_flush_ind_AncestralVarSet


/* Function:  add_anc_AncestralVarSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralVarSet *]
 * Arg:        add [OWNER] Object to add to the list [AncestralIndividual *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_anc_AncestralVarSet(AncestralVarSet * obj,AncestralIndividual * add);
#define add_anc_AncestralVarSet Wise2_add_anc_AncestralVarSet


/* Function:  flush_anc_AncestralVarSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AncestralVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_anc_AncestralVarSet(AncestralVarSet * obj);
#define flush_anc_AncestralVarSet Wise2_flush_anc_AncestralVarSet


/* Function:  AncestralVarSet_alloc_std(void)
 *
 * Descrip:    Equivalent to AncestralVarSet_alloc_len(AncestralVarSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarSet *]
 *
 */
AncestralVarSet * Wise2_AncestralVarSet_alloc_std(void);
#define AncestralVarSet_alloc_std Wise2_AncestralVarSet_alloc_std


/* Function:  AncestralVarSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarSet *]
 *
 */
AncestralVarSet * Wise2_AncestralVarSet_alloc_len(int len);
#define AncestralVarSet_alloc_len Wise2_AncestralVarSet_alloc_len


/* Function:  hard_link_AncestralVarSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AncestralVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarSet *]
 *
 */
AncestralVarSet * Wise2_hard_link_AncestralVarSet(AncestralVarSet * obj);
#define hard_link_AncestralVarSet Wise2_hard_link_AncestralVarSet


/* Function:  AncestralVarSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarSet *]
 *
 */
AncestralVarSet * Wise2_AncestralVarSet_alloc(void);
#define AncestralVarSet_alloc Wise2_AncestralVarSet_alloc


/* Function:  free_AncestralVarSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AncestralVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarSet *]
 *
 */
AncestralVarSet * Wise2_free_AncestralVarSet(AncestralVarSet * obj);
#define free_AncestralVarSet Wise2_free_AncestralVarSet


/* Function:  hard_link_AncestralBlock(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AncestralBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlock *]
 *
 */
AncestralBlock * Wise2_hard_link_AncestralBlock(AncestralBlock * obj);
#define hard_link_AncestralBlock Wise2_hard_link_AncestralBlock


/* Function:  AncestralBlock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlock *]
 *
 */
AncestralBlock * Wise2_AncestralBlock_alloc(void);
#define AncestralBlock_alloc Wise2_AncestralBlock_alloc


/* Function:  free_AncestralBlock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AncestralBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlock *]
 *
 */
AncestralBlock * Wise2_free_AncestralBlock(AncestralBlock * obj);
#define free_AncestralBlock Wise2_free_AncestralBlock


/* Function:  add_AncestralBlockStrain(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralBlockStrain *]
 * Arg:        add [OWNER] Object to add to the list [AncestralBlock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_AncestralBlockStrain(AncestralBlockStrain * obj,AncestralBlock * add);
#define add_AncestralBlockStrain Wise2_add_AncestralBlockStrain


/* Function:  flush_AncestralBlockStrain(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AncestralBlockStrain *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_AncestralBlockStrain(AncestralBlockStrain * obj);
#define flush_AncestralBlockStrain Wise2_flush_AncestralBlockStrain


/* Function:  AncestralBlockStrain_alloc_std(void)
 *
 * Descrip:    Equivalent to AncestralBlockStrain_alloc_len(AncestralBlockStrainLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockStrain *]
 *
 */
AncestralBlockStrain * Wise2_AncestralBlockStrain_alloc_std(void);
#define AncestralBlockStrain_alloc_std Wise2_AncestralBlockStrain_alloc_std


/* Function:  AncestralBlockStrain_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockStrain *]
 *
 */
AncestralBlockStrain * Wise2_AncestralBlockStrain_alloc_len(int len);
#define AncestralBlockStrain_alloc_len Wise2_AncestralBlockStrain_alloc_len


/* Function:  hard_link_AncestralBlockStrain(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AncestralBlockStrain *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockStrain *]
 *
 */
AncestralBlockStrain * Wise2_hard_link_AncestralBlockStrain(AncestralBlockStrain * obj);
#define hard_link_AncestralBlockStrain Wise2_hard_link_AncestralBlockStrain


/* Function:  AncestralBlockStrain_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockStrain *]
 *
 */
AncestralBlockStrain * Wise2_AncestralBlockStrain_alloc(void);
#define AncestralBlockStrain_alloc Wise2_AncestralBlockStrain_alloc


/* Function:  free_AncestralBlockStrain(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AncestralBlockStrain *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockStrain *]
 *
 */
AncestralBlockStrain * Wise2_free_AncestralBlockStrain(AncestralBlockStrain * obj);
#define free_AncestralBlockStrain Wise2_free_AncestralBlockStrain


/* Function:  add_AncestralBlockSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralBlockSet *]
 * Arg:        add [OWNER] Object to add to the list [AncestralBlockStrain *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_AncestralBlockSet(AncestralBlockSet * obj,AncestralBlockStrain * add);
#define add_AncestralBlockSet Wise2_add_AncestralBlockSet


/* Function:  flush_AncestralBlockSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AncestralBlockSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_AncestralBlockSet(AncestralBlockSet * obj);
#define flush_AncestralBlockSet Wise2_flush_AncestralBlockSet


/* Function:  AncestralBlockSet_alloc_std(void)
 *
 * Descrip:    Equivalent to AncestralBlockSet_alloc_len(AncestralBlockSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockSet *]
 *
 */
AncestralBlockSet * Wise2_AncestralBlockSet_alloc_std(void);
#define AncestralBlockSet_alloc_std Wise2_AncestralBlockSet_alloc_std


/* Function:  AncestralBlockSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockSet *]
 *
 */
AncestralBlockSet * Wise2_AncestralBlockSet_alloc_len(int len);
#define AncestralBlockSet_alloc_len Wise2_AncestralBlockSet_alloc_len


/* Function:  hard_link_AncestralBlockSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AncestralBlockSet *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockSet *]
 *
 */
AncestralBlockSet * Wise2_hard_link_AncestralBlockSet(AncestralBlockSet * obj);
#define hard_link_AncestralBlockSet Wise2_hard_link_AncestralBlockSet


/* Function:  AncestralBlockSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockSet *]
 *
 */
AncestralBlockSet * Wise2_AncestralBlockSet_alloc(void);
#define AncestralBlockSet_alloc Wise2_AncestralBlockSet_alloc


/* Function:  free_AncestralBlockSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AncestralBlockSet *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockSet *]
 *
 */
AncestralBlockSet * Wise2_free_AncestralBlockSet(AncestralBlockSet * obj);
#define free_AncestralBlockSet Wise2_free_AncestralBlockSet


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_dump_AncestralBlockSet(AncestralBlockSet * abs,FILE * ofp);
#define dump_AncestralBlockSet Wise2_dump_AncestralBlockSet
AncestralBlockSet * Wise2_AncestralBlockSet_from_AncestralVarSet(AncestralVarSet * avs,GenoVarSet * gvs);
#define AncestralBlockSet_from_AncestralVarSet Wise2_AncestralBlockSet_from_AncestralVarSet
void Wise2_dump_as_pairwise_AncestralVarSet(AncestralVarSet * avs,GenoVarSet * gvs,FILE * ofp);
#define dump_as_pairwise_AncestralVarSet Wise2_dump_as_pairwise_AncestralVarSet
void Wise2_write_simple_AncestralVarSet(AncestralVarSet * avs,FILE * ofp);
#define write_simple_AncestralVarSet Wise2_write_simple_AncestralVarSet
FrameSet * Wise2_AncestralVarChr_as_FrameSet(AncestralVarSet * avs,AncestralVarChr * avc);
#define AncestralVarChr_as_FrameSet Wise2_AncestralVarChr_as_FrameSet


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_AncestralVarChr(AncestralVarLocus ** list,int i,int j) ;
#define swap_AncestralVarChr Wise2_swap_AncestralVarChr
void Wise2_qsort_AncestralVarChr(AncestralVarLocus ** list,int left,int right,int (*comp)(AncestralVarLocus * ,AncestralVarLocus * ));
#define qsort_AncestralVarChr Wise2_qsort_AncestralVarChr
void Wise2_sort_AncestralVarChr(AncestralVarChr * obj,int (*comp)(AncestralVarLocus *, AncestralVarLocus *));
#define sort_AncestralVarChr Wise2_sort_AncestralVarChr
boolean Wise2_expand_AncestralVarChr(AncestralVarChr * obj,int len);
#define expand_AncestralVarChr Wise2_expand_AncestralVarChr
void Wise2_swap_AncestralVarSet(AncestralVarChr ** list,int i,int j) ;
#define swap_AncestralVarSet Wise2_swap_AncestralVarSet
void Wise2_qsort_AncestralVarSet(AncestralVarChr ** list,int left,int right,int (*comp)(AncestralVarChr * ,AncestralVarChr * ));
#define qsort_AncestralVarSet Wise2_qsort_AncestralVarSet
void Wise2_sort_AncestralVarSet(AncestralVarSet * obj,int (*comp)(AncestralVarChr *, AncestralVarChr *));
#define sort_AncestralVarSet Wise2_sort_AncestralVarSet
boolean Wise2_expand_AncestralVarSet(AncestralVarSet * obj,int len);
#define expand_AncestralVarSet Wise2_expand_AncestralVarSet
void Wise2_swap_ind_AncestralVarSet(Individual ** list,int i,int j) ;
#define swap_ind_AncestralVarSet Wise2_swap_ind_AncestralVarSet
void Wise2_qsort_ind_AncestralVarSet(Individual ** list,int left,int right,int (*comp)(Individual * ,Individual * ));
#define qsort_ind_AncestralVarSet Wise2_qsort_ind_AncestralVarSet
void Wise2_sort_ind_AncestralVarSet(AncestralVarSet * obj,int (*comp)(Individual *, Individual *));
#define sort_ind_AncestralVarSet Wise2_sort_ind_AncestralVarSet
boolean Wise2_expand_ind_AncestralVarSet(AncestralVarSet * obj,int len);
#define expand_ind_AncestralVarSet Wise2_expand_ind_AncestralVarSet
void Wise2_swap_anc_AncestralVarSet(AncestralIndividual ** list,int i,int j) ;
#define swap_anc_AncestralVarSet Wise2_swap_anc_AncestralVarSet
void Wise2_qsort_anc_AncestralVarSet(AncestralIndividual ** list,int left,int right,int (*comp)(AncestralIndividual * ,AncestralIndividual * ));
#define qsort_anc_AncestralVarSet Wise2_qsort_anc_AncestralVarSet
void Wise2_sort_anc_AncestralVarSet(AncestralVarSet * obj,int (*comp)(AncestralIndividual *, AncestralIndividual *));
#define sort_anc_AncestralVarSet Wise2_sort_anc_AncestralVarSet
boolean Wise2_expand_anc_AncestralVarSet(AncestralVarSet * obj,int len);
#define expand_anc_AncestralVarSet Wise2_expand_anc_AncestralVarSet
void Wise2_swap_AncestralBlockStrain(AncestralBlock ** list,int i,int j) ;
#define swap_AncestralBlockStrain Wise2_swap_AncestralBlockStrain
void Wise2_qsort_AncestralBlockStrain(AncestralBlock ** list,int left,int right,int (*comp)(AncestralBlock * ,AncestralBlock * ));
#define qsort_AncestralBlockStrain Wise2_qsort_AncestralBlockStrain
void Wise2_sort_AncestralBlockStrain(AncestralBlockStrain * obj,int (*comp)(AncestralBlock *, AncestralBlock *));
#define sort_AncestralBlockStrain Wise2_sort_AncestralBlockStrain
boolean Wise2_expand_AncestralBlockStrain(AncestralBlockStrain * obj,int len);
#define expand_AncestralBlockStrain Wise2_expand_AncestralBlockStrain
void Wise2_swap_AncestralBlockSet(AncestralBlockStrain ** list,int i,int j) ;
#define swap_AncestralBlockSet Wise2_swap_AncestralBlockSet
void Wise2_qsort_AncestralBlockSet(AncestralBlockStrain ** list,int left,int right,int (*comp)(AncestralBlockStrain * ,AncestralBlockStrain * ));
#define qsort_AncestralBlockSet Wise2_qsort_AncestralBlockSet
void Wise2_sort_AncestralBlockSet(AncestralBlockSet * obj,int (*comp)(AncestralBlockStrain *, AncestralBlockStrain *));
#define sort_AncestralBlockSet Wise2_sort_AncestralBlockSet
boolean Wise2_expand_AncestralBlockSet(AncestralBlockSet * obj,int len);
#define expand_AncestralBlockSet Wise2_expand_AncestralBlockSet

#ifdef _cplusplus
}
#endif

#endif
