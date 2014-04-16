#ifndef DYNAMITEpersonHEADERFILE
#define DYNAMITEpersonHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

#define PopulationLISTLENGTH 128
#define PopulationSetLISTLENGTH 128

#ifndef DYNAMITE_DEFINED_Population
typedef struct Wise2_Population Wise2_Population;
#define Population Wise2_Population
#define DYNAMITE_DEFINED_Population
#endif

struct Wise2_Person {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * person_id;    
    Population * pop;    
    } ;  
/* Person defined */ 
#ifndef DYNAMITE_DEFINED_Person
typedef struct Wise2_Person Wise2_Person;
#define Person Wise2_Person
#define DYNAMITE_DEFINED_Person
#endif


struct Wise2_Population {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * pop_name;     
    Person ** person;    
    int len;/* len for above person  */ 
    int maxlen; /* maxlen for above person */ 
    } ;  
/* Population defined */ 
#ifndef DYNAMITE_DEFINED_Population
typedef struct Wise2_Population Wise2_Population;
#define Population Wise2_Population
#define DYNAMITE_DEFINED_Population
#endif


struct Wise2_PopulationSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Population ** pop;   
    int len;/* len for above pop  */ 
    int maxlen; /* maxlen for above pop */ 
    } ;  
/* PopulationSet defined */ 
#ifndef DYNAMITE_DEFINED_PopulationSet
typedef struct Wise2_PopulationSet Wise2_PopulationSet;
#define PopulationSet Wise2_PopulationSet
#define DYNAMITE_DEFINED_PopulationSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_Person(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Person *]
 *
 * Return [UNKN ]  Undocumented return value [Person *]
 *
 */
Person * Wise2_hard_link_Person(Person * obj);
#define hard_link_Person Wise2_hard_link_Person


/* Function:  Person_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Person *]
 *
 */
Person * Wise2_Person_alloc(void);
#define Person_alloc Wise2_Person_alloc


/* Function:  free_Person(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Person *]
 *
 * Return [UNKN ]  Undocumented return value [Person *]
 *
 */
Person * Wise2_free_Person(Person * obj);
#define free_Person Wise2_free_Person


/* Function:  add_Population(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Population *]
 * Arg:        add [OWNER] Object to add to the list [Person *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_Population(Population * obj,Person * add);
#define add_Population Wise2_add_Population


/* Function:  flush_Population(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Population *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_Population(Population * obj);
#define flush_Population Wise2_flush_Population


/* Function:  Population_alloc_std(void)
 *
 * Descrip:    Equivalent to Population_alloc_len(PopulationLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Population *]
 *
 */
Population * Wise2_Population_alloc_std(void);
#define Population_alloc_std Wise2_Population_alloc_std


/* Function:  Population_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Population *]
 *
 */
Population * Wise2_Population_alloc_len(int len);
#define Population_alloc_len Wise2_Population_alloc_len


/* Function:  hard_link_Population(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Population *]
 *
 * Return [UNKN ]  Undocumented return value [Population *]
 *
 */
Population * Wise2_hard_link_Population(Population * obj);
#define hard_link_Population Wise2_hard_link_Population


/* Function:  Population_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Population *]
 *
 */
Population * Wise2_Population_alloc(void);
#define Population_alloc Wise2_Population_alloc


/* Function:  free_Population(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Population *]
 *
 * Return [UNKN ]  Undocumented return value [Population *]
 *
 */
Population * Wise2_free_Population(Population * obj);
#define free_Population Wise2_free_Population


/* Function:  add_PopulationSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PopulationSet *]
 * Arg:        add [OWNER] Object to add to the list [Population *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_PopulationSet(PopulationSet * obj,Population * add);
#define add_PopulationSet Wise2_add_PopulationSet


/* Function:  flush_PopulationSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [PopulationSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_PopulationSet(PopulationSet * obj);
#define flush_PopulationSet Wise2_flush_PopulationSet


/* Function:  PopulationSet_alloc_std(void)
 *
 * Descrip:    Equivalent to PopulationSet_alloc_len(PopulationSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PopulationSet *]
 *
 */
PopulationSet * Wise2_PopulationSet_alloc_std(void);
#define PopulationSet_alloc_std Wise2_PopulationSet_alloc_std


/* Function:  PopulationSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [PopulationSet *]
 *
 */
PopulationSet * Wise2_PopulationSet_alloc_len(int len);
#define PopulationSet_alloc_len Wise2_PopulationSet_alloc_len


/* Function:  hard_link_PopulationSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PopulationSet *]
 *
 * Return [UNKN ]  Undocumented return value [PopulationSet *]
 *
 */
PopulationSet * Wise2_hard_link_PopulationSet(PopulationSet * obj);
#define hard_link_PopulationSet Wise2_hard_link_PopulationSet


/* Function:  PopulationSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PopulationSet *]
 *
 */
PopulationSet * Wise2_PopulationSet_alloc(void);
#define PopulationSet_alloc Wise2_PopulationSet_alloc


/* Function:  free_PopulationSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PopulationSet *]
 *
 * Return [UNKN ]  Undocumented return value [PopulationSet *]
 *
 */
PopulationSet * Wise2_free_PopulationSet(PopulationSet * obj);
#define free_PopulationSet Wise2_free_PopulationSet


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
boolean Wise2_is_in_Population(Population * p,char * person_id);
#define is_in_Population Wise2_is_in_Population
Population * Wise2_find_or_new_Population_in_PopulationSet(PopulationSet * ps,char * pop_name);
#define find_or_new_Population_in_PopulationSet Wise2_find_or_new_Population_in_PopulationSet
Person * Wise2_find_or_new_Person_in_Population(Population * p,char * person_id);
#define find_or_new_Person_in_Population Wise2_find_or_new_Person_in_Population


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_Population(Person ** list,int i,int j) ;
#define swap_Population Wise2_swap_Population
void Wise2_qsort_Population(Person ** list,int left,int right,int (*comp)(Person * ,Person * ));
#define qsort_Population Wise2_qsort_Population
void Wise2_sort_Population(Population * obj,int (*comp)(Person *, Person *));
#define sort_Population Wise2_sort_Population
boolean Wise2_expand_Population(Population * obj,int len);
#define expand_Population Wise2_expand_Population
void Wise2_swap_PopulationSet(Population ** list,int i,int j) ;
#define swap_PopulationSet Wise2_swap_PopulationSet
void Wise2_qsort_PopulationSet(Population ** list,int left,int right,int (*comp)(Population * ,Population * ));
#define qsort_PopulationSet Wise2_qsort_PopulationSet
void Wise2_sort_PopulationSet(PopulationSet * obj,int (*comp)(Population *, Population *));
#define sort_PopulationSet Wise2_sort_PopulationSet
boolean Wise2_expand_PopulationSet(PopulationSet * obj,int len);
#define expand_PopulationSet Wise2_expand_PopulationSet

#ifdef _cplusplus
}
#endif

#endif
