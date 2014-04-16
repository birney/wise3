#ifndef DYNAMITEgraphseqHEADERFILE
#define DYNAMITEgraphseqHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

#define GraphSeqEdgeLISTLENGTH 32

#define GraphSeqLISTLENGTH 64

#ifndef DYNAMITE_DEFINED_GraphSeqEdge
typedef struct Wise2_GraphSeqEdge Wise2_GraphSeqEdge;
#define GraphSeqEdge Wise2_GraphSeqEdge
#define DYNAMITE_DEFINED_GraphSeqEdge
#endif

struct Wise2_GraphSeqEdgeHolder {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GraphSeqEdge * edge;     
    } ;  
/* GraphSeqEdgeHolder defined */ 
#ifndef DYNAMITE_DEFINED_GraphSeqEdgeHolder
typedef struct Wise2_GraphSeqEdgeHolder Wise2_GraphSeqEdgeHolder;
#define GraphSeqEdgeHolder Wise2_GraphSeqEdgeHolder
#define DYNAMITE_DEFINED_GraphSeqEdgeHolder
#endif


struct Wise2_GraphSeqEdge {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * id;   
    int len;     
    Sequence * seq;  
    GraphSeqEdgeHolder ** fivep;     
    int five_len;   /* len for above fivep  */ 
    int five_maxlen;/* maxlen for above fivep */ 
    GraphSeqEdgeHolder ** threep;    
    int three_len;  /* len for above threep  */ 
    int three_maxlen;   /* maxlen for above threep */ 
    } ;  
/* GraphSeqEdge defined */ 
#ifndef DYNAMITE_DEFINED_GraphSeqEdge
typedef struct Wise2_GraphSeqEdge Wise2_GraphSeqEdge;
#define GraphSeqEdge Wise2_GraphSeqEdge
#define DYNAMITE_DEFINED_GraphSeqEdge
#endif


/* Object GraphSeq
 *
 * Descrip: This represents a graph, or rather a series
 *        of independent graphs potentially. 
 *
 *        One drawback is that this system assumes there are free ends.
 *        This is bound to bite me at some point
 *
 *
 */
struct Wise2_GraphSeq {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    GraphSeqEdgeHolder ** free_ends;     
    int free_len;   /* len for above free_ends  */ 
    int free_maxlen;/* maxlen for above free_ends */ 
    GraphSeqEdge ** comp;    
    int len;/* len for above comp  */ 
    int maxlen; /* maxlen for above comp */ 
    } ;  
/* GraphSeq defined */ 
#ifndef DYNAMITE_DEFINED_GraphSeq
typedef struct Wise2_GraphSeq Wise2_GraphSeq;
#define GraphSeq Wise2_GraphSeq
#define DYNAMITE_DEFINED_GraphSeq
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_GraphSeqEdgeHolder(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GraphSeqEdgeHolder *]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdgeHolder *]
 *
 */
GraphSeqEdgeHolder * Wise2_hard_link_GraphSeqEdgeHolder(GraphSeqEdgeHolder * obj);
#define hard_link_GraphSeqEdgeHolder Wise2_hard_link_GraphSeqEdgeHolder


/* Function:  GraphSeqEdgeHolder_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdgeHolder *]
 *
 */
GraphSeqEdgeHolder * Wise2_GraphSeqEdgeHolder_alloc(void);
#define GraphSeqEdgeHolder_alloc Wise2_GraphSeqEdgeHolder_alloc


/* Function:  free_GraphSeqEdgeHolder(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GraphSeqEdgeHolder *]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdgeHolder *]
 *
 */
GraphSeqEdgeHolder * Wise2_free_GraphSeqEdgeHolder(GraphSeqEdgeHolder * obj);
#define free_GraphSeqEdgeHolder Wise2_free_GraphSeqEdgeHolder


/* Function:  add_five_GraphSeqEdge(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GraphSeqEdge *]
 * Arg:        add [OWNER] Object to add to the list [GraphSeqEdgeHolder *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_five_GraphSeqEdge(GraphSeqEdge * obj,GraphSeqEdgeHolder * add);
#define add_five_GraphSeqEdge Wise2_add_five_GraphSeqEdge


/* Function:  flush_five_GraphSeqEdge(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GraphSeqEdge *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_five_GraphSeqEdge(GraphSeqEdge * obj);
#define flush_five_GraphSeqEdge Wise2_flush_five_GraphSeqEdge


/* Function:  add_three_GraphSeqEdge(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GraphSeqEdge *]
 * Arg:        add [OWNER] Object to add to the list [GraphSeqEdgeHolder *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_three_GraphSeqEdge(GraphSeqEdge * obj,GraphSeqEdgeHolder * add);
#define add_three_GraphSeqEdge Wise2_add_three_GraphSeqEdge


/* Function:  flush_three_GraphSeqEdge(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GraphSeqEdge *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_three_GraphSeqEdge(GraphSeqEdge * obj);
#define flush_three_GraphSeqEdge Wise2_flush_three_GraphSeqEdge


/* Function:  GraphSeqEdge_alloc_std(void)
 *
 * Descrip:    Equivalent to GraphSeqEdge_alloc_len(GraphSeqEdgeLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdge *]
 *
 */
GraphSeqEdge * Wise2_GraphSeqEdge_alloc_std(void);
#define GraphSeqEdge_alloc_std Wise2_GraphSeqEdge_alloc_std


/* Function:  GraphSeqEdge_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdge *]
 *
 */
GraphSeqEdge * Wise2_GraphSeqEdge_alloc_len(int len);
#define GraphSeqEdge_alloc_len Wise2_GraphSeqEdge_alloc_len


/* Function:  hard_link_GraphSeqEdge(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GraphSeqEdge *]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdge *]
 *
 */
GraphSeqEdge * Wise2_hard_link_GraphSeqEdge(GraphSeqEdge * obj);
#define hard_link_GraphSeqEdge Wise2_hard_link_GraphSeqEdge


/* Function:  GraphSeqEdge_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdge *]
 *
 */
GraphSeqEdge * Wise2_GraphSeqEdge_alloc(void);
#define GraphSeqEdge_alloc Wise2_GraphSeqEdge_alloc


/* Function:  free_GraphSeqEdge(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GraphSeqEdge *]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdge *]
 *
 */
GraphSeqEdge * Wise2_free_GraphSeqEdge(GraphSeqEdge * obj);
#define free_GraphSeqEdge Wise2_free_GraphSeqEdge


/* Function:  add_free_GraphSeq(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GraphSeq *]
 * Arg:        add [OWNER] Object to add to the list [GraphSeqEdgeHolder *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_free_GraphSeq(GraphSeq * obj,GraphSeqEdgeHolder * add);
#define add_free_GraphSeq Wise2_add_free_GraphSeq


/* Function:  flush_free_GraphSeq(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GraphSeq *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_free_GraphSeq(GraphSeq * obj);
#define flush_free_GraphSeq Wise2_flush_free_GraphSeq


/* Function:  add_GraphSeq(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GraphSeq *]
 * Arg:        add [OWNER] Object to add to the list [GraphSeqEdge *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_GraphSeq(GraphSeq * obj,GraphSeqEdge * add);
#define add_GraphSeq Wise2_add_GraphSeq


/* Function:  flush_GraphSeq(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GraphSeq *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_GraphSeq(GraphSeq * obj);
#define flush_GraphSeq Wise2_flush_GraphSeq


/* Function:  GraphSeq_alloc_std(void)
 *
 * Descrip:    Equivalent to GraphSeq_alloc_len(GraphSeqLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GraphSeq *]
 *
 */
GraphSeq * Wise2_GraphSeq_alloc_std(void);
#define GraphSeq_alloc_std Wise2_GraphSeq_alloc_std


/* Function:  GraphSeq_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeq *]
 *
 */
GraphSeq * Wise2_GraphSeq_alloc_len(int len);
#define GraphSeq_alloc_len Wise2_GraphSeq_alloc_len


/* Function:  hard_link_GraphSeq(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GraphSeq *]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeq *]
 *
 */
GraphSeq * Wise2_hard_link_GraphSeq(GraphSeq * obj);
#define hard_link_GraphSeq Wise2_hard_link_GraphSeq


/* Function:  GraphSeq_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GraphSeq *]
 *
 */
GraphSeq * Wise2_GraphSeq_alloc(void);
#define GraphSeq_alloc Wise2_GraphSeq_alloc


/* Function:  free_GraphSeq(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GraphSeq *]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeq *]
 *
 */
GraphSeq * Wise2_free_GraphSeq(GraphSeq * obj);
#define free_GraphSeq Wise2_free_GraphSeq


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
boolean Wise2_prepare_GraphSeq(GraphSeq * gs);
#define prepare_GraphSeq Wise2_prepare_GraphSeq
GraphSeq * Wise2_read_simple_GraphSeq(FILE * ifp);
#define read_simple_GraphSeq Wise2_read_simple_GraphSeq
GraphSeqEdge * Wise2_find_GraphSeqEdge(GraphSeq * gs,char * id) ;
#define find_GraphSeqEdge Wise2_find_GraphSeqEdge


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_five_GraphSeqEdge(GraphSeqEdgeHolder ** list,int i,int j) ;
#define swap_five_GraphSeqEdge Wise2_swap_five_GraphSeqEdge
void Wise2_qsort_five_GraphSeqEdge(GraphSeqEdgeHolder ** list,int left,int right,int (*comp)(GraphSeqEdgeHolder * ,GraphSeqEdgeHolder * ));
#define qsort_five_GraphSeqEdge Wise2_qsort_five_GraphSeqEdge
void Wise2_sort_five_GraphSeqEdge(GraphSeqEdge * obj,int (*comp)(GraphSeqEdgeHolder *, GraphSeqEdgeHolder *));
#define sort_five_GraphSeqEdge Wise2_sort_five_GraphSeqEdge
boolean Wise2_expand_five_GraphSeqEdge(GraphSeqEdge * obj,int len);
#define expand_five_GraphSeqEdge Wise2_expand_five_GraphSeqEdge
void Wise2_swap_three_GraphSeqEdge(GraphSeqEdgeHolder ** list,int i,int j) ;
#define swap_three_GraphSeqEdge Wise2_swap_three_GraphSeqEdge
void Wise2_qsort_three_GraphSeqEdge(GraphSeqEdgeHolder ** list,int left,int right,int (*comp)(GraphSeqEdgeHolder * ,GraphSeqEdgeHolder * ));
#define qsort_three_GraphSeqEdge Wise2_qsort_three_GraphSeqEdge
void Wise2_sort_three_GraphSeqEdge(GraphSeqEdge * obj,int (*comp)(GraphSeqEdgeHolder *, GraphSeqEdgeHolder *));
#define sort_three_GraphSeqEdge Wise2_sort_three_GraphSeqEdge
boolean Wise2_expand_three_GraphSeqEdge(GraphSeqEdge * obj,int len);
#define expand_three_GraphSeqEdge Wise2_expand_three_GraphSeqEdge
void Wise2_swap_free_GraphSeq(GraphSeqEdgeHolder ** list,int i,int j) ;
#define swap_free_GraphSeq Wise2_swap_free_GraphSeq
void Wise2_qsort_free_GraphSeq(GraphSeqEdgeHolder ** list,int left,int right,int (*comp)(GraphSeqEdgeHolder * ,GraphSeqEdgeHolder * ));
#define qsort_free_GraphSeq Wise2_qsort_free_GraphSeq
void Wise2_sort_free_GraphSeq(GraphSeq * obj,int (*comp)(GraphSeqEdgeHolder *, GraphSeqEdgeHolder *));
#define sort_free_GraphSeq Wise2_sort_free_GraphSeq
boolean Wise2_expand_free_GraphSeq(GraphSeq * obj,int len);
#define expand_free_GraphSeq Wise2_expand_free_GraphSeq
void Wise2_swap_GraphSeq(GraphSeqEdge ** list,int i,int j) ;
#define swap_GraphSeq Wise2_swap_GraphSeq
void Wise2_qsort_GraphSeq(GraphSeqEdge ** list,int left,int right,int (*comp)(GraphSeqEdge * ,GraphSeqEdge * ));
#define qsort_GraphSeq Wise2_qsort_GraphSeq
void Wise2_sort_GraphSeq(GraphSeq * obj,int (*comp)(GraphSeqEdge *, GraphSeqEdge *));
#define sort_GraphSeq Wise2_sort_GraphSeq
boolean Wise2_expand_GraphSeq(GraphSeq * obj,int len);
#define expand_GraphSeq Wise2_expand_GraphSeq

#ifdef _cplusplus
}
#endif

#endif
