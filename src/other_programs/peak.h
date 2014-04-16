#ifndef DYNAMITEpeakHEADERFILE
#define DYNAMITEpeakHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

#define PeakListLISTLENGTH 5000
#define MergedPeakLISTLENGTH 20
#define MergedPeakListLISTLENGTH 5000

#define ListofPeakListLISTLENGTH 20
struct Wise2_PeakPara {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * assay_name;   
    } ;  
/* PeakPara defined */ 
#ifndef DYNAMITE_DEFINED_PeakPara
typedef struct Wise2_PeakPara Wise2_PeakPara;
#define PeakPara Wise2_PeakPara
#define DYNAMITE_DEFINED_PeakPara
#endif


struct Wise2_Peak {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * chr;  
    int start;   
    int end;     
    double score;    
    char * id;   
    PeakPara * para;     
    } ;  
/* Peak defined */ 
#ifndef DYNAMITE_DEFINED_Peak
typedef struct Wise2_Peak Wise2_Peak;
#define Peak Wise2_Peak
#define DYNAMITE_DEFINED_Peak
#endif


struct Wise2_PeakList {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    PeakPara * para;     
    Peak ** peak;    
    int len;/* len for above peak  */ 
    int maxlen; /* maxlen for above peak */ 
    } ;  
/* PeakList defined */ 
#ifndef DYNAMITE_DEFINED_PeakList
typedef struct Wise2_PeakList Wise2_PeakList;
#define PeakList Wise2_PeakList
#define DYNAMITE_DEFINED_PeakList
#endif


struct Wise2_ListofPeakList {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    PeakList ** set;     
    int len;/* len for above set  */ 
    int maxlen; /* maxlen for above set */ 
    } ;  
/* ListofPeakList defined */ 
#ifndef DYNAMITE_DEFINED_ListofPeakList
typedef struct Wise2_ListofPeakList Wise2_ListofPeakList;
#define ListofPeakList Wise2_ListofPeakList
#define DYNAMITE_DEFINED_ListofPeakList
#endif


struct Wise2_MergedPeak {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * chr;  
    int start;   
    int end;     
    Peak ** comp;    
    int len;/* len for above comp  */ 
    int maxlen; /* maxlen for above comp */ 
    } ;  
/* MergedPeak defined */ 
#ifndef DYNAMITE_DEFINED_MergedPeak
typedef struct Wise2_MergedPeak Wise2_MergedPeak;
#define MergedPeak Wise2_MergedPeak
#define DYNAMITE_DEFINED_MergedPeak
#endif


struct Wise2_MergedPeakList {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    MergedPeak ** peak;  
    int len;/* len for above peak  */ 
    int maxlen; /* maxlen for above peak */ 
    } ;  
/* MergedPeakList defined */ 
#ifndef DYNAMITE_DEFINED_MergedPeakList
typedef struct Wise2_MergedPeakList Wise2_MergedPeakList;
#define MergedPeakList Wise2_MergedPeakList
#define DYNAMITE_DEFINED_MergedPeakList
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  write_tagged_MergedPeakList(mpl,ofp)
 *
 * Descrip:    writes tagged peak list
 *
 *
 * Arg:        mpl [UNKN ] Undocumented argument [MergedPeakList *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_write_tagged_MergedPeakList(MergedPeakList * mpl,FILE * ofp);
#define write_tagged_MergedPeakList Wise2_write_tagged_MergedPeakList


/* Function:  make_MergedPeakList(lop,stderr_logging_level)
 *
 * Descrip:    makes a merged peak list (single linkage clustering)
 *             from a list of peaklists
 *
 *
 * Arg:                         lop [UNKN ] Undocumented argument [ListofPeakList *]
 * Arg:        stderr_logging_level [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [MergedPeakList *]
 *
 */
MergedPeakList * Wise2_make_MergedPeakList(ListofPeakList * lop,int stderr_logging_level);
#define make_MergedPeakList Wise2_make_MergedPeakList


/* Function:  sort_PeakList_by_chr_start(pl)
 *
 * Descrip:    sorts a peak list
 *
 *
 * Arg:        pl [UNKN ] Undocumented argument [PeakList *]
 *
 */
void Wise2_sort_PeakList_by_chr_start(PeakList * pl);
#define sort_PeakList_by_chr_start Wise2_sort_PeakList_by_chr_start


/* Function:  comp_Peak(one,two)
 *
 * Descrip:    internal for sorting
 *             !internal
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [Peak *]
 * Arg:        two [UNKN ] Undocumented argument [Peak *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_comp_Peak(Peak * one,Peak * two ) ;
#define comp_Peak Wise2_comp_Peak


/* Function:  read_bed_PeakList(ifp,assay_name)
 *
 * Descrip:    creates a peak list from a 4 column bed file
 *
 *
 * Arg:               ifp [UNKN ] Undocumented argument [FILE *]
 * Arg:        assay_name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [PeakList *]
 *
 */
PeakList * Wise2_read_bed_PeakList(FILE * ifp,char * assay_name);
#define read_bed_PeakList Wise2_read_bed_PeakList


/* Function:  read_npf_PeakList(ifp,assay_name)
 *
 * Descrip:    creates a peak list from a 4 column npf file
 *
 *
 * Arg:               ifp [UNKN ] Undocumented argument [FILE *]
 * Arg:        assay_name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [PeakList *]
 *
 */
PeakList * Wise2_read_npf_PeakList(FILE * ifp,char * assay_name);
#define read_npf_PeakList Wise2_read_npf_PeakList


/* Function:  write_4col_bed_PeakList(pl,out)
 *
 * Descrip:    writes out a 4-column bed
 *
 *
 * Arg:         pl [UNKN ] Undocumented argument [PeakList *]
 * Arg:        out [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_write_4col_bed_PeakList(PeakList * pl,FILE * out);
#define write_4col_bed_PeakList Wise2_write_4col_bed_PeakList


/* Function:  hard_link_PeakPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PeakPara *]
 *
 * Return [UNKN ]  Undocumented return value [PeakPara *]
 *
 */
PeakPara * Wise2_hard_link_PeakPara(PeakPara * obj);
#define hard_link_PeakPara Wise2_hard_link_PeakPara


/* Function:  PeakPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PeakPara *]
 *
 */
PeakPara * Wise2_PeakPara_alloc(void);
#define PeakPara_alloc Wise2_PeakPara_alloc


/* Function:  free_PeakPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PeakPara *]
 *
 * Return [UNKN ]  Undocumented return value [PeakPara *]
 *
 */
PeakPara * Wise2_free_PeakPara(PeakPara * obj);
#define free_PeakPara Wise2_free_PeakPara


/* Function:  hard_link_Peak(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Peak *]
 *
 * Return [UNKN ]  Undocumented return value [Peak *]
 *
 */
Peak * Wise2_hard_link_Peak(Peak * obj);
#define hard_link_Peak Wise2_hard_link_Peak


/* Function:  Peak_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Peak *]
 *
 */
Peak * Wise2_Peak_alloc(void);
#define Peak_alloc Wise2_Peak_alloc


/* Function:  free_Peak(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Peak *]
 *
 * Return [UNKN ]  Undocumented return value [Peak *]
 *
 */
Peak * Wise2_free_Peak(Peak * obj);
#define free_Peak Wise2_free_Peak


/* Function:  add_PeakList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PeakList *]
 * Arg:        add [OWNER] Object to add to the list [Peak *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_PeakList(PeakList * obj,Peak * add);
#define add_PeakList Wise2_add_PeakList


/* Function:  flush_PeakList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [PeakList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_PeakList(PeakList * obj);
#define flush_PeakList Wise2_flush_PeakList


/* Function:  PeakList_alloc_std(void)
 *
 * Descrip:    Equivalent to PeakList_alloc_len(PeakListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PeakList *]
 *
 */
PeakList * Wise2_PeakList_alloc_std(void);
#define PeakList_alloc_std Wise2_PeakList_alloc_std


/* Function:  PeakList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [PeakList *]
 *
 */
PeakList * Wise2_PeakList_alloc_len(int len);
#define PeakList_alloc_len Wise2_PeakList_alloc_len


/* Function:  hard_link_PeakList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PeakList *]
 *
 * Return [UNKN ]  Undocumented return value [PeakList *]
 *
 */
PeakList * Wise2_hard_link_PeakList(PeakList * obj);
#define hard_link_PeakList Wise2_hard_link_PeakList


/* Function:  PeakList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PeakList *]
 *
 */
PeakList * Wise2_PeakList_alloc(void);
#define PeakList_alloc Wise2_PeakList_alloc


/* Function:  free_PeakList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PeakList *]
 *
 * Return [UNKN ]  Undocumented return value [PeakList *]
 *
 */
PeakList * Wise2_free_PeakList(PeakList * obj);
#define free_PeakList Wise2_free_PeakList


/* Function:  add_ListofPeakList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ListofPeakList *]
 * Arg:        add [OWNER] Object to add to the list [PeakList *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_ListofPeakList(ListofPeakList * obj,PeakList * add);
#define add_ListofPeakList Wise2_add_ListofPeakList


/* Function:  flush_ListofPeakList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ListofPeakList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_ListofPeakList(ListofPeakList * obj);
#define flush_ListofPeakList Wise2_flush_ListofPeakList


/* Function:  ListofPeakList_alloc_std(void)
 *
 * Descrip:    Equivalent to ListofPeakList_alloc_len(ListofPeakListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ListofPeakList *]
 *
 */
ListofPeakList * Wise2_ListofPeakList_alloc_std(void);
#define ListofPeakList_alloc_std Wise2_ListofPeakList_alloc_std


/* Function:  ListofPeakList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ListofPeakList *]
 *
 */
ListofPeakList * Wise2_ListofPeakList_alloc_len(int len);
#define ListofPeakList_alloc_len Wise2_ListofPeakList_alloc_len


/* Function:  hard_link_ListofPeakList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ListofPeakList *]
 *
 * Return [UNKN ]  Undocumented return value [ListofPeakList *]
 *
 */
ListofPeakList * Wise2_hard_link_ListofPeakList(ListofPeakList * obj);
#define hard_link_ListofPeakList Wise2_hard_link_ListofPeakList


/* Function:  ListofPeakList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ListofPeakList *]
 *
 */
ListofPeakList * Wise2_ListofPeakList_alloc(void);
#define ListofPeakList_alloc Wise2_ListofPeakList_alloc


/* Function:  free_ListofPeakList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ListofPeakList *]
 *
 * Return [UNKN ]  Undocumented return value [ListofPeakList *]
 *
 */
ListofPeakList * Wise2_free_ListofPeakList(ListofPeakList * obj);
#define free_ListofPeakList Wise2_free_ListofPeakList


/* Function:  add_MergedPeak(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MergedPeak *]
 * Arg:        add [OWNER] Object to add to the list [Peak *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_MergedPeak(MergedPeak * obj,Peak * add);
#define add_MergedPeak Wise2_add_MergedPeak


/* Function:  flush_MergedPeak(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MergedPeak *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_MergedPeak(MergedPeak * obj);
#define flush_MergedPeak Wise2_flush_MergedPeak


/* Function:  MergedPeak_alloc_std(void)
 *
 * Descrip:    Equivalent to MergedPeak_alloc_len(MergedPeakLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MergedPeak *]
 *
 */
MergedPeak * Wise2_MergedPeak_alloc_std(void);
#define MergedPeak_alloc_std Wise2_MergedPeak_alloc_std


/* Function:  MergedPeak_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [MergedPeak *]
 *
 */
MergedPeak * Wise2_MergedPeak_alloc_len(int len);
#define MergedPeak_alloc_len Wise2_MergedPeak_alloc_len


/* Function:  hard_link_MergedPeak(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MergedPeak *]
 *
 * Return [UNKN ]  Undocumented return value [MergedPeak *]
 *
 */
MergedPeak * Wise2_hard_link_MergedPeak(MergedPeak * obj);
#define hard_link_MergedPeak Wise2_hard_link_MergedPeak


/* Function:  MergedPeak_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MergedPeak *]
 *
 */
MergedPeak * Wise2_MergedPeak_alloc(void);
#define MergedPeak_alloc Wise2_MergedPeak_alloc


/* Function:  free_MergedPeak(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MergedPeak *]
 *
 * Return [UNKN ]  Undocumented return value [MergedPeak *]
 *
 */
MergedPeak * Wise2_free_MergedPeak(MergedPeak * obj);
#define free_MergedPeak Wise2_free_MergedPeak


/* Function:  add_MergedPeakList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MergedPeakList *]
 * Arg:        add [OWNER] Object to add to the list [MergedPeak *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_MergedPeakList(MergedPeakList * obj,MergedPeak * add);
#define add_MergedPeakList Wise2_add_MergedPeakList


/* Function:  flush_MergedPeakList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MergedPeakList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_MergedPeakList(MergedPeakList * obj);
#define flush_MergedPeakList Wise2_flush_MergedPeakList


/* Function:  MergedPeakList_alloc_std(void)
 *
 * Descrip:    Equivalent to MergedPeakList_alloc_len(MergedPeakListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MergedPeakList *]
 *
 */
MergedPeakList * Wise2_MergedPeakList_alloc_std(void);
#define MergedPeakList_alloc_std Wise2_MergedPeakList_alloc_std


/* Function:  MergedPeakList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [MergedPeakList *]
 *
 */
MergedPeakList * Wise2_MergedPeakList_alloc_len(int len);
#define MergedPeakList_alloc_len Wise2_MergedPeakList_alloc_len


/* Function:  hard_link_MergedPeakList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MergedPeakList *]
 *
 * Return [UNKN ]  Undocumented return value [MergedPeakList *]
 *
 */
MergedPeakList * Wise2_hard_link_MergedPeakList(MergedPeakList * obj);
#define hard_link_MergedPeakList Wise2_hard_link_MergedPeakList


/* Function:  MergedPeakList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MergedPeakList *]
 *
 */
MergedPeakList * Wise2_MergedPeakList_alloc(void);
#define MergedPeakList_alloc Wise2_MergedPeakList_alloc


/* Function:  free_MergedPeakList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MergedPeakList *]
 *
 * Return [UNKN ]  Undocumented return value [MergedPeakList *]
 *
 */
MergedPeakList * Wise2_free_MergedPeakList(MergedPeakList * obj);
#define free_MergedPeakList Wise2_free_MergedPeakList


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_PeakList(Peak ** list,int i,int j) ;
#define swap_PeakList Wise2_swap_PeakList
void Wise2_qsort_PeakList(Peak ** list,int left,int right,int (*comp)(Peak * ,Peak * ));
#define qsort_PeakList Wise2_qsort_PeakList
void Wise2_sort_PeakList(PeakList * obj,int (*comp)(Peak *, Peak *));
#define sort_PeakList Wise2_sort_PeakList
boolean Wise2_expand_PeakList(PeakList * obj,int len);
#define expand_PeakList Wise2_expand_PeakList
void Wise2_swap_ListofPeakList(PeakList ** list,int i,int j) ;
#define swap_ListofPeakList Wise2_swap_ListofPeakList
void Wise2_qsort_ListofPeakList(PeakList ** list,int left,int right,int (*comp)(PeakList * ,PeakList * ));
#define qsort_ListofPeakList Wise2_qsort_ListofPeakList
void Wise2_sort_ListofPeakList(ListofPeakList * obj,int (*comp)(PeakList *, PeakList *));
#define sort_ListofPeakList Wise2_sort_ListofPeakList
boolean Wise2_expand_ListofPeakList(ListofPeakList * obj,int len);
#define expand_ListofPeakList Wise2_expand_ListofPeakList
void Wise2_swap_MergedPeak(Peak ** list,int i,int j) ;
#define swap_MergedPeak Wise2_swap_MergedPeak
void Wise2_qsort_MergedPeak(Peak ** list,int left,int right,int (*comp)(Peak * ,Peak * ));
#define qsort_MergedPeak Wise2_qsort_MergedPeak
void Wise2_sort_MergedPeak(MergedPeak * obj,int (*comp)(Peak *, Peak *));
#define sort_MergedPeak Wise2_sort_MergedPeak
boolean Wise2_expand_MergedPeak(MergedPeak * obj,int len);
#define expand_MergedPeak Wise2_expand_MergedPeak
void Wise2_swap_MergedPeakList(MergedPeak ** list,int i,int j) ;
#define swap_MergedPeakList Wise2_swap_MergedPeakList
void Wise2_qsort_MergedPeakList(MergedPeak ** list,int left,int right,int (*comp)(MergedPeak * ,MergedPeak * ));
#define qsort_MergedPeakList Wise2_qsort_MergedPeakList
void Wise2_sort_MergedPeakList(MergedPeakList * obj,int (*comp)(MergedPeak *, MergedPeak *));
#define sort_MergedPeakList Wise2_sort_MergedPeakList
boolean Wise2_expand_MergedPeakList(MergedPeakList * obj,int len);
#define expand_MergedPeakList Wise2_expand_MergedPeakList

#ifdef _cplusplus
}
#endif

#endif
