#ifndef DYNAMITEtwo_trackHEADERFILE
#define DYNAMITEtwo_trackHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"


#define TwoTrackLISTLENGTH 80
#define TwoTrackScoreLISTLENGTH 80

#define TwoTrackSetLISTLENGTH 1024


typedef enum TwoTrackUnit_type {
  TwoTrackUnit_main = 67,
  TwoTrackUnit_side
} TwoTrackUnit_type;


typedef long int LongInt;

struct Wise2_TwoTrackUnit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    TwoTrackUnit_type type;  
    Probability emission[5];     
    } ;  
/* TwoTrackUnit defined */ 
#ifndef DYNAMITE_DEFINED_TwoTrackUnit
typedef struct Wise2_TwoTrackUnit Wise2_TwoTrackUnit;
#define TwoTrackUnit Wise2_TwoTrackUnit
#define DYNAMITE_DEFINED_TwoTrackUnit
#endif


struct Wise2_TwoTrack {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    Sequence * seq;  
    TwoTrackUnit ** unit;    
    int len;/* len for above unit  */ 
    int maxlen; /* maxlen for above unit */ 
    } ;  
/* TwoTrack defined */ 
#ifndef DYNAMITE_DEFINED_TwoTrack
typedef struct Wise2_TwoTrack Wise2_TwoTrack;
#define TwoTrack Wise2_TwoTrack
#define DYNAMITE_DEFINED_TwoTrack
#endif


struct Wise2_TwoTrackSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    TwoTrack ** read;    
    int len;/* len for above read  */ 
    int maxlen; /* maxlen for above read */ 
    } ;  
/* TwoTrackSet defined */ 
#ifndef DYNAMITE_DEFINED_TwoTrackSet
typedef struct Wise2_TwoTrackSet Wise2_TwoTrackSet;
#define TwoTrackSet Wise2_TwoTrackSet
#define DYNAMITE_DEFINED_TwoTrackSet
#endif


struct Wise2_TwoTrackScoreUnit {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    TwoTrackUnit_type type;  
    Score emission[5];   
    } ;  
/* TwoTrackScoreUnit defined */ 
#ifndef DYNAMITE_DEFINED_TwoTrackScoreUnit
typedef struct Wise2_TwoTrackScoreUnit Wise2_TwoTrackScoreUnit;
#define TwoTrackScoreUnit Wise2_TwoTrackScoreUnit
#define DYNAMITE_DEFINED_TwoTrackScoreUnit
#endif


struct Wise2_TwoTrackScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    TwoTrackScoreUnit ** unit;   
    int len;/* len for above unit  */ 
    int maxlen; /* maxlen for above unit */ 
    } ;  
/* TwoTrackScore defined */ 
#ifndef DYNAMITE_DEFINED_TwoTrackScore
typedef struct Wise2_TwoTrackScore Wise2_TwoTrackScore;
#define TwoTrackScore Wise2_TwoTrackScore
#define DYNAMITE_DEFINED_TwoTrackScore
#endif


struct Wise2_TwoTrackSetStats {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    LongInt reads;   
    LongInt main_bases;  
    LongInt side_bases;  
    double main_avg_likelihood;  
    double side_avg_likelihood;  
    LongInt called_main_bases[5];    
    LongInt most_likely_main_bases[5];   
    LongInt most_likely_side_bases[5];   
    } ;  
/* TwoTrackSetStats defined */ 
#ifndef DYNAMITE_DEFINED_TwoTrackSetStats
typedef struct Wise2_TwoTrackSetStats Wise2_TwoTrackSetStats;
#define TwoTrackSetStats Wise2_TwoTrackSetStats
#define DYNAMITE_DEFINED_TwoTrackSetStats
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  write_TwoTrackSetStats(st,ofp)
 *
 * Descrip:    writes out stats
 *
 *
 * Arg:         st [UNKN ] Undocumented argument [TwoTrackSetStats *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_write_TwoTrackSetStats(TwoTrackSetStats * st,FILE * ofp);
#define write_TwoTrackSetStats Wise2_write_TwoTrackSetStats


/* Function:  TwoTrackSetStats_from_TwoTrackSet(set)
 *
 * Descrip:    Makes stats from a set
 *
 *
 * Arg:        set [UNKN ] Undocumented argument [TwoTrackSet *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSetStats *]
 *
 */
TwoTrackSetStats * Wise2_TwoTrackSetStats_from_TwoTrackSet(TwoTrackSet * set);
#define TwoTrackSetStats_from_TwoTrackSet Wise2_TwoTrackSetStats_from_TwoTrackSet


/* Function:  TwoTrackScore_from_TwoTrack(t,emission_ratio)
 *
 * Descrip:    complete conversion
 *
 *
 * Arg:                     t [UNKN ] Undocumented argument [TwoTrack *]
 * Arg:        emission_ratio [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScore *]
 *
 */
TwoTrackScore * Wise2_TwoTrackScore_from_TwoTrack(TwoTrack * t,Probability emission_ratio);
#define TwoTrackScore_from_TwoTrack Wise2_TwoTrackScore_from_TwoTrack


/* Function:  TwoTrackScoreUnit_from_TwoTrackUnit(u,ratio)
 *
 * Descrip:    conversion
 *
 *
 * Arg:            u [UNKN ] Undocumented argument [TwoTrackUnit *]
 * Arg:        ratio [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScoreUnit *]
 *
 */
TwoTrackScoreUnit * Wise2_TwoTrackScoreUnit_from_TwoTrackUnit(TwoTrackUnit * u,Probability ratio);
#define TwoTrackScoreUnit_from_TwoTrackUnit Wise2_TwoTrackScoreUnit_from_TwoTrackUnit


/* Function:  write_plain_TwoTrackSet(set,ofp)
 *
 * Descrip:    plain output
 *
 *
 * Arg:        set [UNKN ] Undocumented argument [TwoTrackSet *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_write_plain_TwoTrackSet(TwoTrackSet * set,FILE * ofp);
#define write_plain_TwoTrackSet Wise2_write_plain_TwoTrackSet


/* Function:  TwoTrackUnit_type_to_string(type)
 *
 * Descrip:    type string
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [TwoTrackUnit_type]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_TwoTrackUnit_type_to_string(TwoTrackUnit_type type);
#define TwoTrackUnit_type_to_string Wise2_TwoTrackUnit_type_to_string


/* Function:  write_plain_TwoTrack(two,ofp)
 *
 * Descrip:    output in columns - base, and 4 likes
 *
 *
 *
 * Arg:        two [UNKN ] Undocumented argument [TwoTrack *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_write_plain_TwoTrack(TwoTrack * two,FILE * ofp);
#define write_plain_TwoTrack Wise2_write_plain_TwoTrack


/* Function:  read_two_track_Tim_style(seq_file,like_file)
 *
 * Descrip:    Reads a two track file with a default reporting of every
 *             1000 sequences read and the entire file
 *
 *
 * Arg:         seq_file [UNKN ] Undocumented argument [char *]
 * Arg:        like_file [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSet *]
 *
 */
TwoTrackSet * Wise2_read_two_track_Tim_style(char * seq_file,char * like_file);
#define read_two_track_Tim_style Wise2_read_two_track_Tim_style


/* Function:  read_two_track_Tim_style_report(seq_file,like_file,report_freq,truncate_after)
 *
 * Descrip:    Reads a two track file from Solexa style information
 *             from Tim Massinghams lik files
 *
 *             Assummes all positions are main emissions
 *
 *             report freq means the function will report every xx lines
 *
 *
 * Arg:              seq_file [UNKN ] Undocumented argument [char *]
 * Arg:             like_file [UNKN ] Undocumented argument [char *]
 * Arg:           report_freq [UNKN ] Undocumented argument [int]
 * Arg:        truncate_after [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSet *]
 *
 */
TwoTrackSet * Wise2_read_two_track_Tim_style_report(char * seq_file,char * like_file,int report_freq, int truncate_after);
#define read_two_track_Tim_style_report Wise2_read_two_track_Tim_style_report


/* Function:  hard_link_TwoTrackUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TwoTrackUnit *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackUnit *]
 *
 */
TwoTrackUnit * Wise2_hard_link_TwoTrackUnit(TwoTrackUnit * obj);
#define hard_link_TwoTrackUnit Wise2_hard_link_TwoTrackUnit


/* Function:  TwoTrackUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackUnit *]
 *
 */
TwoTrackUnit * Wise2_TwoTrackUnit_alloc(void);
#define TwoTrackUnit_alloc Wise2_TwoTrackUnit_alloc


/* Function:  free_TwoTrackUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TwoTrackUnit *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackUnit *]
 *
 */
TwoTrackUnit * Wise2_free_TwoTrackUnit(TwoTrackUnit * obj);
#define free_TwoTrackUnit Wise2_free_TwoTrackUnit


/* Function:  add_TwoTrack(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TwoTrack *]
 * Arg:        add [OWNER] Object to add to the list [TwoTrackUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_TwoTrack(TwoTrack * obj,TwoTrackUnit * add);
#define add_TwoTrack Wise2_add_TwoTrack


/* Function:  flush_TwoTrack(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TwoTrack *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_TwoTrack(TwoTrack * obj);
#define flush_TwoTrack Wise2_flush_TwoTrack


/* Function:  TwoTrack_alloc_std(void)
 *
 * Descrip:    Equivalent to TwoTrack_alloc_len(TwoTrackLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrack *]
 *
 */
TwoTrack * Wise2_TwoTrack_alloc_std(void);
#define TwoTrack_alloc_std Wise2_TwoTrack_alloc_std


/* Function:  TwoTrack_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrack *]
 *
 */
TwoTrack * Wise2_TwoTrack_alloc_len(int len);
#define TwoTrack_alloc_len Wise2_TwoTrack_alloc_len


/* Function:  hard_link_TwoTrack(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TwoTrack *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrack *]
 *
 */
TwoTrack * Wise2_hard_link_TwoTrack(TwoTrack * obj);
#define hard_link_TwoTrack Wise2_hard_link_TwoTrack


/* Function:  TwoTrack_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrack *]
 *
 */
TwoTrack * Wise2_TwoTrack_alloc(void);
#define TwoTrack_alloc Wise2_TwoTrack_alloc


/* Function:  free_TwoTrack(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TwoTrack *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrack *]
 *
 */
TwoTrack * Wise2_free_TwoTrack(TwoTrack * obj);
#define free_TwoTrack Wise2_free_TwoTrack


/* Function:  add_TwoTrackSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TwoTrackSet *]
 * Arg:        add [OWNER] Object to add to the list [TwoTrack *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_TwoTrackSet(TwoTrackSet * obj,TwoTrack * add);
#define add_TwoTrackSet Wise2_add_TwoTrackSet


/* Function:  flush_TwoTrackSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TwoTrackSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_TwoTrackSet(TwoTrackSet * obj);
#define flush_TwoTrackSet Wise2_flush_TwoTrackSet


/* Function:  TwoTrackSet_alloc_std(void)
 *
 * Descrip:    Equivalent to TwoTrackSet_alloc_len(TwoTrackSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSet *]
 *
 */
TwoTrackSet * Wise2_TwoTrackSet_alloc_std(void);
#define TwoTrackSet_alloc_std Wise2_TwoTrackSet_alloc_std


/* Function:  TwoTrackSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSet *]
 *
 */
TwoTrackSet * Wise2_TwoTrackSet_alloc_len(int len);
#define TwoTrackSet_alloc_len Wise2_TwoTrackSet_alloc_len


/* Function:  hard_link_TwoTrackSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TwoTrackSet *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSet *]
 *
 */
TwoTrackSet * Wise2_hard_link_TwoTrackSet(TwoTrackSet * obj);
#define hard_link_TwoTrackSet Wise2_hard_link_TwoTrackSet


/* Function:  TwoTrackSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSet *]
 *
 */
TwoTrackSet * Wise2_TwoTrackSet_alloc(void);
#define TwoTrackSet_alloc Wise2_TwoTrackSet_alloc


/* Function:  free_TwoTrackSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TwoTrackSet *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSet *]
 *
 */
TwoTrackSet * Wise2_free_TwoTrackSet(TwoTrackSet * obj);
#define free_TwoTrackSet Wise2_free_TwoTrackSet


/* Function:  hard_link_TwoTrackScoreUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TwoTrackScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScoreUnit *]
 *
 */
TwoTrackScoreUnit * Wise2_hard_link_TwoTrackScoreUnit(TwoTrackScoreUnit * obj);
#define hard_link_TwoTrackScoreUnit Wise2_hard_link_TwoTrackScoreUnit


/* Function:  TwoTrackScoreUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScoreUnit *]
 *
 */
TwoTrackScoreUnit * Wise2_TwoTrackScoreUnit_alloc(void);
#define TwoTrackScoreUnit_alloc Wise2_TwoTrackScoreUnit_alloc


/* Function:  free_TwoTrackScoreUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TwoTrackScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScoreUnit *]
 *
 */
TwoTrackScoreUnit * Wise2_free_TwoTrackScoreUnit(TwoTrackScoreUnit * obj);
#define free_TwoTrackScoreUnit Wise2_free_TwoTrackScoreUnit


/* Function:  add_TwoTrackScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TwoTrackScore *]
 * Arg:        add [OWNER] Object to add to the list [TwoTrackScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_TwoTrackScore(TwoTrackScore * obj,TwoTrackScoreUnit * add);
#define add_TwoTrackScore Wise2_add_TwoTrackScore


/* Function:  flush_TwoTrackScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TwoTrackScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_TwoTrackScore(TwoTrackScore * obj);
#define flush_TwoTrackScore Wise2_flush_TwoTrackScore


/* Function:  TwoTrackScore_alloc_std(void)
 *
 * Descrip:    Equivalent to TwoTrackScore_alloc_len(TwoTrackScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScore *]
 *
 */
TwoTrackScore * Wise2_TwoTrackScore_alloc_std(void);
#define TwoTrackScore_alloc_std Wise2_TwoTrackScore_alloc_std


/* Function:  TwoTrackScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScore *]
 *
 */
TwoTrackScore * Wise2_TwoTrackScore_alloc_len(int len);
#define TwoTrackScore_alloc_len Wise2_TwoTrackScore_alloc_len


/* Function:  hard_link_TwoTrackScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TwoTrackScore *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScore *]
 *
 */
TwoTrackScore * Wise2_hard_link_TwoTrackScore(TwoTrackScore * obj);
#define hard_link_TwoTrackScore Wise2_hard_link_TwoTrackScore


/* Function:  TwoTrackScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScore *]
 *
 */
TwoTrackScore * Wise2_TwoTrackScore_alloc(void);
#define TwoTrackScore_alloc Wise2_TwoTrackScore_alloc


/* Function:  free_TwoTrackScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TwoTrackScore *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScore *]
 *
 */
TwoTrackScore * Wise2_free_TwoTrackScore(TwoTrackScore * obj);
#define free_TwoTrackScore Wise2_free_TwoTrackScore


/* Function:  hard_link_TwoTrackSetStats(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TwoTrackSetStats *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSetStats *]
 *
 */
TwoTrackSetStats * Wise2_hard_link_TwoTrackSetStats(TwoTrackSetStats * obj);
#define hard_link_TwoTrackSetStats Wise2_hard_link_TwoTrackSetStats


/* Function:  TwoTrackSetStats_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSetStats *]
 *
 */
TwoTrackSetStats * Wise2_TwoTrackSetStats_alloc(void);
#define TwoTrackSetStats_alloc Wise2_TwoTrackSetStats_alloc


/* Function:  free_TwoTrackSetStats(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TwoTrackSetStats *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSetStats *]
 *
 */
TwoTrackSetStats * Wise2_free_TwoTrackSetStats(TwoTrackSetStats * obj);
#define free_TwoTrackSetStats Wise2_free_TwoTrackSetStats


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
base Wise2_most_likely_base_TwoTrackUnit(TwoTrackUnit * u);
#define most_likely_base_TwoTrackUnit Wise2_most_likely_base_TwoTrackUnit


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_TwoTrack(TwoTrackUnit ** list,int i,int j) ;
#define swap_TwoTrack Wise2_swap_TwoTrack
void Wise2_qsort_TwoTrack(TwoTrackUnit ** list,int left,int right,int (*comp)(TwoTrackUnit * ,TwoTrackUnit * ));
#define qsort_TwoTrack Wise2_qsort_TwoTrack
void Wise2_sort_TwoTrack(TwoTrack * obj,int (*comp)(TwoTrackUnit *, TwoTrackUnit *));
#define sort_TwoTrack Wise2_sort_TwoTrack
boolean Wise2_expand_TwoTrack(TwoTrack * obj,int len);
#define expand_TwoTrack Wise2_expand_TwoTrack
void Wise2_swap_TwoTrackSet(TwoTrack ** list,int i,int j) ;
#define swap_TwoTrackSet Wise2_swap_TwoTrackSet
void Wise2_qsort_TwoTrackSet(TwoTrack ** list,int left,int right,int (*comp)(TwoTrack * ,TwoTrack * ));
#define qsort_TwoTrackSet Wise2_qsort_TwoTrackSet
void Wise2_sort_TwoTrackSet(TwoTrackSet * obj,int (*comp)(TwoTrack *, TwoTrack *));
#define sort_TwoTrackSet Wise2_sort_TwoTrackSet
boolean Wise2_expand_TwoTrackSet(TwoTrackSet * obj,int len);
#define expand_TwoTrackSet Wise2_expand_TwoTrackSet
void Wise2_swap_TwoTrackScore(TwoTrackScoreUnit ** list,int i,int j) ;
#define swap_TwoTrackScore Wise2_swap_TwoTrackScore
void Wise2_qsort_TwoTrackScore(TwoTrackScoreUnit ** list,int left,int right,int (*comp)(TwoTrackScoreUnit * ,TwoTrackScoreUnit * ));
#define qsort_TwoTrackScore Wise2_qsort_TwoTrackScore
void Wise2_sort_TwoTrackScore(TwoTrackScore * obj,int (*comp)(TwoTrackScoreUnit *, TwoTrackScoreUnit *));
#define sort_TwoTrackScore Wise2_sort_TwoTrackScore
boolean Wise2_expand_TwoTrackScore(TwoTrackScore * obj,int len);
#define expand_TwoTrackScore Wise2_expand_TwoTrackScore

#ifdef _cplusplus
}
#endif

#endif
