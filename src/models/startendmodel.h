#ifndef DYNAMITEstartendmodelHEADERFILE
#define DYNAMITEstartendmodelHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"

struct Wise2_StartEndModelMode {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    char * style;    
    int i_ramp;  
    int j_ramp;  
    } ;  
/* StartEndModelMode defined */ 
#ifndef DYNAMITE_DEFINED_StartEndModelMode
typedef struct Wise2_StartEndModelMode Wise2_StartEndModelMode;
#define StartEndModelMode Wise2_StartEndModelMode
#define DYNAMITE_DEFINED_StartEndModelMode
#endif


struct Wise2_StartEndModel {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int i_len;   
    int j_len;   
    Probability * i_entry;   
    Probability * j_entry;   
    } ;  
/* StartEndModel defined */ 
#ifndef DYNAMITE_DEFINED_StartEndModel
typedef struct Wise2_StartEndModel Wise2_StartEndModel;
#define StartEndModel Wise2_StartEndModel
#define DYNAMITE_DEFINED_StartEndModel
#endif


struct Wise2_StartEndScore {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    int i_len;   
    int j_len;   
    Score * i_entry;     
    Score * j_entry;     
    } ;  
/* StartEndScore defined */ 
#ifndef DYNAMITE_DEFINED_StartEndScore
typedef struct Wise2_StartEndScore Wise2_StartEndScore;
#define StartEndScore Wise2_StartEndScore
#define DYNAMITE_DEFINED_StartEndScore
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_StartEndModelMode(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [StartEndModelMode *]
 *
 * Return [UNKN ]  Undocumented return value [StartEndModelMode *]
 *
 */
StartEndModelMode * Wise2_hard_link_StartEndModelMode(StartEndModelMode * obj);
#define hard_link_StartEndModelMode Wise2_hard_link_StartEndModelMode


/* Function:  StartEndModelMode_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StartEndModelMode *]
 *
 */
StartEndModelMode * Wise2_StartEndModelMode_alloc(void);
#define StartEndModelMode_alloc Wise2_StartEndModelMode_alloc


/* Function:  free_StartEndModelMode(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StartEndModelMode *]
 *
 * Return [UNKN ]  Undocumented return value [StartEndModelMode *]
 *
 */
StartEndModelMode * Wise2_free_StartEndModelMode(StartEndModelMode * obj);
#define free_StartEndModelMode Wise2_free_StartEndModelMode


/* Function:  hard_link_StartEndModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [StartEndModel *]
 *
 * Return [UNKN ]  Undocumented return value [StartEndModel *]
 *
 */
StartEndModel * Wise2_hard_link_StartEndModel(StartEndModel * obj);
#define hard_link_StartEndModel Wise2_hard_link_StartEndModel


/* Function:  StartEndModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StartEndModel *]
 *
 */
StartEndModel * Wise2_StartEndModel_alloc(void);
#define StartEndModel_alloc Wise2_StartEndModel_alloc


/* Function:  free_StartEndModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StartEndModel *]
 *
 * Return [UNKN ]  Undocumented return value [StartEndModel *]
 *
 */
StartEndModel * Wise2_free_StartEndModel(StartEndModel * obj);
#define free_StartEndModel Wise2_free_StartEndModel


/* Function:  hard_link_StartEndScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [StartEndScore *]
 *
 * Return [UNKN ]  Undocumented return value [StartEndScore *]
 *
 */
StartEndScore * Wise2_hard_link_StartEndScore(StartEndScore * obj);
#define hard_link_StartEndScore Wise2_hard_link_StartEndScore


/* Function:  StartEndScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StartEndScore *]
 *
 */
StartEndScore * Wise2_StartEndScore_alloc(void);
#define StartEndScore_alloc Wise2_StartEndScore_alloc


/* Function:  free_StartEndScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StartEndScore *]
 *
 * Return [UNKN ]  Undocumented return value [StartEndScore *]
 *
 */
StartEndScore * Wise2_free_StartEndScore(StartEndScore * obj);
#define free_StartEndScore Wise2_free_StartEndScore


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
void Wise2_free_Probability(Probability * p);
#define free_Probability Wise2_free_Probability
void Wise2_free_Score(Score * s);
#define free_Score Wise2_free_Score
StartEndModelMode * Wise2_StartEndModelMode_from_argc(int * argc,char ** argv);
#define StartEndModelMode_from_argc Wise2_StartEndModelMode_from_argc
StartEndScore * Wise2_StartEndScore_from_StartEndModel(StartEndModel * sem);
#define StartEndScore_from_StartEndModel Wise2_StartEndScore_from_StartEndModel


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
