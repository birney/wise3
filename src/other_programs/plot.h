#ifndef DYNAMITEplotHEADERFILE
#define DYNAMITEplotHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "dyna.h"


typedef enum BoxGlyphType_enum {
	BoxGlyphType_filled_rect = 56,
	BoxGlyphType_circle,
	BoxGlyphType_text
	} BoxGlyphType;

#define FrameSetLISTLENGTH 128
struct Wise2_Colour {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    double red;  
    double green;    
    double blue;     
    } ;  
/* Colour defined */ 
#ifndef DYNAMITE_DEFINED_Colour
typedef struct Wise2_Colour Wise2_Colour;
#define Colour Wise2_Colour
#define DYNAMITE_DEFINED_Colour
#endif


struct Wise2_BoundingBox {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    double min_x;    
    double min_y;    
    double max_x;    
    double max_y;    
    } ;  
/* BoundingBox defined */ 
#ifndef DYNAMITE_DEFINED_BoundingBox
typedef struct Wise2_BoundingBox Wise2_BoundingBox;
#define BoundingBox Wise2_BoundingBox
#define DYNAMITE_DEFINED_BoundingBox
#endif


struct Wise2_BoxGlyph {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BoxGlyphType type;   
    BoundingBox * box;   
    char * string;   
    Colour * fill_colour;    
    Colour * stroke_colour;  
    double point_size;   
    } ;  
/* BoxGlyph defined */ 
#ifndef DYNAMITE_DEFINED_BoxGlyph
typedef struct Wise2_BoxGlyph Wise2_BoxGlyph;
#define BoxGlyph Wise2_BoxGlyph
#define DYNAMITE_DEFINED_BoxGlyph
#endif


struct Wise2_FrameSet {  
    int dynamite_hard_link;  
#ifdef PTHREAD   
    pthread_mutex_t dynamite_mutex;  
#endif   
    BoundingBox * frame;     
    BoxGlyph ** boxglyph;    
    int len;/* len for above boxglyph  */ 
    int maxlen; /* maxlen for above boxglyph */ 
    } ;  
/* FrameSet defined */ 
#ifndef DYNAMITE_DEFINED_FrameSet
typedef struct Wise2_FrameSet Wise2_FrameSet;
#define FrameSet Wise2_FrameSet
#define DYNAMITE_DEFINED_FrameSet
#endif




    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  hard_link_Colour(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Colour *]
 *
 * Return [UNKN ]  Undocumented return value [Colour *]
 *
 */
Colour * Wise2_hard_link_Colour(Colour * obj);
#define hard_link_Colour Wise2_hard_link_Colour


/* Function:  Colour_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Colour *]
 *
 */
Colour * Wise2_Colour_alloc(void);
#define Colour_alloc Wise2_Colour_alloc


/* Function:  free_Colour(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Colour *]
 *
 * Return [UNKN ]  Undocumented return value [Colour *]
 *
 */
Colour * Wise2_free_Colour(Colour * obj);
#define free_Colour Wise2_free_Colour


/* Function:  hard_link_BoundingBox(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [BoundingBox *]
 *
 * Return [UNKN ]  Undocumented return value [BoundingBox *]
 *
 */
BoundingBox * Wise2_hard_link_BoundingBox(BoundingBox * obj);
#define hard_link_BoundingBox Wise2_hard_link_BoundingBox


/* Function:  BoundingBox_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BoundingBox *]
 *
 */
BoundingBox * Wise2_BoundingBox_alloc(void);
#define BoundingBox_alloc Wise2_BoundingBox_alloc


/* Function:  free_BoundingBox(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [BoundingBox *]
 *
 * Return [UNKN ]  Undocumented return value [BoundingBox *]
 *
 */
BoundingBox * Wise2_free_BoundingBox(BoundingBox * obj);
#define free_BoundingBox Wise2_free_BoundingBox


/* Function:  hard_link_BoxGlyph(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [BoxGlyph *]
 *
 * Return [UNKN ]  Undocumented return value [BoxGlyph *]
 *
 */
BoxGlyph * Wise2_hard_link_BoxGlyph(BoxGlyph * obj);
#define hard_link_BoxGlyph Wise2_hard_link_BoxGlyph


/* Function:  BoxGlyph_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BoxGlyph *]
 *
 */
BoxGlyph * Wise2_BoxGlyph_alloc(void);
#define BoxGlyph_alloc Wise2_BoxGlyph_alloc


/* Function:  free_BoxGlyph(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [BoxGlyph *]
 *
 * Return [UNKN ]  Undocumented return value [BoxGlyph *]
 *
 */
BoxGlyph * Wise2_free_BoxGlyph(BoxGlyph * obj);
#define free_BoxGlyph Wise2_free_BoxGlyph


/* Function:  add_FrameSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FrameSet *]
 * Arg:        add [OWNER] Object to add to the list [BoxGlyph *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_add_FrameSet(FrameSet * obj,BoxGlyph * add);
#define add_FrameSet Wise2_add_FrameSet


/* Function:  flush_FrameSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [FrameSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_flush_FrameSet(FrameSet * obj);
#define flush_FrameSet Wise2_flush_FrameSet


/* Function:  FrameSet_alloc_std(void)
 *
 * Descrip:    Equivalent to FrameSet_alloc_len(FrameSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FrameSet *]
 *
 */
FrameSet * Wise2_FrameSet_alloc_std(void);
#define FrameSet_alloc_std Wise2_FrameSet_alloc_std


/* Function:  FrameSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [FrameSet *]
 *
 */
FrameSet * Wise2_FrameSet_alloc_len(int len);
#define FrameSet_alloc_len Wise2_FrameSet_alloc_len


/* Function:  hard_link_FrameSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [FrameSet *]
 *
 * Return [UNKN ]  Undocumented return value [FrameSet *]
 *
 */
FrameSet * Wise2_hard_link_FrameSet(FrameSet * obj);
#define hard_link_FrameSet Wise2_hard_link_FrameSet


/* Function:  FrameSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FrameSet *]
 *
 */
FrameSet * Wise2_FrameSet_alloc(void);
#define FrameSet_alloc Wise2_FrameSet_alloc


/* Function:  free_FrameSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FrameSet *]
 *
 * Return [UNKN ]  Undocumented return value [FrameSet *]
 *
 */
FrameSet * Wise2_free_FrameSet(FrameSet * obj);
#define free_FrameSet Wise2_free_FrameSet


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */
Colour * Wise2_std_Colour(char * colour_string);
#define std_Colour Wise2_std_Colour
BoxGlyph * Wise2_filled_rect_BoxGlyph(double min_x,double min_y,double max_x,double max_y,Colour * fill,Colour * stroke);
#define filled_rect_BoxGlyph Wise2_filled_rect_BoxGlyph
BoxGlyph * Wise2_text_BoxGlyph(double x,double y,char * string,double point_size);
#define text_BoxGlyph Wise2_text_BoxGlyph
BoundingBox * Wise2_new_BoundingBox(double min_x,double min_y,double max_x,double max_y);
#define new_BoundingBox Wise2_new_BoundingBox
void Wise2_postscript_BoxGlyph(BoxGlyph * b,FILE * ofp);
#define postscript_BoxGlyph Wise2_postscript_BoxGlyph
void Wise2_postscript_BoxGlyph_filled_rectangle(BoxGlyph * b,FILE * ofp);
#define postscript_BoxGlyph_filled_rectangle Wise2_postscript_BoxGlyph_filled_rectangle
void Wise2_postscript_BoxGlyph_text(BoxGlyph * b,FILE * ofp);
#define postscript_BoxGlyph_text Wise2_postscript_BoxGlyph_text
void Wise2_flat_no_frame_postscript_FrameSet(FrameSet * fs,FILE * ofp);
#define flat_no_frame_postscript_FrameSet Wise2_flat_no_frame_postscript_FrameSet


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
void Wise2_swap_FrameSet(BoxGlyph ** list,int i,int j) ;
#define swap_FrameSet Wise2_swap_FrameSet
void Wise2_qsort_FrameSet(BoxGlyph ** list,int left,int right,int (*comp)(BoxGlyph * ,BoxGlyph * ));
#define qsort_FrameSet Wise2_qsort_FrameSet
void Wise2_sort_FrameSet(FrameSet * obj,int (*comp)(BoxGlyph *, BoxGlyph *));
#define sort_FrameSet Wise2_sort_FrameSet
boolean Wise2_expand_FrameSet(FrameSet * obj,int len);
#define expand_FrameSet Wise2_expand_FrameSet

#ifdef _cplusplus
}
#endif

#endif
