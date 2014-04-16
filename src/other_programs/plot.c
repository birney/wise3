#ifdef _cplusplus
extern "C" {
#endif
#include "plot.h"


# line 44 "plot.dy"
Colour * std_Colour(char * colour_string)
{
  Colour * out;

  out = Colour_alloc();
  
  if( strcmp(colour_string,"red") == 0 ) {
    out->red = 1.0;
    out->green = 0.0;
    out->blue = 0.0;
  } else if ( strcmp(colour_string,"green") == 0 ) {
    out->red = 0.0;
    out->green = 1.0;
    out->blue = 0.0;
  } else if ( strcmp(colour_string,"blue") == 0 ) {
    out->red = 0.0;
    out->green = 0.0;
    out->blue = 1.0;
  } else {
    out->red = 1.0;
    out->green = 1.0;
    out->blue = 1.0;
  }


  return(out);
}

# line 72 "plot.dy"
BoxGlyph * filled_rect_BoxGlyph(double min_x,double min_y,double max_x,double max_y,Colour * fill,Colour * stroke)
{
  BoxGlyph * out;

  out = BoxGlyph_alloc();

  out->fill_colour = fill;
  out->stroke_colour = stroke;

  out->box  = new_BoundingBox(min_x,min_y,max_x,max_y);
  out->type = BoxGlyphType_filled_rect;
  return(out);
}

# line 86 "plot.dy"
BoxGlyph * text_BoxGlyph(double x,double y,char * string,double point_size)
{
  BoxGlyph * out;

  out = BoxGlyph_alloc();

  out->string = stringalloc(string);
  out->point_size = point_size;
  out->box = new_BoundingBox(x,y,x,y);
  out->type = BoxGlyphType_text;
  return(out);
}



# line 101 "plot.dy"
BoundingBox * new_BoundingBox(double min_x,double min_y,double max_x,double max_y)
{
  BoundingBox * out;

  out = BoundingBox_alloc();
  out->min_x = min_x;
  out->min_y = min_y;
  out->max_x = max_x;
  out->max_y = max_y;

  return(out);
}

# line 114 "plot.dy"
void postscript_BoxGlyph(BoxGlyph * b,FILE * ofp)
{
  assert(b != NULL);
  assert(ofp != NULL);

  switch(b->type) {
  case BoxGlyphType_filled_rect :
    postscript_BoxGlyph_filled_rectangle(b,ofp);
    return;
  case BoxGlyphType_text :
    postscript_BoxGlyph_text(b,ofp);
    return;
  default :
    assert("No Box glyph of this type!" == 0);
    
  }

  assert("Dropped through switch" == 0);
  return;
}

# line 135 "plot.dy"
void postscript_BoxGlyph_filled_rectangle(BoxGlyph * b,FILE * ofp)
{
  assert(b != NULL);
  assert(ofp != NULL);

  fprintf(ofp,"gsave\n");
  fprintf(ofp,"%f %f %f setrgbcolor\n",b->fill_colour->red,
	  b->fill_colour->green,
	  b->fill_colour->blue);
  fprintf(ofp,"%f %f moveto\n",b->box->min_x,b->box->min_y);
  fprintf(ofp,"%f %f lineto\n",b->box->min_x,b->box->max_y);
  fprintf(ofp,"%f %f lineto\n",b->box->max_x,b->box->max_y);
  fprintf(ofp,"%f %f lineto\n",b->box->max_x,b->box->min_y);
  fprintf(ofp,"closepath fill\n");
  
  fprintf(ofp,"grestore\n");

}

# line 154 "plot.dy"
void postscript_BoxGlyph_text(BoxGlyph * b,FILE * ofp)
{

  assert(b != NULL);
  assert(ofp != NULL);


  fprintf(ofp,"/Helvetica findfont %f scalefont setfont\n",b->point_size);

  fprintf(ofp,"%f %f moveto (%s) show\n",b->box->min_x,b->box->min_y,b->string);


}


# line 169 "plot.dy"
void flat_no_frame_postscript_FrameSet(FrameSet * fs,FILE * ofp)
{
  int i;

  assert(fs != NULL);
  assert(ofp != NULL);

  fprintf(ofp,"%%!PS-Adobe-2.0\n");
  fprintf(ofp,"% Created from flat plot\n");
  

  for(i=0;i<fs->len;i++) {
    postscript_BoxGlyph(fs->boxglyph[i],ofp);
  }
   
}



# line 158 "plot.c"
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
Colour * hard_link_Colour(Colour * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Colour object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Colour_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Colour *]
 *
 */
Colour * Colour_alloc(void) 
{
    Colour * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Colour *) ckalloc (sizeof(Colour))) == NULL)    {  
      warn("Colour_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->red = 0;    
    out->green = 0;  
    out->blue = 0;   


    return out;  
}    


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
Colour * free_Colour(Colour * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Colour obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   


    ckfree(obj); 
    return NULL; 
}    


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
BoundingBox * hard_link_BoundingBox(BoundingBox * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a BoundingBox object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  BoundingBox_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BoundingBox *]
 *
 */
BoundingBox * BoundingBox_alloc(void) 
{
    BoundingBox * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(BoundingBox *) ckalloc (sizeof(BoundingBox))) == NULL)  {  
      warn("BoundingBox_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->min_x = 0;  
    out->min_y = 0;  
    out->max_x = 0;  
    out->max_y = 0;  


    return out;  
}    


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
BoundingBox * free_BoundingBox(BoundingBox * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a BoundingBox obj. Should be trappable");   
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   


    ckfree(obj); 
    return NULL; 
}    


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
BoxGlyph * hard_link_BoxGlyph(BoxGlyph * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a BoxGlyph object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  BoxGlyph_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BoxGlyph *]
 *
 */
BoxGlyph * BoxGlyph_alloc(void) 
{
    BoxGlyph * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(BoxGlyph *) ckalloc (sizeof(BoxGlyph))) == NULL)    {  
      warn("BoxGlyph_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = 0;   
    out->box = NULL; 
    out->string = NULL;  
    out->fill_colour = NULL; 
    out->stroke_colour = NULL;   
    out->point_size = 0; 


    return out;  
}    


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
BoxGlyph * free_BoxGlyph(BoxGlyph * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a BoxGlyph obj. Should be trappable");  
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->box != NULL)    
      free_BoundingBox(obj->box);    
    if( obj->string != NULL) 
      ckfree(obj->string);   
    if( obj->fill_colour != NULL)    
      free_Colour(obj->fill_colour);     
    if( obj->stroke_colour != NULL)  
      free_Colour(obj->stroke_colour);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_FrameSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_FrameSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [BoxGlyph **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_FrameSet(BoxGlyph ** list,int i,int j)  
{
    BoxGlyph * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_FrameSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_FrameSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [BoxGlyph **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_FrameSet(BoxGlyph ** list,int left,int right,int (*comp)(BoxGlyph * ,BoxGlyph * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_FrameSet(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_FrameSet (list,++last,i);   
      }  
    swap_FrameSet (list,left,last);  
    qsort_FrameSet(list,left,last-1,comp);   
    qsort_FrameSet(list,last+1,right,comp);  
}    


/* Function:  sort_FrameSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_FrameSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [FrameSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_FrameSet(FrameSet * obj,int (*comp)(BoxGlyph *, BoxGlyph *)) 
{
    qsort_FrameSet(obj->boxglyph,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_FrameSet(obj,len)
 *
 * Descrip:    Really an internal function for add_FrameSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [FrameSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_FrameSet(FrameSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_FrameSet called with no need");   
      return TRUE;   
      }  


    if( (obj->boxglyph = (BoxGlyph ** ) ckrealloc (obj->boxglyph,sizeof(BoxGlyph *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_FrameSet, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


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
/* will expand function if necessary */ 
boolean add_FrameSet(FrameSet * obj,BoxGlyph * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_FrameSet(obj,obj->len + FrameSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->boxglyph[obj->len++]=add;   
    return TRUE; 
}    


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
int flush_FrameSet(FrameSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->boxglyph[i] != NULL)  {  
        free_BoxGlyph(obj->boxglyph[i]); 
        obj->boxglyph[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  FrameSet_alloc_std(void)
 *
 * Descrip:    Equivalent to FrameSet_alloc_len(FrameSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FrameSet *]
 *
 */
FrameSet * FrameSet_alloc_std(void) 
{
    return FrameSet_alloc_len(FrameSetLISTLENGTH);   
}    


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
FrameSet * FrameSet_alloc_len(int len) 
{
    FrameSet * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = FrameSet_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->boxglyph = (BoxGlyph ** ) ckcalloc (len,sizeof(BoxGlyph *))) == NULL)   {  
      warn("Warning, ckcalloc failed in FrameSet_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


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
FrameSet * hard_link_FrameSet(FrameSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a FrameSet object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  FrameSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FrameSet *]
 *
 */
FrameSet * FrameSet_alloc(void) 
{
    FrameSet * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FrameSet *) ckalloc (sizeof(FrameSet))) == NULL)    {  
      warn("FrameSet_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->frame = NULL;   
    out->boxglyph = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


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
FrameSet * free_FrameSet(FrameSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a FrameSet obj. Should be trappable");  
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->frame != NULL)  
      free_BoundingBox(obj->frame);  
    if( obj->boxglyph != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->boxglyph[i] != NULL)    
          free_BoxGlyph(obj->boxglyph[i]);   
        }  
      ckfree(obj->boxglyph); 
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
