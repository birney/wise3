#ifdef _cplusplus
extern "C" {
#endif
#include "dnamapping.h"


/* Function:  map_to_basepair_numbers(seq,len)
 *
 * Descrip:    Maps sequence 
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [char *]
 * Arg:        len [UNKN ] Undocumented argument [long int]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 24 "dnamapping.dy"
char * map_to_basepair_numbers(char * seq,long int len)
{
  long int i;
  char * out;

  out = calloc(len,sizeof(char));

  for(i=0;i<len;i++) {
    out[i] = base_from_char(toupper(seq[i]));
  }

  return out;
}

/* Function:  reverse_complement_dna_number(number,kmer_size)
 *
 * Descrip:    makes the reverse complement of a kmer number
 *
 *
 * Arg:           number [UNKN ] Undocumented argument [kmer_t]
 * Arg:        kmer_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [kmer_t]
 *
 */
# line 41 "dnamapping.dy"
kmer_t  reverse_complement_dna_number(kmer_t number,int kmer_size)
{
  kmer_t base = 1;
  int base_number[128];
  int i;
  kmer_t total =0;

  for(i=0;i<kmer_size-1;i++) {
    base = base* 4;
  }

  for(i=kmer_size-1;i >= 0;i--) {
    base_number[i] = number / base;
    number -= base*base_number[i];
    base = base / 4;
  }

  base = 1;
  total = 0;
  for(i=0;i<kmer_size;i++) {
    total += complement_base(base_number[kmer_size-i-1]) * base;
    base = base * 4;
  }

  return total;

}

/* Function:  reverse_map_dna_number(number,kmer_size,buffer)
 *
 * Descrip:    puts a sequence into the buffer of the number
 *
 *
 * Arg:           number [UNKN ] Undocumented argument [kmer_t]
 * Arg:        kmer_size [UNKN ] Undocumented argument [int]
 * Arg:           buffer [UNKN ] Undocumented argument [char *]
 *
 */
# line 72 "dnamapping.dy"
void reverse_map_dna_number(kmer_t number,int kmer_size,char * buffer)
{
  int i;
  kmer_t base = 1;
  kmer_t no;
  
  for(i=0;i<kmer_size-1;i++) {
    base = base* 4;
  }

  for(i=kmer_size-1;i >= 0;i--) {
    no = number / base;
    buffer[i] = char_from_base(no);
    number -= base*no;
    base = base / 4;
  }


  return;
}

/* Function:  lexical_last_base_from_kmer(kmer,kmer_size,lexical_for)
 *
 * Descrip:    Find the lexically last base from kmer
 *
 *
 * Arg:               kmer [UNKN ] Undocumented argument [kmer_t]
 * Arg:          kmer_size [UNKN ] Undocumented argument [int]
 * Arg:        lexical_for [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
# line 96 "dnamapping.dy"
char lexical_last_base_from_kmer(kmer_t kmer,int kmer_size,int lexical_for)
{
  int i;
  kmer_t base = 1;
  int last_base;

  if( lexical_for == 1) {
    for(i=0;i<kmer_size-1;i++) {
      base = base* 4;
    }
    last_base = (int) (kmer / base);
    return char_from_base(last_base);
  } else {
    last_base = kmer % 4;
    return char_from_base(complement_base(last_base));
  }


}

/* Function:  lexical_dna_number_from_string(str,kmer_size)
 *
 * Descrip:    makes a lexical kmer from a string
 *
 *
 * Arg:              str [UNKN ] Undocumented argument [char *]
 * Arg:        kmer_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [struct lexical_kmer]
 *
 */
# line 119 "dnamapping.dy"
struct lexical_kmer lexical_dna_number_from_string(char * str,int kmer_size)
{
  int i;
  int no;
  int for_buffer[50];
  int rev_buffer[50];
  kmer_t temp;
  long base = 1;
  
  struct lexical_kmer out;

  assert(kmer_size < 50 && kmer_size > 0);

  for(i=0;i<kmer_size;i++) {
    if( (no = base_from_char(str[i])) == BASE_N ) {
      out.lexical_forward = 3;
      out.kmer = 0;
      return out;
    }
    for_buffer[i] = no;
    rev_buffer[kmer_size-i-1] = complement_base(no);
  }


   
  out.lexical_forward = 1;
  for(i=0;i<kmer_size;i++) {
    if( for_buffer[i] > rev_buffer[i] ) {
      out.lexical_forward = 1;
      break;
    }
    if( for_buffer[i] < rev_buffer[i] ) {
      out.lexical_forward = 0;
      break;
    }
  }

  temp = 0;
  base = 1;

  if( out.lexical_forward == 1 ) {
    for(i=0;i<kmer_size;i++) {
      temp += for_buffer[i]*base;
      base = base * 4;
    }
  } else {
    for(i=0;i<kmer_size;i++) {
      temp += rev_buffer[i]*base;
      base = base * 4;
    }
  }

  out.kmer = temp;

  return out;
}

/* Function:  forward_dna_number_from_string(str,kmer_size)
 *
 * Descrip:    makes a kmer number from a string
 *
 *
 * Arg:              str [UNKN ] Undocumented argument [char *]
 * Arg:        kmer_size [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [kmer_t]
 *
 */
# line 179 "dnamapping.dy"
kmer_t forward_dna_number_from_string(char * str,int kmer_size)
{
  int i;
  long base = 1;
  kmer_t out = 0;

  for(i=0;i<kmer_size;i++) {
    out += base_from_char(toupper(str[i]))*base;
    base = base * 4;

    /* fprintf(stderr,"At %d, got %ld base, %ld out...\n",i,base,out); */
  }

	
  /*  fprintf(stderr," returning for %.*s %ld\n",kmer_size,str,out);*/

  return out;
}

# line 225 "dnamapping.c"

#ifdef _cplusplus
}
#endif
