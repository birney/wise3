#ifndef DYNAMITEdnamappingHEADERFILE
#define DYNAMITEdnamappingHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "codon.h"

/*typedef long int kmer_t;*/
typedef long long kmer_t;

struct lexical_kmer {
  kmer_t kmer;
  char lexical_forward;
};



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



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
char * Wise2_map_to_basepair_numbers(char * seq,long int len);
#define map_to_basepair_numbers Wise2_map_to_basepair_numbers


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
kmer_t Wise2_reverse_complement_dna_number(kmer_t number,int kmer_size);
#define reverse_complement_dna_number Wise2_reverse_complement_dna_number


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
void Wise2_reverse_map_dna_number(kmer_t number,int kmer_size,char * buffer);
#define reverse_map_dna_number Wise2_reverse_map_dna_number


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
char Wise2_lexical_last_base_from_kmer(kmer_t kmer,int kmer_size,int lexical_for);
#define lexical_last_base_from_kmer Wise2_lexical_last_base_from_kmer


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
struct lexical_kmer Wise2_lexical_dna_number_from_string(char * str,int kmer_size);
#define lexical_dna_number_from_string Wise2_lexical_dna_number_from_string


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
kmer_t Wise2_forward_dna_number_from_string(char * str,int kmer_size);
#define forward_dna_number_from_string Wise2_forward_dna_number_from_string


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/

#ifdef _cplusplus
}
#endif

#endif
