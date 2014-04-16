
#include "sequence.h"
#include "dnamapping.h"




static kmer_t ** number_list;

static int memory_point;
static int active_row;
static int cursor;

static int buffer_length = 16777216;
static int kmer_size = 21;

int add_number(kmer_t kmer)
{
  if( cursor > buffer_length) {
    if( active_row > memory_point ) {
      memory_point *= 2;
      number_list = realloc(number_list,sizeof(kmer_t **)*memory_point);
      assert(number_list != NULL);
    }
    active_row++;
    cursor = 0;
  }

  number_list[active_row][cursor] = kmer;
}


void print_buffer_debug(void)
{
  int i;
  int j;

  for(i=0;i<active_row;i++) {

  }
  

}


void initialise_buffer(void)
{
  int i;

  memory_point = 512;
  active_row = 0;
  cursor = 0;

  number_list = calloc(memory_point,sizeof(kmer_t *));
  
  for(i=0;i<memory_point;i++) {
    number_list[i] = calloc(buffer_length,sizeof(kmer_t));
    assert(number_list[i] != NULL);
  }
}
  

int main(int argc,char ** argv)
{
  long long base_buffer[400];
  int i;
  int j;
  kmer_t kmer;
  long long base = 1;
  int skip;
  Sequence * seq;

  initialise_buffer();

  for(i=0;i<kmer_size;i++) {
    base_buffer[i] = base;
    base = base * 4;
  }

  while( (seq = read_fasta_Sequence(stdin)) != NULL ) {
    for(i=0;i<seq->len-kmer_size;i++) {
      kmer = 0;
      skip = 0;
      for(j=0;j<kmer_size;j++) {
	if( toupper(seq->seq[j]) == 'N' ) {
	  skip = 1;
	  break;
	}
	kmer += base_from_char(toupper(seq->seq[j])) * base_buffer[j];
      }

      if( skip == 0 ) {
	add_number(kmer);
      }
    }
  }


}




