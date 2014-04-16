/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* msf.c
 * SRE, Sun Jul 11 16:17:32 1993
 * 
 * Export of GCG MSF multiple sequence alignment
 * formatted files.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/* Function: WriteMSF()
 * 
 * Purpose:  Write aseqs, names, weights to an open fp,
 *           in GCG MSF format. The alignment must
 *           be flushed (all aseqs the same length, padded
 *           with gaps)
 * 
 * Return:   (void)
 */
void
WriteMSF(FILE   *fp,            /* open fp for writing           */
	 char  **aseqs,         /* aligned sequences             */
	 AINFO  *ainfo)
{
  int    still_going;		/* True if writing another block */
  int    idx;			/* counter for sequences         */
  int    pos;			/* position counter              */
  int    namelen;		/* maximum name length used      */
  int    len;			/* tmp variable for name lengths */
  char   buffer[51];		/* buffer for writing seq        */
  char **sqptr;                 /* ptrs into each sequence       */
  int    charcount;		/* num. symbols we're writing    */

				/* allocate seq pointers that we'll
				   move across each sequence */
  sqptr = (char **) MallocOrDie (ainfo->nseq * sizeof(char *));

				/* set sqptrs to start of each seq */
  for (idx = 0; idx < ainfo->nseq; idx++)
    sqptr[idx] = aseqs[idx];
				/* calculate max namelen used */
  namelen = 0;
  for (idx = 0; idx < ainfo->nseq; idx++)
    if ((len = strlen(ainfo->sqinfo[idx].name)) > namelen) 
      namelen = len;

  /*****************************************************
   * Write the title line
   *****************************************************/
  fprintf(fp, "\n");
				/* ack! we're writing bullshit here */
  fprintf(fp, "    MSF:  000  Type: X  Check: 0000  ..\n");
  fprintf(fp, "\n");

  /*****************************************************
   * Write the names
   *****************************************************/

  for (idx = 0; idx < ainfo->nseq; idx++)
    {
      fprintf(fp, "  Name: %-*.*s  Len:  %5d  Check:  %5d  Weight: %.4f\n",
	      namelen, namelen,
	      ainfo->sqinfo[idx].name,
	      ainfo->alen,
	      GCGchecksum(aseqs[idx], ainfo->alen),
	      ainfo->wgt[idx]);
    }
  fprintf(fp, "\n");
  fprintf(fp, "//\n");
  fprintf(fp, "\n");

  /*****************************************************
   * Write the sequences
   *****************************************************/

  still_going = 1;
  while (still_going)
    {
      still_going = 0;
      for (idx = 0; idx < ainfo->nseq; idx++)
	{
	  fprintf(fp, "%-*.*s  ", namelen, namelen, 
		  ainfo->sqinfo[idx].name);

				/* get next line's worth of 50 from seq */
	  strncpy(buffer, sqptr[idx], 50);
	  buffer[50] = '\0';
	  charcount = strlen(buffer);

				/* is there still more to go? */
	  if (charcount == 50 && sqptr[idx][50] != '\0')
	    still_going = 1;

				/* shift the seq ptr by a line */
	  sqptr[idx] += charcount;

				/* draw the sequence line */
	  pos = 0; 
	  while (pos < charcount)
	    {
	      if (isgap(buffer[pos])) fputc('.', fp);
	      else fputc(buffer[pos], fp);
	      pos++;
	      if (!(pos % 10)) fputc(' ', fp);
	    }
	  fputc('\n', fp);
	}
				/* put blank line between blocks */
      fputc('\n', fp);
    }

  free(sqptr);
}



