/* overprint.c
   A subroutine to overprint screen output of a character string.
   Take any string and number combination in a Fortran program,
   make an internal write to a 60 character-long buffer, and 
   send it to overprint(buf). This is very useful for calling
   from Fortran programs in order to keep track of model progress.
   Jed O. Kaplan, June, 2004
*/
 
#include <stdio.h>

#ifdef ifort
void overprint_(char *text)
#else 
void overprint(char *text)
#endif
{
	fprintf(stderr,"%s\r",text);
	fflush(stderr);
}
