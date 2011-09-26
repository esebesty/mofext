/*
MoFext Motif Find and extract
Copyright (C) 2006 Tibor Nagy

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*  This module contain a lots of useful string rutines */
#include "mystrings.h"
#include <stdlib.h>



int substr(char *instr,unsigned int start,unsigned int length, char *retstr){
  unsigned int     i; /*Loop variable*/
  char            *startstr;

/*Error check*/
if ( (start+length) > strlen(instr)) return(-1);
for (i=0;i<start;i++) startstr = instr++;

strncpy(retstr,instr,length);
/* Some times if this line is missing, the retstr is not null terminated*/
retstr[length] = '\0';
return(0);
}/* substr */

/*Suprise !!!*/
