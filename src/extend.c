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
*/
#include <string.h>
#include <stdio.h>
/*#include <stdlib.h>*/
#include "extend.h"
#include "matrix.h"

void GetCorrectPos(int *sA, int *sB, int *lA, int *lB){

int startA,startB,lenA;

if((sA == NULL) || (sB == NULL) || 
   (lA == NULL) || (lB == NULL)) return;

startA = (*sA) - (*sB);
if (startA < 0) startA = 0;
startB = (*sB) - (*sA);
if (startB < 0) startB = 0;

lenA = *lA - startA;
*lB  = *lB - startB;

if (lenA > *lB) lenA = *lB;

*sA = startA;
*sB = startB;
*lA = lenA;
*lB = lenA;
return;
}

int ExtScore(char *seqA, char *seqB, int ws, int *newStart, int *newLen){
struct max {
	int    sc;
	int    st;
	int    ln;
};

struct max m;
int    i, j, k, seqAl;
int score = 0;

seqAl = strlen(seqA);

m.sc = 0;
m.st = 0;
m.ln = 0;

for (i=0; i <= (seqAl - ws); i++){
	for (j=ws; j<=(seqAl-i); j++){
		score = 0;
		for (k=i; k<(i+j); k++)
			score += getscore(seqA[k],seqB[k]);

		if (score > m.sc) {
			m.st = i;
			m.ln = j;
			m.sc = score;
		}
	}

}

*newStart = m.st;
*newLen   = m.ln;
return(m.sc);
}














