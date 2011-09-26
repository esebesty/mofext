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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"

struct sym_matrix mymatrix;

int load_matrix(char *file_name) {

  FILE *matrix_file;
  char line[100];
  int  i,j;/*Loop variables*/
  char *tok;
  int   item_num; /* The items count in the matrix */

if ((matrix_file = fopen(file_name,"r"))==NULL){
	fprintf(stderr,"I can not open the matrix file\n");
	return(-1);
}

/*Read the header*/
fgets(line,100,matrix_file);
item_num = atoi(strtok(line," "));

for (i=0; i<item_num; i++){
	tok=strtok(NULL," ");
	if (tok == NULL){
		fprintf(stderr,"Missing character from the matrix header\n");
		exit(4);
	}
	strncat(mymatrix.header,tok,1);
}

/*while( (tok=strtok(NULL," "))!=NULL){
	strncat( mymatrix.header,tok,1);
}*/

/*Cut the newline*/
mymatrix.header[strlen(mymatrix.header)-1]='\0';

/*Fill the matrix*/
for(i=0;i<item_num;i++){ 

	fgets(line,100,matrix_file);
	line[strlen(line)-1]='\0';
	strtok(line," ");
	for(j=0;j<=i;j++){
		tok=strtok(NULL," ");
		if (tok == NULL) {
			fprintf(stderr,"Error in matrix %d. column %d. row\n",j,i);
			exit(4);
		}
		mymatrix.score[i][j]=atoi(tok);
		if(i!=j) mymatrix.score[j][i]=atoi(tok);
	}
}
mymatrix.item_num = item_num;
fclose(matrix_file);
return(0);
}

double getscore(char A, char B){

  int  posA = -1, posB = -1;
  int  i;

for(i=0; i < mymatrix.item_num; i++){
	if (A == mymatrix.header[i]) posA=i;
	if (B == mymatrix.header[i]) posB=i;
	if ((posA != -1) && (posB != -1)) break;
}

return(mymatrix.score[posA][posB]);
}

int validseq(char *seq){
unsigned int i,lenseq;
int          j,matchpos;

lenseq = strlen(seq);

for (i=0; i<lenseq; i++){
   matchpos = -1;
   for (j=0; j < mymatrix.item_num; j++)
      if (seq[i] == mymatrix.header[j]) matchpos = i;
   if (matchpos == -1) return(i+1);  /* Return the bad character position */

}

return(-1);
}


/*End of module*/
