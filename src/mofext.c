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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"
#include "similar.h"
#include "getline.h"


#define LINESIZE 100
#define BASIC_HELP 1
#define FULL_HELP 3
#define FEATURE_NUM 50
#define TRUE 1
#define FALSE 0
#define OUTFORMATLEN 12


void showhelp(int which){
   if (which > 0) {
      printf("        **** Mofext v1.0.3 ****\n");
      printf("Usage: mofext -d mypatterns1.list mypatterns2.list -q GGATCC TTGANTGA -m matrix -w 10 -c 0.95\n\n");
      printf("Options:\n");
      printf("~~~~~~~~\n");
      printf("-h     Display full help.\n");
      printf("-d     Databases to search. Space separated, maximum 50.\n");
      printf("-q     Query patterns. Space separated, maximum 50.\n");
      printf("-m     The similarity matrix.\n");
      printf("-w     Wordsize. Default: 6.\n");
      printf("-c     The similarity percentage limit (cutoff). Default: 0.8.\n");
      printf("-o     Output format. See below.\n");
      printf("Example:\n");
      printf("mofext -d mypatterns1.list mypatterns2.list -q GGATCC TTGANTGA -m matrix -w 4 -c 0.5\n\n");
   }

   if (which > BASIC_HELP) {
      printf("Matrix:\n");
      printf("~~~~~~~\n");
      printf(" The matrix used by mofext is a plain text file. Each row ends with a newline character.\n");
      printf(" Columns are separated by one or more spaces. The matrix is quadratic, so A[j,i] = A[i,j].\n");
      printf(" The number in the upper left corner denotes the matrix size. See default matrix included\n");
      printf(" in the source.\n");
   }

   if (which > 2){
      printf("Motif list:\n");
      printf("~~~~~~~~~~~\n");
      printf(" The motif list file used by mofext is a plain text file. Each row ends with a newline character.\n");
      printf(" The motif list file must have at least 2 columns, separated by spaces or tabs. The first column\n");
      printf(" is a UNIQUE motif id. It can contain any characters except space, newline or tab. The second\n");
      printf(" column is the sequence (pattern) of the motif. The available characters are the ones defined in\n");
      printf(" the matrix file. Additional columns are ignored by the program.\n");
      printf("Output:\n");
      printf("~~~~~~~\n");
      printf(" The -o option sets the output format of the program. The following characters represent an output\n");
      printf(" element. The default is ieqd.\n");
      printf(" i: ID column.\n");
      printf(" s: Score.\n");
      printf(" e: Extended score.\n");
      printf(" p: Probability value.\n"); /*It is an undocumented feature */
      printf(" d: Hit subsequence.\n");
      printf(" D: Full hit sequence.\n");
      printf(" q: Query subsequence.\n");
      printf(" Q: Full query sequence.\n");
      printf(" F: The position of the first base of the hit subsequence in the full hit sequence.\n");
      printf(" f: The position of the first base of the query subsequence in the full query sequence.\n");
      printf(" %%: The similarity percent score of the hit/query subsequence pair.\n");
      printf(" Example: mofext -d mypatterns1.list -q GGATCC -m matrix -o ied\n");
      printf(" The program prints out the ID, extended score and the hit subsequence.\n");
    }
    exit(0);
}

int main(int argc, char *argv[]){
  int               i, j = -1;        /*Loop variables*/
  char             *list_array[FEATURE_NUM];  /*The database filenames*/
  char             *patt_array[FEATURE_NUM];  /*The pattern seqs*/
  int               c_la = 0;        /*The list array count*/
  int               c_pa = 0;        /*The patt array count*/
  char             *line;            /*The line*/
  size_t            ls;              /*The line size*/
  char             *id;              /*The sequence's id*/
  char             *seq;             /*The sequence*/
  unsigned int      ws = 6;          /* Word size */
  double            limit=0.9;       /* The print score percentage limit */
  FILE             *list_file;
  double            bases[4];        /* Frequences of bases */
  char              list=0;
  char              pattern=0;
  char              output[OUTFORMATLEN] = {'i','e','q','d'};      /* The output format */
  int               mofext_error=FALSE;  /* The error presence*/
  int               error_pos;


/*The argument handling*/
if (argc==1) showhelp(BASIC_HELP);

for (i=1;i<argc;i++){
	if (argv[i][0] == '-') {
                switch(argv[i][1]){
                   case 'd':
                      list    = TRUE;
                      pattern = FALSE;
                      break;
                   case 'h':
                      showhelp(FULL_HELP);
                   case 'q':
                      pattern = TRUE;
                      list    = FALSE;
                      break;
                   case 'f':
                      c_pa = -1;
                      pattern = FALSE;
                      list    = FALSE;
                      break;
                   case 'm':
                      j       = load_matrix(argv[i+1]);
                      pattern = FALSE;
                      list    = FALSE;
                      break;
                   case 'w':
                      ws      = atoi(argv[i+1]);
                      if(ws < 5){
                         fprintf(stderr, "Wordsize too small\n");
                         mofext_error = TRUE;
                      }
                      pattern = FALSE;
                      list    = FALSE;
                      break;
                   case 'c':
                      limit = atof(argv[i+1]);
                      if( (limit <= 0.0) || (limit > 1.0)){
                         fprintf(stderr, "Cutoff parameter should be between 0.0 - 1.0\n");
                         mofext_error = TRUE;
                      }
                      pattern = FALSE;
                      list    = FALSE;
                      break;
                   case 'o':
                      strncpy(output, argv[i+1], OUTFORMATLEN - 1);
                      pattern = FALSE;
                      list    = FALSE;
                      break;
                }
	}
	else {
		if (list) {
			if (c_la == FEATURE_NUM - 1){
				fprintf(stderr,"I can't handle more database\n");
				continue;
			}
			list_array[c_la] = malloc(strlen(argv[i])+1);
			strcpy(list_array[c_la],argv[i]);
			c_la++;
		}
		if (pattern) {
			if (c_pa == FEATURE_NUM - 1){
				fprintf(stderr,"I can't handle more pattern\n");
				continue;
			}
			patt_array[c_pa] = malloc(strlen(argv[i])+1);
			strcpy(patt_array[c_pa],argv[i]);
			c_pa++;
		}
	}
}
/* Check the presence of the matrix argument */
if (j == -1) {
	fprintf(stderr,"Can't open matrix file\n");
	mofext_error = TRUE;
}
/* Walk through all the patterns and all the files */
line = malloc(LINESIZE);
if (line == NULL) {
	fprintf(stderr,"Memory not enought\n");
	exit(2);
}

if (c_pa == 0) {
	fprintf(stderr,"Query pattern did not specified\n");
	mofext_error = TRUE;
}
if (c_la == 0){
	fprintf(stderr,"Motif database did not specified\n");
	mofext_error = TRUE;
}

if (mofext_error == TRUE) showhelp(BASIC_HELP);


ls   = LINESIZE;
for (j=0; j<c_pa;j++){

  if (strlen(patt_array[j]) < ws) {
         fprintf(stderr,"%s pattern smaller than the wordsize\n",patt_array[j]);
         continue;
  }
  if ( (error_pos=validseq(patt_array[j])) != -1) {
         fprintf(stderr,"Invalid character in the query (%d. position)\n",error_pos);
         continue;
  }
   for (i=0;i<c_la;i++){
	/* Open the file and check */
	list_file = fopen(list_array[i],"r");

	/* process the lines */
	while( getline(&line,&ls,list_file) != -1 ){
	   id   = strtok(line," \t");
	   seq  = strtok(NULL," \t\n");

           if (seq == NULL) {
              fprintf(stderr,"Invalid row in file %s\n",list_array[i]);
              continue;
           }

           basefreq(seq,bases);

           if ((error_pos=validseq(seq))!= -1){
              fprintf(stderr,"Invalid character in the database element:%s (%d. position)\n",seq,error_pos);
              continue;
           }

	   /* Find the similarity and try to extend */
   	   similarity(id,seq,patt_array[j],ws,limit,bases,output);
	}/* while ( getline )*/

	fclose(list_file);
   }
}



/*Ende*/
return(0);
}
