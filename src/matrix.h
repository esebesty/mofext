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
struct sym_matrix{
	int    item_num;
	double score[255][255];
	char   header[255];
};

extern struct sym_matrix mymatrix;


extern int         load_matrix(char *file_name);
extern double      getscore(char A, char B);
int                validseq(char *seq);
