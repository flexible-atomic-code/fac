#ifndef _PARSER_H_
#define _PARSER_H_ 1

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

char StrTrim(char *s, char c);
int StrSplit(char *s, char sep);
int SetParserQuote(char *qbegin, char *qend);
int SetParserBreak(char *brkch);
int SetParserEscape(char escape);
int SetParserWhite(char *white);
int Parse(char *token, int tokmax, char *line, 
	  int *brkpos, int *next, int *quoted);

#endif
