#ifndef _SFAC_H_
#define _SFAC_H_ 1

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "init.h"

#define MAXNARGS 128
#define MAXLINELENGTH 1024
#define MAXMETHODNAME 64

#define ERR_LINEUNTERMINATED (-1)
#define ERR_LINETOOLONG      (-2)
#define ERR_NOVARIABLE       (-3)
#define ERR_ARGSTOOMANY      (-4)
#define ERR_SYNTAX           (-5)
#define ERR_EVAL             (-6)

#define NUMBER  (0)
#define STRING  (1)
#define TUPLE   (2)
#define LIST    (3)
#define KEYWORD (4)

typedef struct _METHOD_ {
  char name[MAXMETHODNAME];
  int (*func)(int argc, char *argv[], int argt[], ARRAY *variables);
  int arg_style;
} METHOD;

typedef struct _VARIABLE_ {
  char *name;
  char *value;
  int type;
} VARIABLE;

typedef struct _STATEMENT_ {
  int nline;
  int imethod;
  int argc;
  char *argv[MAXNARGS];
  int argt[MAXNARGS];
} STATEMENT;

int GetValidLine(FILE *f, char *line, int *nlines);
int TokenizeLine(int nline, char *line, ARRAY *statements, ARRAY *variables);
int EvalFile(FILE *f, int exebyline);
int EvalStatement(STATEMENT *st, ARRAY *variables);
void FreeStatementData(void *p);
void FreeVariableData(void *p);
void ErrorOcurred(int ierr, int loc);

#endif

