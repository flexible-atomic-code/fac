/*
 *   FAC - Flexible Atomic Code
 *   Copyright (C) 2001-2015 Ming Feng Gu
 * 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 * 
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *   GNU General Public License for more details.
 * 
 *   You should have received a copy of the GNU General Public License
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

static char *rcsid="$Id: parser.c,v 1.9 2006/08/04 07:43:53 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#include "parser.h"

/*
** I found the original version of this string parser (the Parse() function) 
** on the internet.
** I modified it so that it can handel quotes with different beginning and
** ending tokens, and added some utility functions.
** I have since lost the track of where this code came from, although I believe 
** that it was in the public domain. Please let me know if anyone knows otherwise.
*/

#define IN_WHITE 0
#define IN_TOKEN 1
#define IN_QUOTE 2
#define IN_OZONE 3

#define NN 32

static int state;
static int tokpos;
static char white[NN] = " \t";
static char delim[NN] = ",";
static char quote_begin[NN] = "\"'";
static char quote_end[NN]   = "\"'";
static int  quote_open[NN];
static char escape = '\0';

int SetParserQuote(char *qbegin, char *qend) {
  int n1, n2;
  n1 = strlen(qbegin);
  n2 = strlen(qend);
  if (n1 != n2) return -1;
  if (n1 >= NN) return -1;
  strcpy(quote_begin, qbegin);
  strcpy(quote_end, qend);
  return 0;
}

int SetParserWhite(char *w) {
  if (strlen(w) >= NN) return -1;
  strcpy(white, w);
  return 0;
}

int SetParserBreak(char *d) {
  if (strlen(d) >= NN) return -1;
  strcpy(delim, d);
  return 0;
}

int SetParserEscape(char e) {
  escape = e;
  return 0;
}

char StrTrim(char *s, char c) {
  char *p;
  int i;
  char r;
  
  i = strspn(s, " \t");
  p = s+i;
  i = 0;
  while (p[i] && p[i] != c) {
    s[i] = p[i];
    i++;
  }
  if (i == 0) {
    s[i] = '\0';
    r = '\0';
  } else {
    i--;
    r = p[i];
    while (i >= 0 && (s[i] == ' '  || 
		      s[i] == '\t' || 
		      s[i] == '\n' || 
		      s[i] == '\r' ||
		      s[i] == EOF)) {
      i--;
    }
    s[i+1] = '\0';
  }
  return r;
}
 
int QuotedStrSplit(char *s, char sep, char qb, char qe) {
  char *p;
  int ns, qopen;

  ns = 0;
  p = s;
  qopen = 0;
  while (*p == sep) {
    *p = ' ';
    p++;
  }
  while (*p) {
    if (*p == qb) qopen++;
    else if (*p == qe) {
      qopen--;
      if (qopen < 0) {
	printf("The quoted string %s does not have matched quotes\n", s);
	return -1;
      }
    }
    if (qopen > 0) {
      p++;
      continue;
    }
    if (*p == sep) {
      ns++;
      *p = '\0';
      p++;
      while (*p == sep) {	
	*p = ' ';
	p++;
      }
    } else {
      p++;
    }
  }
  
  if (qopen != 0) {
    printf("The quoted string %s does not have matched quotes\n", s);
    return -1;
  }

  if (p == s) return 0;
  p--;
  while (*p) {
    if (*p != sep) {
      ns++;
      break;
    }
    p--;
  }

  return ns;
}

int StrReplace(int n, char *s, char r0, char r1, char r2, char r3) {
  int i = 0;
  int m = 0;
  char *c;
  c = s;
  while (i < n && *c) {
    if (*c == r0) {      
      if (i == n-1 || (i < n-1 && c[1] == '\0')) *c = r3;
      else if (m == 0) *c = r1;
      else *c = r2;
      m++;
    }
    i++;
    c++;
  }
  return m;
}

int StrSplit(char *s, char sep) {
  char *p;
  int ns;

  ns = 0;
  p = s;
  while (*p == sep) {
    *p = ' ';
    p++;
  }
  while (*p) {
    if (*p == sep) {
      ns++;
      *p = '\0';
      p++;
      while (*p == sep) {
	*p = ' ';
	p++;
      }
    } else {
      p++;
    }
  }
  if (p == s) return 0;
  p--;
  while (*p) {
    if (*p != sep) {
      ns++;
      break;
    }
    p--;
  }

  return ns;
}

static int sindex(char ch, char *s) {
  char *cp;

  for (cp = s; *cp; cp++) {
    if (*cp == ch) {
      return (int)(cp-s);
    }
  }
  return -1;
}

static void chstore(char *s, int max, char ch) {
  if (tokpos >= 0 && tokpos < max-1) {
    s[tokpos++] = ch;
  }
}

int Parse(char *token, int tokmax, char *line, 
	  int *next, int *brkpos, int *quotepos) {
  int qp;
  int i, c, nc;

  *brkpos = -1;
  *quotepos = -1;
  for (i = 0; i < NN; i++) quote_open[i] = 0;
  
  if (!line[*next]) return 1;

  state = IN_WHITE;
  
  for (tokpos = 0; (c = line[*next]); (*next)++) {
    if ((qp = sindex(c, delim)) >= 0) {
      switch(state) {
      case IN_WHITE:
      case IN_TOKEN:
      case IN_OZONE:
	(*next)++;
	*brkpos = qp;
	goto OUT;

      case IN_QUOTE:
	chstore(token, tokmax, c);
	break;
      }
    } else if ((qp = sindex(c, quote_begin)) >= 0) {
      switch (state) {
      case IN_WHITE:
	state = IN_QUOTE;
	*quotepos = qp;
	quote_open[qp]++;
	break;

      case IN_QUOTE:
	if (qp == *quotepos) {
	  if (c == quote_end[*quotepos]) {
	    state = IN_OZONE;
	    quote_open[qp]--;
	  } else {
	    quote_open[qp]++;
	    chstore(token, tokmax, c);
	  }
	} else {
	  chstore(token, tokmax, c);
	}
	break;
      case IN_TOKEN:
      case IN_OZONE:
	goto OUT;
      }
    } else if ((qp = sindex(c, quote_end)) >= 0) {
      switch (state) {
      case IN_QUOTE:
	if (quote_open[qp] > 1) {
	  chstore(token, tokmax, c);
	  quote_open[qp]--;
	} else if (quote_open[qp] < 1) {
	  chstore(token, tokmax, c);
	} else {
	  state = IN_OZONE;
	  quote_open[qp]--;
	}
	break;
      default:
	*quotepos = qp;
	return -1;
      }
    } else if ((qp = sindex(c, white)) >= 0) {
      switch (state) {
      case IN_WHITE:
      case IN_OZONE:
	break;
      case IN_TOKEN:
	state = IN_OZONE;
	break;
      case IN_QUOTE:
	chstore(token, tokmax, c);
	break;
      }
    } else if (c == escape) {
      nc = line[(*next)+1];
      if (nc == '\0') {
	chstore(token, tokmax, c);
	(*next)++;
	goto OUT;
      }
      switch(state) {
      case IN_WHITE:
	(*next)--;
	state = IN_TOKEN;
	break;
      case IN_TOKEN:
      case IN_QUOTE:
	(*next)++;
	chstore(token, tokmax, nc);
	break;
      case IN_OZONE:
	goto OUT;
      }
    } else {
      switch (state) {
      case IN_WHITE:
	state = IN_TOKEN;
      case IN_TOKEN:
      case IN_QUOTE:
	chstore(token, tokmax, c);
	break;
      case IN_OZONE:
	goto OUT;
      }
    }
  }

 OUT:
  token[tokpos] = '\0';
  
  return 0;
}


