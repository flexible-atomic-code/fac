#include "dbase.h"

static char *rcsid="$Id: dbase.c,v 1.4 2001/12/14 21:35:32 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif


static F_HEADER fheader;
static EN_HEADER en_header;
static TR_HEADER tr_header;

int InitDBase() {
  fheader.tsession = time(0);
  fheader.version = VERSION;
  fheader.sversion = SUBVERSION;
  fheader.ssversion = SUBSUBVERSION;
  fheader.type = 0;
  fheader.atom = 0;
  fheader.nblocks = 0;
}

FILE *InitFile(char *fn, F_HEADER *fhdr, void *rhdr) {
  EN_HEADER *en_hdr;
  TR_HEADER *tr_hdr;
  FILE *f;
  off_t p;
  size_t n;

  f = fopen(fn, "r+");
  if (f == NULL) {
    f = fopen(fn, "w");
    if (f == NULL) return NULL;
  }

  fheader.type = fhdr->type;
  strcpy(fheader.symbol, fhdr->symbol);
  fheader.atom = fhdr->atom;
  fheader.nblocks += 1;
  n = fwrite(&fheader, sizeof(F_HEADER), 1, f);

  fseek(f, 0, SEEK_END);
  p = ftell(f);

  switch (fhdr->type) {
  case DB_EN:
    en_hdr = (EN_HEADER *) rhdr;
    en_header.position = p;
    en_header.length = sizeof(EN_HEADER);
    en_header.nele = en_hdr->nele;
    en_header.nlevels = 0;
    n = fwrite(&en_header, sizeof(EN_HEADER), 1, f);
    break;
  case DB_TR:
    tr_hdr = (TR_HEADER *) rhdr;
    tr_header.position = p;
    tr_header.length = sizeof(TR_HEADER);
    tr_header.nele = tr_hdr->nele;
    tr_header.gauge = tr_hdr->gauge;
    tr_header.mode = tr_hdr->mode;
    tr_header.multipole = tr_hdr->multipole;
    tr_header.ntransitions = 0;
    n = fwrite(&tr_header, sizeof(TR_HEADER), 1, f);
    break;
  case DB_CE:
    break;
  case DB_RR:
    break;
  case DB_AI:
    break;
  case DB_CI:
    break;
  default:
    break;
  }

  return f;
}

int CloseFile(FILE *f) {
  int n;

  if (f == NULL || fheader.type <= 0) return 0;

  switch (fheader.type) {
  case DB_EN:
    fseek(f, en_header.position, SEEK_SET);  
    n = fwrite(&en_header, sizeof(EN_HEADER), 1, f);
    break;
  case DB_TR:
    fseek(f, tr_header.position, SEEK_SET);
    n = fwrite(&tr_header, sizeof(TR_HEADER), 1, f);
    break;
  case DB_CE:
    break;
  case DB_RR:
    break;
  case DB_AI:
    break;
  case DB_CI:
    break;
  default:
    break;
  }

  fheader.type = 0;
  fclose(f);

  return 0;
}

int PrintTable(char *ifn, char *ofn) {
  F_HEADER fh;
  FILE *f1, *f2;
  int n;
  
  f1 = fopen(ifn, "r");
  if (f1 == NULL) return -1;

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }

  n = fread(&fh, sizeof(F_HEADER), 1, f1);
  if (n != 1) return 0;
  fprintf(f2, "FAC %d.%d.%d\n", fh.version, fh.sversion, fh.ssversion);
  fprintf(f2, "TSess\t= %uL\n", fh.tsession);
  fprintf(f2, "Type\t= %d\n", fh.type);
  fprintf(f2, "%s Z\t= %d\n", fh.symbol, fh.atom);
  fprintf(f2, "NBlocks\t= %d\n", fh.nblocks);
  
  switch (fh.type) {
  case DB_EN:
    n = PrintENTable(f1, f2);
    break;
  case DB_TR:
    n = PrintTRTable(f1, f2);
    break;
  case DB_CE:
    break;
  case DB_RR:
    break;
  case DB_AI:
    break;
  case DB_CI:
    break;
  default:
    break;
  }
  
  fclose(f1);
  if (f2 != stdout) fclose(f2);
  return n;
}

int WriteENRecord(FILE *f, EN_RECORD *r) {
  int n;
  en_header.nlevels += 1;
  en_header.length += sizeof(EN_RECORD);
  n = fwrite(r, sizeof(EN_RECORD), 1, f);  
  return n;
}

int PrintENTable(FILE *f1, FILE *f2) {
  EN_HEADER h;
  EN_RECORD r;
  int n, i;
  int nb;
  float e;

  nb = 0;
  while (1) {
    n = fread(&h, sizeof(EN_HEADER), 1, f1);
    if (n != 1) break;
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NLEV\t= %d\n", h.nlevels);
    fprintf(f2, "         Energy(eV)   P 2J \n");
    for (i = 0; i < h.nlevels; i++) {
      n = fread(&r, sizeof(EN_RECORD), 1, f1);
      if (n != 1) break;
      e = r.energy;
      e *= HARTREE_EV;
      fprintf(f2, "%5d %15.8E %1d %2d %-20s %-20s %-s\n",
	      r.ilev, e, r.p, r.j, r.ncomplex, r.sname, r.name);
    }
    nb += 1;
  }
  
  return nb;
}

int WriteTRRecord(FILE *f, TR_RECORD *r) {
  int n;
  tr_header.ntransitions += 1;
  tr_header.length += sizeof(TR_RECORD);
  n = fwrite(r, sizeof(TR_RECORD), 1, f);  
  return n;
}

int PrintTRTable(FILE *f1, FILE *f2) {
  TR_HEADER h;
  TR_RECORD r;
  int n, i;
  int nb;

  nb = 0;
  
  while (1) {
    n = fread(&h, sizeof(TR_HEADER), 1, f1);
    if (n != 1) break;
    
    fprintf(f2, "\n");
    fprintf(f2, "NELE\t= %d\n", h.nele);
    fprintf(f2, "NTRS\t= %d\n", h.ntransitions);
    fprintf(f2, "Multi\t= %d\n", (int)h.multipole);
    fprintf(f2, "Gauge\t= %d\n", (int)h.gauge);
    fprintf(f2, "Mode\t= %d\n", (int)h.mode);
    fprintf(f2, "Upper\tLower\t    gf\n");
    
    for (i = 0; i < h.ntransitions; i++) {
      n = fread(&r, sizeof(TR_RECORD), 1, f1);
      if (n != 1) break;
      fprintf(f2, "%5d\t%5d\t%15.8E\n", 
	      r.upper, r.lower, r.strength);
    }
    nb += 1;
  }

  return nb;
}

	    
