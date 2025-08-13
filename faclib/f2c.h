extern double f_erf(double);
#define ERF(a0) f_erf((a0))
extern double f_fm5p(double);
#define FM5P(a0) f_fm5p((a0))
extern double f_fm5pi(double);
#define FM5PI(a0) f_fm5pi((a0))
extern double f_fm3p(double);
#define FM3P(a0) f_fm3p((a0))
extern double f_fm3pi(double);
#define FM3PI(a0) f_fm3pi((a0))
extern double f_fm1p(double);
#define FM1P(a0) f_fm1p((a0))
extern double f_fm1pi(double);
#define FM1PI(a0) f_fm1pi((a0))
extern double f_fm1m(double);
#define FM1M(a0) f_fm1m((a0))
extern double f_fm1mi(double);
#define FM1MI(a0) f_fm1mi((a0))
extern double f_expint(double,int);
#define EXPINT(a0,a1) f_expint((a0),(a1))
extern double f_eione(int,double);
#define EIONE(a0,a1) f_eione((a0),(a1))
extern double f_d1mach(int);
#define D1MACH(a0) f_d1mach((a0))
extern void f_radfnd(int,double *);
#define RADFND(a0,a1) f_radfnd((a0),(a1))
extern void f_locsep(int,int,double *,double *);
#define LOCSEP(a0,a1,a2,a3) f_locsep((a0),(a1),(a2),(a3))
extern void f_rgmqed(double *,double *);
#define RGMQED(a0,a1) f_rgmqed((a0),(a1))
extern void f_fsedat(int,int,int,int,double *,double *);
#define FSEDAT(a0,a1,a2,a3,a4,a5) f_fsedat((a0),(a1),(a2),(a3),(a4),(a5))
extern void f_genqed(int,int,int,int,int,double *,double *,double *,double *,double *,double *);
#define GENQED(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) f_genqed((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10))
extern void f_iniqed(double,int,int,double);
#define INIQED(a0,a1,a2,a3) f_iniqed((a0),(a1),(a2),(a3))
extern void f_mohrfin(int,int,double,double,double *,double *,double *,double *,double *);
#define MOHRFIN(a0,a1,a2,a3,a4,a5,a6,a7,a8) f_mohrfin((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8))
extern void f_uvip3p(int,int,double *,double *,int,double *,double *);
#define UVIP3P(a0,a1,a2,a3,a4,a5,a6) f_uvip3p((a0),(a1),(a2),(a3),(a4),(a5),(a6))
extern void f_uvip3i(int,int,double *,double *,int,double *,double *,double *);
#define UVIP3I(a0,a1,a2,a3,a4,a5,a6,a7) f_uvip3i((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7))
extern void f_uvip3c(int,int,double *,double *,double *,double *,double *);
#define UVIP3C(a0,a1,a2,a3,a4,a5,a6) f_uvip3c((a0),(a1),(a2),(a3),(a4),(a5),(a6))
extern void f_sdbi3p(int,int,double *,double *,double *,int,double *,double *,double *,int *,double *,int *);
#define SDBI3P(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11) f_sdbi3p((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10),(a11))
extern void f_rgbi3p(int,int,int,double *,double *,double *,int,double *,double *,double *,int *,double *);
#define RGBI3P(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11) f_rgbi3p((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10),(a11))
extern void f_lmqn(int *,int,double *,double *,double *,double *,int,void (*)(int *,double *,double *,double *),int,int,int,double,double,double,double);
#define LMQN(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14) f_lmqn((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10),(a11),(a12),(a13),(a14))
extern void f_lmqnbc(int *,int,double *,double *,double *,double *,int,void (*)(int *,double *,double *,double *),double *,double *,int *,int,int,int,double,double,double,double);
#define LMQNBC(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17) f_lmqnbc((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10),(a11),(a12),(a13),(a14),(a15),(a16),(a17))
extern void f_subplx(double (*)(int *, double *),int,double,int,int,double *,double *,double *,int *,double *,int *,int *);
#define SUBPLX(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11) f_subplx((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10),(a11))
extern double f_argam(double,double);
#define ARGAM(a0,a1) f_argam((a0),(a1))
extern double f_dlogam(double);
#define DLOGAM(a0) f_dlogam((a0))
extern void f_beslik(double,double,double *,double *);
#define BESLIK(a0,a1,a2,a3) f_beslik((a0),(a1),(a2),(a3))
extern double f_besljn(int,int,double);
#define BESLJN(a0,a1,a2) f_besljn((a0),(a1),(a2))
extern void f_y5n(double,double,double,double,double *,double *,double *,double *,int *);
#define Y5N(a0,a1,a2,a3,a4,a5,a6,a7,a8) f_y5n((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8))
extern void f_dcoul(double,double,int,double,double *,double *,double *,double *,int *);
#define DCOUL(a0,a1,a2,a3,a4,a5,a6,a7,a8) f_dcoul((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8))
extern void f_dcoul1(double,double,int,double,double *,double *,double *,double *,int *);
#define DCOUL1(a0,a1,a2,a3,a4,a5,a6,a7,a8) f_dcoul1((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8))
extern void f_cmultip(double,double,double,double,double,int,int,int,double *,int,int *);
#define CMULTIP(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) f_cmultip((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10))
extern void f_acofz1(double,double,int,int,double *,double *,int,int);
#define ACOFZ1(a0,a1,a2,a3,a4,a5,a6,a7) f_acofz1((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7))
extern void f_pixz1(double,double,int,int,double *,double *,double *,double *,int,int,int,int);
#define PIXZ1(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11) f_pixz1((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10),(a11))
extern void f_lmder(void (*)(int *,int *,double *,double *,double *,int *,int *),int,int,double *,double *,double *,int,double,double,double,int,double *,int,double,int,int *,int *,int *,int *,double *,double *,double *,double *,double *);
#define LMDER(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23) f_lmder((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10),(a11),(a12),(a13),(a14),(a15),(a16),(a17),(a18),(a19),(a20),(a21),(a22),(a23))
extern void f_chkder(int,int,double *,double *,double *,int,double *,double *,int,double *);
#define CHKDER(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9) f_chkder((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9))
extern double f_ddot(int,double *,int,double *,int);
#define DDOT(a0,a1,a2,a3,a4) f_ddot((a0),(a1),(a2),(a3),(a4))
extern void f_dscal(int,double,double *,int);
#define DSCAL(a0,a1,a2,a3) f_dscal((a0),(a1),(a2),(a3))
extern void f_dgemv(char *,int,int,double,double *,int,double *,int,double,double *,int);
#define DGEMV(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) f_dgemv((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10))
extern void f_dsbev(char *,char *,int,int,double *,int,double *,double *,int,double *,int *);
#define DSBEV(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) f_dsbev((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10))
extern void f_dspevd(char *,char *,int,double *,double *,double *,int,double *,int,int *,int,int *);
#define DSPEVD(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11) f_dspevd((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10),(a11))
extern void f_dspev(char *,char *,int,double *,double *,double *,int,double *,int *);
#define DSPEV(a0,a1,a2,a3,a4,a5,a6,a7,a8) f_dspev((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8))
extern void f_dgeev(char *,char *,int,double *,int,double *,double *,double *,int,double *,int,double *,int,int *);
#define DGEEV(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13) f_dgeev((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10),(a11),(a12),(a13))
extern void f_dgesv(int,int,double *,int,int *,double *,int,int *);
#define DGESV(a0,a1,a2,a3,a4,a5,a6,a7) f_dgesv((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7))
extern void f_dgesdd(char *,int,int,double *,int,double *,double *,int,double *,int,double *,int,int *,int *);
#define DGESDD(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13) f_dgesdd((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10),(a11),(a12),(a13))
extern void f_dgbsv(int,int,int,int,double *,int,int *,double *,int,int *);
#define DGBSV(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9) f_dgbsv((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9))
extern void f_lsode(void (*)(int *,double *,double *,double *),int *,double *,double *,double,int,double,double *,int,int *,int,double *,int,int *,int,void (*)(int *,double *,double *,int *,int *,double *,int *),int);
#define LSODE(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16) f_lsode((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10),(a11),(a12),(a13),(a14),(a15),(a16))
extern void f_dqagi(double (*)(double *),double,int,double,double,double *,double *,int *,int *,int,int,int *,int *,double *);
#define DQAGI(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13) f_dqagi((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10),(a11),(a12),(a13))
extern void f_dqng(double (*)(double *),double,double,double,double,double *,double *,int *,int *);
#define DQNG(a0,a1,a2,a3,a4,a5,a6,a7,a8) f_dqng((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8))
extern void f_dqags(double (*)(double *),double,double,double,double,double *,double *,int *,int *,int,int,int *,int *,double *);
#define DQAGS(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13) f_dqags((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8),(a9),(a10),(a11),(a12),(a13))
extern void f_fermid(double,double,double,double *,int *);
#define FERMID(a0,a1,a2,a3,a4) f_fermid((a0),(a1),(a2),(a3),(a4))
extern void f_ferinc(double,double,double,double,double *,int *);
#define FERINC(a0,a1,a2,a3,a4,a5) f_ferinc((a0),(a1),(a2),(a3),(a4),(a5))
extern void f_ionis(int,int,double,double *,double *,double *);
#define IONIS(a0,a1,a2,a3,a4,a5) f_ionis((a0),(a1),(a2),(a3),(a4),(a5))
extern void f_recomb(int,int,double,double *,double *);
#define RECOMB(a0,a1,a2,a3,a4) f_recomb((a0),(a1),(a2),(a3),(a4))
extern void f_recombfe(int,int,double,double *,double *);
#define RECOMBFE(a0,a1,a2,a3,a4) f_recombfe((a0),(a1),(a2),(a3),(a4))
extern void f_nrrfit(int,int,double,double *);
#define NRRFIT(a0,a1,a2,a3) f_nrrfit((a0),(a1),(a2),(a3))
extern void f_ndrfit(int,int,double,double *);
#define NDRFIT(a0,a1,a2,a3) f_ndrfit((a0),(a1),(a2),(a3))
extern void f_rrfit(int,int,double,double *);
#define RRFIT(a0,a1,a2,a3) f_rrfit((a0),(a1),(a2),(a3))
extern void f_drfit(int,int,double,double *);
#define DRFIT(a0,a1,a2,a3) f_drfit((a0),(a1),(a2),(a3))
extern void f_phfit2(int,int,int,double,double *);
#define PHFIT2(a0,a1,a2,a3,a4) f_phfit2((a0),(a1),(a2),(a3),(a4))
extern void f_cbeli(int,int,double,double *,double *,double *);
#define CBELI(a0,a1,a2,a3,a4,a5) f_cbeli((a0),(a1),(a2),(a3),(a4),(a5))
extern void f_rbeli(int,int,double,double *,double *);
#define RBELI(a0,a1,a2,a3,a4) f_rbeli((a0),(a1),(a2),(a3),(a4))
extern void f_cfit(int,int,double,double *);
#define CFIT(a0,a1,a2,a3) f_cfit((a0),(a1),(a2),(a3))
extern void f_colfit(int,int,int,double,double *,double *);
#define COLFIT(a0,a1,a2,a3,a4,a5) f_colfit((a0),(a1),(a2),(a3),(a4),(a5))
extern void f_ccolfit(int,int,int,double,double *,double *);
#define CCOLFIT(a0,a1,a2,a3,a4,a5) f_ccolfit((a0),(a1),(a2),(a3),(a4),(a5))
extern void f_ephfit2(int,int,int,double *);
#define EPHFIT2(a0,a1,a2,a3) f_ephfit2((a0),(a1),(a2),(a3))
extern void f_ecolfit(int,int,int,double *);
#define ECOLFIT(a0,a1,a2,a3) f_ecolfit((a0),(a1),(a2),(a3))
extern void f_ebeli(int,int,double *);
#define EBELI(a0,a1,a2) f_ebeli((a0),(a1),(a2))
extern double f_fu(double);
#define FU(a0) f_fu((a0))
extern void f_dxlegf(double,int,int,int,double,int,double *,int *,int *);
#define DXLEGF(a0,a1,a2,a3,a4,a5,a6,a7,a8) f_dxlegf((a0),(a1),(a2),(a3),(a4),(a5),(a6),(a7),(a8))
extern void f_njform(int,int,int *,int *,int *);
#define NJFORM(a0,a1,a2,a3,a4) f_njform((a0),(a1),(a2),(a3),(a4))
extern void f_njsum(int *,double *);
#define NJSUM(a0,a1) f_njsum((a0),(a1))
extern void f_cpydat(int,int *,int);
#define CPYDAT(a0,a1,a2) f_cpydat((a0),(a1),(a2))
extern void f_factt();
#define FACTT f_factt
extern void f_njdbug(int);
#define NJDBUG(a0) f_njdbug((a0))
