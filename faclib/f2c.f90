module f2c
use iso_c_binding

contains

function f_erf(a0) bind(C, name='f_erf')
real(c_double), value :: a0
real(c_double) :: f_erf, ERF
f_erf = ERF(a0)
end function

function f_fm5p(a0) bind(C, name='f_fm5p')
real(c_double), value :: a0
real(c_double) :: f_fm5p, FM5P
f_fm5p = FM5P(a0)
end function

function f_fm5pi(a0) bind(C, name='f_fm5pi')
real(c_double), value :: a0
real(c_double) :: f_fm5pi, FM5PI
f_fm5pi = FM5PI(a0)
end function

function f_fm3p(a0) bind(C, name='f_fm3p')
real(c_double), value :: a0
real(c_double) :: f_fm3p, FM3P
f_fm3p = FM3P(a0)
end function

function f_fm3pi(a0) bind(C, name='f_fm3pi')
real(c_double), value :: a0
real(c_double) :: f_fm3pi, FM3PI
f_fm3pi = FM3PI(a0)
end function

function f_fm1p(a0) bind(C, name='f_fm1p')
real(c_double), value :: a0
real(c_double) :: f_fm1p, FM1P
f_fm1p = FM1P(a0)
end function

function f_fm1pi(a0) bind(C, name='f_fm1pi')
real(c_double), value :: a0
real(c_double) :: f_fm1pi, FM1PI
f_fm1pi = FM1PI(a0)
end function

function f_fm1m(a0) bind(C, name='f_fm1m')
real(c_double), value :: a0
real(c_double) :: f_fm1m, FM1M
f_fm1m = FM1M(a0)
end function

function f_fm1mi(a0) bind(C, name='f_fm1mi')
real(c_double), value :: a0
real(c_double) :: f_fm1mi, FM1MI
f_fm1mi = FM1MI(a0)
end function

function f_expint(a0,a1) bind(C, name='f_expint')
real(c_double), value :: a0
integer(c_int), value :: a1
real(c_double) :: f_expint, EXPINT
f_expint = EXPINT(a0,a1)
end function

function f_eione(a0,a1) bind(C, name='f_eione')
integer(c_int), value :: a0
real(c_double), value :: a1
real(c_double) :: f_eione, EIONE
f_eione = EIONE(a0,a1)
end function

function f_d1mach(a0) bind(C, name='f_d1mach')
integer(c_int), value :: a0
real(c_double) :: f_d1mach, D1MACH
f_d1mach = D1MACH(a0)
end function

subroutine f_radfnd(a0,a1) bind(C, name='f_radfnd')
integer(c_int), value :: a0
real(c_double) :: a1(*)
call RADFND(a0,a1)
end subroutine

subroutine f_locsep(a0,a1,a2,a3) bind(C, name='f_locsep')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double) :: a2(*)
real(c_double) :: a3(*)
call LOCSEP(a0,a1,a2,a3)
end subroutine

subroutine f_rgmqed(a0,a1) bind(C, name='f_rgmqed')
real(c_double) :: a0(*)
real(c_double) :: a1(*)
call RGMQED(a0,a1)
end subroutine

subroutine f_fsedat(a0,a1,a2,a3,a4,&
                    a5) bind(C, name='f_fsedat')
integer(c_int), value :: a0
integer(c_int), value :: a1
integer(c_int), value :: a2
integer(c_int), value :: a3
real(c_double) :: a4(*)
real(c_double) :: a5(*)
call FSEDAT(a0,a1,a2,a3,a4,&
            a5)
end subroutine

subroutine f_genqed(a0,a1,a2,a3,a4,&
                    a5,a6,a7,a8,a9,&
                    a10) bind(C, name='f_genqed')
integer(c_int), value :: a0
integer(c_int), value :: a1
integer(c_int), value :: a2
integer(c_int), value :: a3
integer(c_int), value :: a4
real(c_double) :: a5(*)
real(c_double) :: a6(*)
real(c_double) :: a7(*)
real(c_double) :: a8(*)
real(c_double) :: a9(*)
real(c_double) :: a10(*)
call GENQED(a0,a1,a2,a3,a4,&
            a5,a6,a7,a8,a9,&
            a10)
end subroutine

subroutine f_iniqed(a0,a1,a2,a3) bind(C, name='f_iniqed')
real(c_double), value :: a0
integer(c_int), value :: a1
integer(c_int), value :: a2
real(c_double), value :: a3
call INIQED(a0,a1,a2,a3)
end subroutine

subroutine f_mohrfin(a0,a1,a2,a3,a4,&
                     a5,a6,a7,a8) bind(C, name='f_mohrfin')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double), value :: a2
real(c_double), value :: a3
real(c_double) :: a4(*)
real(c_double) :: a5(*)
real(c_double) :: a6(*)
real(c_double) :: a7(*)
real(c_double) :: a8(*)
call MOHRFIN(a0,a1,a2,a3,a4,&
             a5,a6,a7,a8)
end subroutine

subroutine f_uvip3p(a0,a1,a2,a3,a4,&
                    a5,a6) bind(C, name='f_uvip3p')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double) :: a2(*)
real(c_double) :: a3(*)
integer(c_int), value :: a4
real(c_double) :: a5(*)
real(c_double) :: a6(*)
call UVIP3P(a0,a1,a2,a3,a4,&
            a5,a6)
end subroutine

subroutine f_uvip3i(a0,a1,a2,a3,a4,&
                    a5,a6,a7) bind(C, name='f_uvip3i')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double) :: a2(*)
real(c_double) :: a3(*)
integer(c_int), value :: a4
real(c_double) :: a5(*)
real(c_double) :: a6(*)
real(c_double) :: a7(*)
call UVIP3I(a0,a1,a2,a3,a4,&
            a5,a6,a7)
end subroutine

subroutine f_uvip3c(a0,a1,a2,a3,a4,&
                    a5,a6) bind(C, name='f_uvip3c')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double) :: a2(*)
real(c_double) :: a3(*)
real(c_double) :: a4(*)
real(c_double) :: a5(*)
real(c_double) :: a6(*)
call UVIP3C(a0,a1,a2,a3,a4,&
            a5,a6)
end subroutine

subroutine f_sdbi3p(a0,a1,a2,a3,a4,&
                    a5,a6,a7,a8,a9,&
                    a10,a11) bind(C, name='f_sdbi3p')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double) :: a2(*)
real(c_double) :: a3(*)
real(c_double) :: a4(*)
integer(c_int), value :: a5
real(c_double) :: a6(*)
real(c_double) :: a7(*)
real(c_double) :: a8(*)
integer(c_int) :: a9(*)
real(c_double) :: a10(*)
integer(c_int) :: a11(*)
call SDBI3P(a0,a1,a2,a3,a4,&
            a5,a6,a7,a8,a9,&
            a10,a11)
end subroutine

subroutine f_rgbi3p(a0,a1,a2,a3,a4,&
                    a5,a6,a7,a8,a9,&
                    a10,a11) bind(C, name='f_rgbi3p')
integer(c_int), value :: a0
integer(c_int), value :: a1
integer(c_int), value :: a2
real(c_double) :: a3(*)
real(c_double) :: a4(*)
real(c_double) :: a5(*)
integer(c_int), value :: a6
real(c_double) :: a7(*)
real(c_double) :: a8(*)
real(c_double) :: a9(*)
integer(c_int) :: a10(*)
real(c_double) :: a11(*)
call RGBI3P(a0,a1,a2,a3,a4,&
            a5,a6,a7,a8,a9,&
            a10,a11)
end subroutine

subroutine f_lmqn(a0,a1,a2,a3,a4,&
                  a5,a6,a7,a8,a9,&
                  a10,a11,a12,a13,a14) bind(C, name='f_lmqn')
integer(c_int) :: a0(*)
integer(c_int), value :: a1
real(c_double) :: a2(*)
real(c_double) :: a3(*)
real(c_double) :: a4(*)
real(c_double) :: a5(*)
integer(c_int), value :: a6
type(c_funptr), value :: a7
integer(c_int), value :: a8
integer(c_int), value :: a9
integer(c_int), value :: a10
real(c_double), value :: a11
real(c_double), value :: a12
real(c_double), value :: a13
real(c_double), value :: a14
interface
subroutine c7(n,x,f,g) bind(C)
import :: c_int, c_double
integer(c_int) :: n
real(c_double), dimension(n) :: x, g
real(c_double) :: f
end subroutine
end interface
procedure(c7), pointer :: p7
call c_f_procpointer(a7,p7)
call LMQN(a0,a1,a2,a3,a4,&
          a5,a6,p7,a8,a9,&
          a10,a11,a12,a13,a14)
end subroutine

subroutine f_lmqnbc(a0,a1,a2,a3,a4,&
                    a5,a6,a7,a8,a9,&
                    a10,a11,a12,a13,a14,&
                    a15,a16,a17) bind(C, name='f_lmqnbc')
integer(c_int) :: a0(*)
integer(c_int), value :: a1
real(c_double) :: a2(*)
real(c_double) :: a3(*)
real(c_double) :: a4(*)
real(c_double) :: a5(*)
integer(c_int), value :: a6
type(c_funptr), value :: a7
real(c_double) :: a8(*)
real(c_double) :: a9(*)
integer(c_int) :: a10(*)
integer(c_int), value :: a11
integer(c_int), value :: a12
integer(c_int), value :: a13
real(c_double), value :: a14
real(c_double), value :: a15
real(c_double), value :: a16
real(c_double), value :: a17
interface
subroutine c7(n,x,f,g) bind(C)
import :: c_int, c_double
integer(c_int) :: n
real(c_double), dimension(n) :: x, g
real(c_double) :: f
end subroutine
end interface
procedure(c7), pointer :: p7
call c_f_procpointer(a7,p7)
call LMQNBC(a0,a1,a2,a3,a4,&
            a5,a6,p7,a8,a9,&
            a10,a11,a12,a13,a14,&
            a15,a16,a17)
end subroutine

subroutine f_subplx(a0,a1,a2,a3,a4,&
                    a5,a6,a7,a8,a9,&
                    a10,a11) bind(C, name='f_subplx')
type(c_funptr), value :: a0
integer(c_int), value :: a1
real(c_double), value :: a2
integer(c_int), value :: a3
integer(c_int), value :: a4
real(c_double) :: a5(*)
real(c_double) :: a6(*)
real(c_double) :: a7(*)
integer(c_int) :: a8(*)
real(c_double) :: a9(*)
integer(c_int) :: a10(*)
integer(c_int) :: a11(*)
interface
function c0(n,x) bind(C)
import :: c_int, c_double
integer(c_int) :: n
real(c_double) :: x
real(c_double) :: c0
end function
end interface
procedure(c0), pointer :: p0
call c_f_procpointer(a0,p0)
call SUBPLX(p0,a1,a2,a3,a4,&
            a5,a6,a7,a8,a9,&
            a10,a11)
end subroutine

function f_argam(a0,a1) bind(C, name='f_argam')
real(c_double), value :: a0
real(c_double), value :: a1
real(c_double) :: f_argam, ARGAM
f_argam = ARGAM(a0,a1)
end function

function f_dlogam(a0) bind(C, name='f_dlogam')
real(c_double), value :: a0
real(c_double) :: f_dlogam, DLOGAM
f_dlogam = DLOGAM(a0)
end function

subroutine f_beslik(a0,a1,a2,a3) bind(C, name='f_beslik')
real(c_double), value :: a0
real(c_double), value :: a1
real(c_double) :: a2(*)
real(c_double) :: a3(*)
call BESLIK(a0,a1,a2,a3)
end subroutine

function f_besljn(a0,a1,a2) bind(C, name='f_besljn')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double), value :: a2
real(c_double) :: f_besljn, BESLJN
f_besljn = BESLJN(a0,a1,a2)
end function

subroutine f_y5n(a0,a1,a2,a3,a4,&
                 a5,a6,a7,a8) bind(C, name='f_y5n')
real(c_double), value :: a0
real(c_double), value :: a1
real(c_double), value :: a2
real(c_double), value :: a3
real(c_double) :: a4(*)
real(c_double) :: a5(*)
real(c_double) :: a6(*)
real(c_double) :: a7(*)
integer(c_int) :: a8(*)
call Y5N(a0,a1,a2,a3,a4,&
         a5,a6,a7,a8)
end subroutine

subroutine f_dcoul(a0,a1,a2,a3,a4,&
                   a5,a6,a7,a8) bind(C, name='f_dcoul')
real(c_double), value :: a0
real(c_double), value :: a1
integer(c_int), value :: a2
real(c_double), value :: a3
real(c_double) :: a4(*)
real(c_double) :: a5(*)
real(c_double) :: a6(*)
real(c_double) :: a7(*)
integer(c_int) :: a8(*)
call DCOUL(a0,a1,a2,a3,a4,&
           a5,a6,a7,a8)
end subroutine

subroutine f_dcoul1(a0,a1,a2,a3,a4,&
                    a5,a6,a7,a8) bind(C, name='f_dcoul1')
real(c_double), value :: a0
real(c_double), value :: a1
integer(c_int), value :: a2
real(c_double), value :: a3
real(c_double) :: a4(*)
real(c_double) :: a5(*)
real(c_double) :: a6(*)
real(c_double) :: a7(*)
integer(c_int) :: a8(*)
call DCOUL1(a0,a1,a2,a3,a4,&
            a5,a6,a7,a8)
end subroutine

subroutine f_cmultip(a0,a1,a2,a3,a4,&
                     a5,a6,a7,a8,a9,&
                     a10) bind(C, name='f_cmultip')
real(c_double), value :: a0
real(c_double), value :: a1
real(c_double), value :: a2
real(c_double), value :: a3
real(c_double), value :: a4
integer(c_int), value :: a5
integer(c_int), value :: a6
integer(c_int), value :: a7
real(c_double) :: a8(*)
integer(c_int), value :: a9
integer(c_int) :: a10(*)
call CMULTIP(a0,a1,a2,a3,a4,&
             a5,a6,a7,a8,a9,&
             a10)
end subroutine

subroutine f_acofz1(a0,a1,a2,a3,a4,&
                    a5,a6,a7) bind(C, name='f_acofz1')
real(c_double), value :: a0
real(c_double), value :: a1
integer(c_int), value :: a2
integer(c_int), value :: a3
real(c_double) :: a4(*)
real(c_double) :: a5(*)
integer(c_int), value :: a6
integer(c_int), value :: a7
call ACOFZ1(a0,a1,a2,a3,a4,&
            a5,a6,a7)
end subroutine

subroutine f_pixz1(a0,a1,a2,a3,a4,&
                   a5,a6,a7,a8,a9,&
                   a10,a11) bind(C, name='f_pixz1')
real(c_double), value :: a0
real(c_double), value :: a1
integer(c_int), value :: a2
integer(c_int), value :: a3
real(c_double) :: a4(*)
real(c_double) :: a5(*)
real(c_double) :: a6(*)
real(c_double) :: a7(*)
integer(c_int), value :: a8
integer(c_int), value :: a9
integer(c_int), value :: a10
integer(c_int), value :: a11
call PIXZ1(a0,a1,a2,a3,a4,&
           a5,a6,a7,a8,a9,&
           a10,a11)
end subroutine

subroutine f_lmder(a0,a1,a2,a3,a4,&
                   a5,a6,a7,a8,a9,&
                   a10,a11,a12,a13,a14,&
                   a15,a16,a17,a18,a19,&
                   a20,a21,a22,a23) bind(C, name='f_lmder')
type(c_funptr), value :: a0
integer(c_int), value :: a1
integer(c_int), value :: a2
real(c_double) :: a3(*)
real(c_double) :: a4(*)
real(c_double) :: a5(*)
integer(c_int), value :: a6
real(c_double), value :: a7
real(c_double), value :: a8
real(c_double), value :: a9
integer(c_int), value :: a10
real(c_double) :: a11(*)
integer(c_int), value :: a12
real(c_double), value :: a13
integer(c_int), value :: a14
integer(c_int) :: a15(*)
integer(c_int) :: a16(*)
integer(c_int) :: a17(*)
integer(c_int) :: a18(*)
real(c_double) :: a19(*)
real(c_double) :: a20(*)
real(c_double) :: a21(*)
real(c_double) :: a22(*)
real(c_double) :: a23(*)
interface
subroutine c0(m,n,x,fv,fj,ldfj,iflag)
import :: c_int, c_double
integer(c_int) :: m, n, ldfj, iflag
real(c_double), dimension(n) :: x
real(c_double), dimension(m) :: fv
real(c_double), dimension(ldfj,n) :: fj
end subroutine
end interface
procedure(c0), pointer :: p0
call c_f_procpointer(a0,p0)
call LMDER(p0,a1,a2,a3,a4,&
           a5,a6,a7,a8,a9,&
           a10,a11,a12,a13,a14,&
           a15,a16,a17,a18,a19,&
           a20,a21,a22,a23)
end subroutine

subroutine f_chkder(a0,a1,a2,a3,a4,&
                    a5,a6,a7,a8,a9) bind(C, name='f_chkder')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double) :: a2(*)
real(c_double) :: a3(*)
real(c_double) :: a4(*)
integer(c_int), value :: a5
real(c_double) :: a6(*)
real(c_double) :: a7(*)
integer(c_int), value :: a8
real(c_double) :: a9(*)
call CHKDER(a0,a1,a2,a3,a4,&
            a5,a6,a7,a8,a9)
end subroutine

function f_ddot(a0,a1,a2,a3,a4) bind(C, name='f_ddot')
integer(c_int), value :: a0
real(c_double) :: a1(*)
integer(c_int), value :: a2
real(c_double) :: a3(*)
integer(c_int), value :: a4
real(c_double) :: f_ddot, DDOT
f_ddot = DDOT(a0,a1,a2,a3,a4)
end function

subroutine f_dscal(a0,a1,a2,a3) bind(C, name='f_dscal')
integer(c_int), value :: a0
real(c_double), value :: a1
real(c_double) :: a2(*)
integer(c_int), value :: a3
call DSCAL(a0,a1,a2,a3)
end subroutine

subroutine f_dgemv(a0,a1,a2,a3,a4,&
                   a5,a6,a7,a8,a9,&
                   a10) bind(C, name='f_dgemv')
character(c_char) :: a0(*)
integer(c_int), value :: a1
integer(c_int), value :: a2
real(c_double), value :: a3
real(c_double) :: a4(*)
integer(c_int), value :: a5
real(c_double) :: a6(*)
integer(c_int), value :: a7
real(c_double), value :: a8
real(c_double) :: a9(*)
integer(c_int), value :: a10
call DGEMV(a0,a1,a2,a3,a4,&
           a5,a6,a7,a8,a9,&
           a10)
end subroutine

subroutine f_dsbev(a0,a1,a2,a3,a4,&
                   a5,a6,a7,a8,a9,&
                   a10) bind(C, name='f_dsbev')
character(c_char) :: a0(*)
character(c_char) :: a1(*)
integer(c_int), value :: a2
integer(c_int), value :: a3
real(c_double) :: a4(*)
integer(c_int), value :: a5
real(c_double) :: a6(*)
real(c_double) :: a7(*)
integer(c_int), value :: a8
real(c_double) :: a9(*)
integer(c_int) :: a10(*)
call DSBEV(a0,a1,a2,a3,a4,&
           a5,a6,a7,a8,a9,&
           a10)
end subroutine

subroutine f_dspevd(a0,a1,a2,a3,a4,&
                    a5,a6,a7,a8,a9,&
                    a10,a11) bind(C, name='f_dspevd')
character(c_char) :: a0(*)
character(c_char) :: a1(*)
integer(c_int), value :: a2
real(c_double) :: a3(*)
real(c_double) :: a4(*)
real(c_double) :: a5(*)
integer(c_int), value :: a6
real(c_double) :: a7(*)
integer(c_int), value :: a8
integer(c_int) :: a9(*)
integer(c_int), value :: a10
integer(c_int) :: a11(*)
call DSPEVD(a0,a1,a2,a3,a4,&
            a5,a6,a7,a8,a9,&
            a10,a11)
end subroutine

subroutine f_dspev(a0,a1,a2,a3,a4,&
                   a5,a6,a7,a8) bind(C, name='f_dspev')
character(c_char) :: a0(*)
character(c_char) :: a1(*)
integer(c_int), value :: a2
real(c_double) :: a3(*)
real(c_double) :: a4(*)
real(c_double) :: a5(*)
integer(c_int), value :: a6
real(c_double) :: a7(*)
integer(c_int) :: a8(*)
call DSPEV(a0,a1,a2,a3,a4,&
           a5,a6,a7,a8)
end subroutine

subroutine f_dgeev(a0,a1,a2,a3,a4,&
                   a5,a6,a7,a8,a9,&
                   a10,a11,a12,a13) bind(C, name='f_dgeev')
character(c_char) :: a0(*)
character(c_char) :: a1(*)
integer(c_int), value :: a2
real(c_double) :: a3(*)
integer(c_int), value :: a4
real(c_double) :: a5(*)
real(c_double) :: a6(*)
real(c_double) :: a7(*)
integer(c_int), value :: a8
real(c_double) :: a9(*)
integer(c_int), value :: a10
real(c_double) :: a11(*)
integer(c_int), value :: a12
integer(c_int) :: a13(*)
call DGEEV(a0,a1,a2,a3,a4,&
           a5,a6,a7,a8,a9,&
           a10,a11,a12,a13)
end subroutine

subroutine f_dgesv(a0,a1,a2,a3,a4,&
                   a5,a6,a7) bind(C, name='f_dgesv')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double) :: a2(*)
integer(c_int), value :: a3
integer(c_int) :: a4(*)
real(c_double) :: a5(*)
integer(c_int), value :: a6
integer(c_int) :: a7(*)
call DGESV(a0,a1,a2,a3,a4,&
           a5,a6,a7)
end subroutine

subroutine f_dgesdd(a0,a1,a2,a3,a4,&
                    a5,a6,a7,a8,a9,&
                    a10,a11,a12,a13) bind(C, name='f_dgesdd')
character(c_char) :: a0(*)
integer(c_int), value :: a1
integer(c_int), value :: a2
real(c_double) :: a3(*)
integer(c_int), value :: a4
real(c_double) :: a5(*)
real(c_double) :: a6(*)
integer(c_int), value :: a7
real(c_double) :: a8(*)
integer(c_int), value :: a9
real(c_double) :: a10(*)
integer(c_int), value :: a11
integer(c_int) :: a12(*)
integer(c_int) :: a13(*)
call DGESDD(a0,a1,a2,a3,a4,&
            a5,a6,a7,a8,a9,&
            a10,a11,a12,a13)
end subroutine

subroutine f_dgbsv(a0,a1,a2,a3,a4,&
                   a5,a6,a7,a8,a9) bind(C, name='f_dgbsv')
integer(c_int), value :: a0
integer(c_int), value :: a1
integer(c_int), value :: a2
integer(c_int), value :: a3
real(c_double) :: a4(*)
integer(c_int), value :: a5
integer(c_int) :: a6(*)
real(c_double) :: a7(*)
integer(c_int), value :: a8
integer(c_int) :: a9(*)
call DGBSV(a0,a1,a2,a3,a4,&
           a5,a6,a7,a8,a9)
end subroutine

subroutine f_lsode(a0,a1,a2,a3,a4,&
                   a5,a6,a7,a8,a9,&
                   a10,a11,a12,a13,a14,&
                   a15,a16) bind(C, name='f_lsode')
type(c_funptr), value :: a0
integer(c_int) :: a1(*)
real(c_double) :: a2(*)
real(c_double) :: a3(*)
real(c_double), value :: a4
integer(c_int), value :: a5
real(c_double), value :: a6
real(c_double) :: a7(*)
integer(c_int), value :: a8
integer(c_int) :: a9(*)
integer(c_int), value :: a10
real(c_double) :: a11(*)
integer(c_int), value :: a12
integer(c_int) :: a13(*)
integer(c_int), value :: a14
type(c_funptr), value :: a15
integer(c_int), value :: a16
interface
subroutine c0(n, t, y, yd) bind(C)
import :: c_int, c_double
integer(c_int) :: n(*)
real(c_double) :: t(*),y(*),yd(*)
end subroutine
end interface

interface
subroutine c15(n,t,y,ml,mu,pd,nrpd) bind(C)
import :: c_int, c_double
integer(c_int) :: n(*), ml, mu, nrpd
real(c_double) :: t(*),y(*),pd(*)
end subroutine
end interface

procedure(c0), pointer :: p0
procedure(c15), pointer :: p15
call c_f_procpointer(a0, p0)
if (a16 .eq. 21 .or. a16 .eq. 24) then
  call c_f_procpointer(a15, p15)
endif
call LSODE(p0,a1,a2,a3,a4,&
           a5,a6,a7,a8,a9,&
           a10,a11,a12,a13,a14,&
           p15,a16)
end subroutine

subroutine f_dqagi(a0,a1,a2,a3,a4,&
                   a5,a6,a7,a8,a9,&
                   a10,a11,a12,a13) bind(C, name='f_dqagi')
type(c_funptr), value :: a0
real(c_double), value :: a1
integer(c_int), value :: a2
real(c_double), value :: a3
real(c_double), value :: a4
real(c_double) :: a5(*)
real(c_double) :: a6(*)
integer(c_int) :: a7(*)
integer(c_int) :: a8(*)
integer(c_int), value :: a9
integer(c_int), value :: a10
integer(c_int) :: a11(*)
integer(c_int) :: a12(*)
real(c_double) :: a13(*)
interface
function c0(x) bind(C)
import :: c_double
real(c_double) :: x
real(c_double) :: c0
end function
end interface

procedure(c0), pointer :: p0
call c_f_procpointer(a0, p0)
call DQAGI(p0,a1,a2,a3,a4,&
           a5,a6,a7,a8,a9,&
           a10,a11,a12,a13)
end subroutine

subroutine f_dqng(a0,a1,a2,a3,a4,&
                  a5,a6,a7,a8) bind(C, name='f_dqng')
type(c_funptr), value :: a0
real(c_double), value :: a1
real(c_double), value :: a2
real(c_double), value :: a3
real(c_double), value :: a4
real(c_double) :: a5(*)
real(c_double) :: a6(*)
integer(c_int) :: a7(*)
integer(c_int) :: a8(*)
interface
function c0(x) bind(C)
import :: c_double
real(c_double) :: x
real(c_double) :: c0
end function
end interface

procedure(c0), pointer :: p0
call c_f_procpointer(a0, p0)
call DQNG(p0,a1,a2,a3,a4,&
          a5,a6,a7,a8)
end subroutine

subroutine f_dqags(a0,a1,a2,a3,a4,&
                   a5,a6,a7,a8,a9,&
                   a10,a11,a12,a13) bind(C, name='f_dqags')
type(c_funptr), value :: a0
real(c_double), value :: a1
real(c_double), value :: a2
real(c_double), value :: a3
real(c_double), value :: a4
real(c_double) :: a5(*)
real(c_double) :: a6(*)
integer(c_int) :: a7(*)
integer(c_int) :: a8(*)
integer(c_int), value :: a9
integer(c_int), value :: a10
integer(c_int) :: a11(*)
integer(c_int) :: a12(*)
real(c_double) :: a13(*)
interface
function c0(x) bind(C)
import :: c_double
real(c_double) :: x
real(c_double) :: c0
end function
end interface

procedure(c0), pointer :: p0
call c_f_procpointer(a0, p0)
call DQAGS(p0,a1,a2,a3,a4,&
           a5,a6,a7,a8,a9,&
           a10,a11,a12,a13)
end subroutine

subroutine f_fermid(a0,a1,a2,a3,a4) bind(C, name='f_fermid')
real(c_double), value :: a0
real(c_double), value :: a1
real(c_double), value :: a2
real(c_double) :: a3(*)
integer(c_int) :: a4(*)
call FERMID(a0,a1,a2,a3,a4)
end subroutine

subroutine f_ferinc(a0,a1,a2,a3,a4,&
                    a5) bind(C, name='f_ferinc')
real(c_double), value :: a0
real(c_double), value :: a1
real(c_double), value :: a2
real(c_double), value :: a3
real(c_double) :: a4(*)
integer(c_int) :: a5(*)
call FERINC(a0,a1,a2,a3,a4,&
            a5)
end subroutine

subroutine f_ionis(a0,a1,a2,a3,a4,&
                   a5) bind(C, name='f_ionis')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double), value :: a2
real(c_double) :: a3(*)
real(c_double) :: a4(*)
real(c_double) :: a5(*)
call IONIS(a0,a1,a2,a3,a4,&
           a5)
end subroutine

subroutine f_recomb(a0,a1,a2,a3,a4) bind(C, name='f_recomb')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double), value :: a2
real(c_double) :: a3(*)
real(c_double) :: a4(*)
call RECOMB(a0,a1,a2,a3,a4)
end subroutine

subroutine f_recombfe(a0,a1,a2,a3,a4) bind(C, name='f_recombfe')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double), value :: a2
real(c_double) :: a3(*)
real(c_double) :: a4(*)
call RECOMBFE(a0,a1,a2,a3,a4)
end subroutine

subroutine f_nrrfit(a0,a1,a2,a3) bind(C, name='f_nrrfit')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double), value :: a2
real(c_double) :: a3(*)
call NRRFIT(a0,a1,a2,a3)
end subroutine

subroutine f_ndrfit(a0,a1,a2,a3) bind(C, name='f_ndrfit')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double), value :: a2
real(c_double) :: a3(*)
call NDRFIT(a0,a1,a2,a3)
end subroutine

subroutine f_rrfit(a0,a1,a2,a3) bind(C, name='f_rrfit')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double), value :: a2
real(c_double) :: a3(*)
call RRFIT(a0,a1,a2,a3)
end subroutine

subroutine f_drfit(a0,a1,a2,a3) bind(C, name='f_drfit')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double), value :: a2
real(c_double) :: a3(*)
call DRFIT(a0,a1,a2,a3)
end subroutine

subroutine f_phfit2(a0,a1,a2,a3,a4) bind(C, name='f_phfit2')
integer(c_int), value :: a0
integer(c_int), value :: a1
integer(c_int), value :: a2
real(c_double), value :: a3
real(c_double) :: a4(*)
call PHFIT2(a0,a1,a2,a3,a4)
end subroutine

subroutine f_cbeli(a0,a1,a2,a3,a4,&
                   a5) bind(C, name='f_cbeli')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double), value :: a2
real(c_double) :: a3(*)
real(c_double) :: a4(*)
real(c_double) :: a5(*)
call CBELI(a0,a1,a2,a3,a4,&
           a5)
end subroutine

subroutine f_rbeli(a0,a1,a2,a3,a4) bind(C, name='f_rbeli')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double), value :: a2
real(c_double) :: a3(*)
real(c_double) :: a4(*)
call RBELI(a0,a1,a2,a3,a4)
end subroutine

subroutine f_cfit(a0,a1,a2,a3) bind(C, name='f_cfit')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double), value :: a2
real(c_double) :: a3(*)
call CFIT(a0,a1,a2,a3)
end subroutine

subroutine f_colfit(a0,a1,a2,a3,a4,&
                    a5) bind(C, name='f_colfit')
integer(c_int), value :: a0
integer(c_int), value :: a1
integer(c_int), value :: a2
real(c_double), value :: a3
real(c_double) :: a4(*)
real(c_double) :: a5(*)
call COLFIT(a0,a1,a2,a3,a4,&
            a5)
end subroutine

subroutine f_ccolfit(a0,a1,a2,a3,a4,&
                     a5) bind(C, name='f_ccolfit')
integer(c_int), value :: a0
integer(c_int), value :: a1
integer(c_int), value :: a2
real(c_double), value :: a3
real(c_double) :: a4(*)
real(c_double) :: a5(*)
call CCOLFIT(a0,a1,a2,a3,a4,&
             a5)
end subroutine

subroutine f_ephfit2(a0,a1,a2,a3) bind(C, name='f_ephfit2')
integer(c_int), value :: a0
integer(c_int), value :: a1
integer(c_int), value :: a2
real(c_double) :: a3(*)
call EPHFIT2(a0,a1,a2,a3)
end subroutine

subroutine f_ecolfit(a0,a1,a2,a3) bind(C, name='f_ecolfit')
integer(c_int), value :: a0
integer(c_int), value :: a1
integer(c_int), value :: a2
real(c_double) :: a3(*)
call ECOLFIT(a0,a1,a2,a3)
end subroutine

subroutine f_ebeli(a0,a1,a2) bind(C, name='f_ebeli')
integer(c_int), value :: a0
integer(c_int), value :: a1
real(c_double) :: a2(*)
call EBELI(a0,a1,a2)
end subroutine

function f_fu(a0) bind(C, name='f_fu')
real(c_double), value :: a0
real(c_double) :: f_fu, FU
f_fu = FU(a0)
end function

subroutine f_dxlegf(a0,a1,a2,a3,a4,&
                    a5,a6,a7,a8) bind(C, name='f_dxlegf')
real(c_double), value :: a0
integer(c_int), value :: a1
integer(c_int), value :: a2
integer(c_int), value :: a3
real(c_double), value :: a4
integer(c_int), value :: a5
real(c_double) :: a6(*)
integer(c_int) :: a7(*)
integer(c_int) :: a8(*)
call DXLEGF(a0,a1,a2,a3,a4,&
            a5,a6,a7,a8)
end subroutine

subroutine f_njform(a0,a1,a2,a3,a4) bind(C, name='f_njform')
integer(c_int), value :: a0
integer(c_int), value :: a1
integer(c_int) :: a2(*)
integer(c_int) :: a3(*)
integer(c_int) :: a4(*)
call NJFORM(a0,a1,a2,a3,a4)
end subroutine

subroutine f_njsum(a0,a1) bind(C, name='f_njsum')
integer(c_int) :: a0(*)
real(c_double) :: a1(*)
call NJSUM(a0,a1)
end subroutine

subroutine f_cpydat(a0,a1,a2) bind(C, name='f_cpydat')
integer(c_int), value :: a0
integer(c_int) :: a1(*)
integer(c_int), value :: a2
call CPYDAT(a0,a1,a2)
end subroutine

subroutine f_factt() bind(C, name='f_factt')
call FACTT()
end subroutine

subroutine f_njdbug(a0) bind(C, name='f_njdbug')
integer(c_int), value :: a0
call NJDBUG(a0)
end subroutine

end module f2c
