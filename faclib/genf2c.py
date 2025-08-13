
ifn = 'cf77.h'
ofn0 = 'f2c.f90'
ofn1 = 'f2c.h'
f0 = open(ofn0, 'w')
f1 = open(ofn1, 'w')

f0.write('module f2c\n')
f0.write('use iso_c_binding\n')

with open(ifn, 'r') as f:
    r = f.readlines()
    ip = -1
    for d in r:
        d = d.strip()
        if len(d) < 3:
            continue
        if d[:3] == 'PRO':
            if ip < 0:
                f0.write('\ncontains\n\n')
            ip = 1
            x = d.split('(')
            a = x[1].split(',')
            na = 0
            atp = []
            if x[0][:15] == 'PROTOCCALLSFFUN':
                sn = a[1].strip()
                s = 'extern %s f_%s('%(a[0].lower(), sn.lower())
                i0 = 3
                ifun = 1
                rtp = a[0].strip()
            else:
                sn = a[0].strip()
                s = 'extern void f_%s('%(sn.lower())
                i0 = 2
                ifun = 0
        if ip > 0 and ip < 100:
            if d[:3] == '#de':
                if s[-1] == ',':
                    s = s[:-1]
                f1.write(s+');\n')
                if na == 0:
                    f1.write('#define %s f_%s\n'%(sn, sn.lower()))
                else:
                    s = '#define %s('%sn
                    for i in range(na):
                        s = s + 'a%d,'%i
                    if s[-1] == ',':
                        s = s[:-1]
                    s = s+') f_%s('%(sn.lower())
                    for i in range(na):
                        s = s + '(a%d),'%i
                    if s[-1] == ',':
                        s = s[:-1]
                    f1.write(s+')\n')

                if ifun == 0:
                    s = 'subroutine f_%s('%sn.lower()
                else:
                    s = 'function f_%s('%sn.lower()
                nfs = len(s)
                f0.write(s)
                s = ''
                for i in range(na):
                    s = s + 'a%d,'%i
                    if i < na-1 and (i+1)%5 == 0:
                        s = s + '&\n' + (' '*nfs)
                            
                if len(s) > 0 and s[-1] == ',':
                    s = s[:-1]
                f0.write(s + ") bind(C, name='f_%s')\n"%sn.lower())
                for i in range(na):
                    b = atp[i]
                    if b == 'INT':
                        f0.write('integer(c_int), value :: a%d\n'%i)
                    elif b == 'PINT':
                        f0.write('integer(c_int) :: a%d\n'%i)
                    elif b == 'INTV':
                        f0.write('integer(c_int) :: a%d(*)\n'%i)
                    if b == 'DOUBLE':
                        f0.write('real(c_double), value :: a%d\n'%i)
                    elif b == 'PDOUBLE':
                        f0.write('real(c_double) :: a%d\n'%i)
                    elif b == 'DOUBLEV':
                        f0.write('real(c_double) :: a%d(*)\n'%i)
                    elif b == 'STRING':
                        f0.write('character(c_char) :: a%d(*)\n'%i)
                    elif b == 'ROUTINE':
                        f0.write('type(c_funptr), value :: a%d\n'%i)
                ifa = []
                if sn in ['LMQN','LMQNBC']:
                    f0.write('interface\n')
                    f0.write('subroutine c7(n,x,f,g) bind(C)\n')
                    f0.write('import :: c_int, c_double\n')
                    f0.write('integer(c_int) :: n\n')
                    f0.write('real(c_double), dimension(n) :: x, g\n')
                    f0.write('real(c_double) :: f\n')
                    f0.write('end subroutine\n')
                    f0.write('end interface\n')
                    f0.write('procedure(c7), pointer :: p7\n')
                    f0.write('call c_f_procpointer(a7,p7)\n')
                    ifa = [7]
                elif sn in ['SUBPLX']:
                    f0.write('interface\n')
                    f0.write('function c0(n,x) bind(C)\n')
                    f0.write('import :: c_int, c_double\n')
                    f0.write('integer(c_int) :: n\n')
                    f0.write('real(c_double) :: x\n')
                    f0.write('real(c_double) :: c0\n')
                    f0.write('end function\n')
                    f0.write('end interface\n')
                    f0.write('procedure(c0), pointer :: p0\n')
                    f0.write('call c_f_procpointer(a0,p0)\n')
                    ifa = [0]
                elif sn in ['LMDER']:
                    f0.write('interface\n')
                    f0.write('subroutine c0(m,n,x,fv,fj,ldfj,iflag)\n')
                    f0.write('import :: c_int, c_double\n')
                    f0.write('integer(c_int) :: m, n, ldfj, iflag\n')
                    f0.write('real(c_double), dimension(n) :: x\n')
                    f0.write('real(c_double), dimension(m) :: fv\n')
                    f0.write('real(c_double), dimension(ldfj,n) :: fj\n')
                    f0.write('end subroutine\n')
                    f0.write('end interface\n')
                    f0.write('procedure(c0), pointer :: p0\n')
                    f0.write('call c_f_procpointer(a0,p0)\n')
                    ifa = [0]
                elif sn in ['LSODE']:
                    f0.write('interface\n')
                    f0.write('subroutine c0(n, t, y, yd) bind(C)\n')
                    f0.write('import :: c_int, c_double\n')
                    f0.write('integer(c_int) :: n(*)\n')
                    f0.write('real(c_double) :: t(*),y(*),yd(*)\n')
                    f0.write('end subroutine\n')
                    f0.write('end interface\n\n')
                    f0.write('interface\n')
                    f0.write('subroutine c15(n,t,y,ml,mu,pd,nrpd) bind(C)\n')
                    f0.write('import :: c_int, c_double\n')
                    f0.write('integer(c_int) :: n(*), ml, mu, nrpd\n')
                    f0.write('real(c_double) :: t(*),y(*),pd(*)\n')
                    f0.write('end subroutine\n')
                    f0.write('end interface\n\n')
                    f0.write('procedure(c0), pointer :: p0\n')
                    f0.write('procedure(c15), pointer :: p15\n')
                    f0.write('call c_f_procpointer(a0, p0)\n')
                    f0.write('if (a16 .eq. 21 .or. a16 .eq. 24) then\n')
                    f0.write('  call c_f_procpointer(a15, p15)\n')
                    f0.write('endif\n')
                    ifa = [0, 15]
                elif sn in ['DQAGI','DQNG','DQAGS']:
                    f0.write('interface\n')
                    f0.write('function c0(x) bind(C)\n')
                    f0.write('import :: c_double\n')
                    f0.write('real(c_double) :: x\n')
                    f0.write('real(c_double) :: c0\n')
                    f0.write('end function\n')
                    f0.write('end interface\n\n')
                    f0.write('procedure(c0), pointer :: p0\n')
                    f0.write('call c_f_procpointer(a0, p0)\n')
                    ifa = [0]
                if ifun == 0:
                    s = 'call %s('%sn
                    nfs = len(s)
                    for i in range(na):
                        if i in ifa:
                            s = s + 'p%d,'%i
                        else:
                            s = s + 'a%d,'%i
                        if i < na-1 and (i+1)%5 == 0:
                            s = s + '&\n' + (' '*nfs)
                    if s[-1] == ',':
                        s = s[:-1]
                    f0.write(s+')\n')
                    f0.write('end subroutine\n\n')
                else:
                    if rtp == 'INT':
                        f0.write('integer(c_int) :: f_%s, %s\n'%(sn.lower(),sn))
                    elif rtp == 'DOUBLE':
                        f0.write('real(c_double) :: f_%s, %s\n'%(sn.lower(),sn))
                    s = 'f_%s = %s('%(sn.lower(),sn)
                    nfs = len(s)
                    for i in range(na):
                        s = s + 'a%d,'%i
                        if i < na-1 and (i+1)%5 == 0:
                            s = s + '&\n' + (' '*nfs)
                    if s[-1] == ',':
                        s = s[:-1]
                    f0.write(s+')\n')
                    f0.write('end function\n\n')
                ip = 0
                continue
            if ip > 1:
                i0 = 0
                a = d.split(',')
            for i in range(i0, len(a)):
                b = a[i].strip()                
                if len(b) == 0:
                    break
                if b[0] == '\\':
                    continue
                if b[-1] == ')':
                    b = b[:-1]
                if b == 'INT':
                    s = s + 'int,'
                elif b == 'DOUBLE':
                    s = s + 'double,'
                elif b == 'INTV' or b == 'PINT':
                    s = s + 'int *,'
                elif b == 'DOUBLEV' or b == 'PDOUBLE':
                    s = s + 'double *,'
                elif b == 'STRING':
                    s = s + 'char *,'
                elif b == 'ROUTINE':
                    if sn in ['LMQN', 'LMQNBC']:
                        fs = 'int *,double *,double *,double *'
                        fr = 'void'
                    elif sn in ['SUBPLX']:
                        fs = 'int *, double *'
                        fr = 'double'
                    elif sn in ['LMDER']:
                        fs = 'int *,int *,double *,double *,double *,int *,int *'
                        fr = 'void'
                    elif sn in ['LSODE']:
                        if i == i0:
                            fs = 'int *,double *,double *,double *'
                            fr = 'void'
                        else:
                            fs = 'int *,double *,double *,int *,int *,double *,int *'
                            fr = 'void'
                    elif sn in ['DQAGI','DQNG','DQAGS']:
                        fs = 'double *'
                        fr = 'double'
                    else:
                        fs = ''
                        fr = 'void'
                    s = s + '%s (*)(%s),'%(fr,fs)
                na = na+1
                atp.append(b)
            ip = ip+1
            continue
        if d[:3] == 'CCA':
            ip = 0
            continue
        if d[:3] == 'FCA':
            ip = -1
            continue
        if ip <= 0:
            continue

f0.write('end module f2c\n')
f0.close()
f1.close()
