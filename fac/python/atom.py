from fac import *
from types import *
import time

# reinitialize FAC
def reinit(m = 0):
    Reinit(m)

    return

# the index of n-th shell in a complex
def get_index(n, complex):
    for i in range(len(complex)):
        if (complex[i][0] == n):
            return i
        elif (complex[i][0] > n):
            return -i-1
        
    return -i-2

def get_terms(n, nq):
    if (nq == 1):
        t = ['%d*1'%n]
    elif (nq == 2):
        if (n == 2):
            t = ['2s2', '2s1 2p1', '2p2']
        else:
            t = ['%d*2'%n]
    else:
        if (n == 2):
            t = ['2s2 2p%d'%(nq-2)]
            if (nq < 8):
                t.append('2s1 2p%d'%(nq-1))
                if (nq < 7):
                    t.append('2p%d'%nq)
            if (nq > 8):
                raise 'nq <= 8 for n = 2'
        elif (n == 3):
            if (nq < 8):
                t = ['3s2 3p%d'%(nq-2)]
                t.append('3s1 3p%d'%(nq-1))
                t.append('3s1 3p%d 3d1'%(nq-2))
                if (nq == 3):
                    t.append('3s2 3d1')
                else:
                    t.append('3s2 3p%d 3d1'%(nq-3))
            else:
                if (nq != 8):
                    t = ['3s2 3p6 3d%d'%(nq-8)]
                else:
                    t = ['3s2 2p6']
                t.append('3s2 3p5 3d%d'%(nq-7))
                t.append('3s1 3p6 3d%d'%(nq-7))
        else:
            raise 'only n <= 3 supports nq > 2'

    return t

    
# a configuration complex
class COMPLEX:
    def __init__(self, bname, nele = 0):
        self.complex = []
        self.terms = []
        self.nrec = 0
        self.nrec_ext = 0
        self.bname = bname
        self.name = []

        if (nele > 0):
            self.set_ground(nele)
            
        return

    def set_terms(self):
        t = []
        for c in self.complex:
            t.append(get_terms(c[0], c[1]))
        if (len(t) == 0):
            self.terms = ['']
            return
        a = t[0]
        for b in t[1:]:
            c = []
            for i in a:
                for j in b:
                    c.append(i+' '+j)
            a = c
        self.terms = a
                       
        return

    def set_name(self, name = []):
        if (len(name) > 0):
            self.name = name
            return
        
        a = []
        for i in range(len(self.terms)):
            a.append('%s%d'%(self.bname, i))
        self.name = a

        return

    def set_ground(self, nelectrons):
        n = 1
        nele = nelectrons
        g = []
        while (nele > 0):
            nqm = 2*n*n
            if (nele < nqm):
                g.append((n, nele))
            else:
                g.append((n, nqm))
            nele = nele - nqm
            n = n+1

        self.complex = g
        self.set_terms()
        self.set_name()
        
        return 

    def set_excited(self, i0, i1, cbase):
        base = cbase.complex
        if (i0 >= i1):
            return -1
        i = get_index(i0, base)
        if (i < 0):
            return -1
        j = get_index(i1, base)
        ex = base[:]
        if (j >= 0):
            nq1 = ex[j][1]
            if (nq1 == 2*i1*i1):
                return 0
            else:
                ex[j] = (i1, nq1+1)
        else:
            j = -j - 1
            ex.insert(j, (i1, 1))

        nq0 = ex[i][1]
        if (nq0 == 1):
            ex.remove(ex[i])
        else:
            ex[i] = (i0, nq0-1)
    
        self.complex = ex
        self.set_terms()
        self.set_name()

        return 0

    def set_ionized(self, i0, cbase):
        base = cbase.complex
        i = get_index(i0, base)
        if (i < 0):
            return -1
        ion = base[:]
        nq0 = ion[i][1]
        if (nq0 == 1):
            ion.remove(ion[i])
        else:
            ion[i] = (i0, nq0-1)

        self.complex = ion
        self.set_terms()
        self.set_name()
        
        return 0

    def set_recombined(self, i0, bname):
        self.set_name(bname)
        self.nrec = i0
        self.complex = (bname, i0)
        
        return 0
    

# a group of complexes
class CGROUP:
    def __init__(self, type):
        self.cgroup = []
        self.type = type

        return

    def add_complex(self, complex):
        self.cgroup.append(complex)

        return


# a generic atomic ion
class ATOM:    
    def __init__(self, nele, asym=0, dir=''):
        if (nele <= 0):
            raise 'NELE must > 0'
        
	self.process = {'ce': 1, 
		        'tr': 1, 
			'rr': 1, 
			'ci': 1,
			'ai': 1}

        self.ecorrections = []
        self.nele = nele
        self.nele_max = [0, 2, 10, 28]
        self.nele_sim = [6,7,8,9]
        
        self.grd_complex = COMPLEX('grd.', nele)
        self.exc_complex = []
        self.ion_complex = CGROUP('ion')

        self.angz_cut0 = 1E-6
        self.tr_cut0 = 1E-3
        self.ai_cut0 = 1E-3
        if (self.nele <= self.nele_max[1]):
            self.n_shells = 1
            self.nterms = [-1,-1,-1]
            self.nexc_max = [8, 8, 6]
            self.nfrozen = [10, 10, 10]
            self.nexc_rec = [10, 6, 5]
            self.nrec_max = [10, 10, 10]
            self.rec_pw_max = [9, 9, 6]
            self.nrec_ext = 25
            self.n_decay = [10, 3, -1]
            self.angz_cut1 = 1E-6
            self.tr_cut1 = 1E-3
            self.ai_cut1 = 1E-3
            self.angz_cut2 = 1E-4
            self.tr_cut2 = 1E-3
            self.ai_cut2 = 1E-3
        elif (self.nele <= self.nele_max[2]):
            self.n_shells = 2
            self.nterms = [-1,-1,-1,-1]
            if (self.nele > 5):
                self.nterms = [2, 2, 2, 2]
            self.nexc_max = [4, 4, 4, 4]
            self.nfrozen = [10, 10, 10, 10]
            self.nexc_rec = [7, 4, 4, 4]
            self.nrec_max = [10, 10, 10, 8]
            self.rec_pw_max = [9, 6, 8, 5]
            self.nrec_ext = 20
            self.n_decay = [10, 3, 3, -1]
            self.angz_cut1 = 1E-5
            self.tr_cut1 = 1E-3
            self.ai_cut1 = 1E-3
            self.angz_cut2 = 1E-4
            self.tr_cut2 = 1E-2
            self.ai_cut2 = 1E-2
        elif (self.nele <= self.nele_max[3]):
            self.n_shells = 3
            self.nterms = [-1,-1,-1,-1,-1]
            self.nexc_max = [4, 4, 4, 4, 0]
            self.nfrozen = [10, 10, 10, 10, 10]
            self.nexc_rec = [4, 4, 4, 4, 0]
            self.nrec_max = [10, 10, 10, 8, 0]
            self.rec_pw_max = [9, 8, 6, 5, 5]
            self.nrec_ext = 20
            self.n_decay = [10, 3, 3, -1, -1]
            self.angz_cut1 = 1E-5
            self.tr_cut1 = 1E-3
            self.ai_cut1 = 1E-3
            self.angz_cut2 = 1E-3
            self.tr_cut2 = 1E-2
            self.ai_cut2 = 1E-2
        else:
            raise 'ion with NELE >= %d not supported'%(self.nele_max[3])

        self.rec_pw_max = []
        for i in range(len(self.nrec_max)):
            self.rec_pw_max.append(self.nrec_max[i]-1)
            
        if (type(asym) == StringType):
            self.set_atom(asym, dir=dir)
            
        return

    
    def set_atom(self, asym, dir=''):
        self.asym = asym
        self.bfiles = {'en': dir+'%s%02db.en'%(asym, self.nele),
                       'tr': dir+'%s%02db.tr'%(asym, self.nele),
                       'ce': dir+'%s%02db.ce'%(asym, self.nele),
                       'rr': dir+'%s%02db.rr'%(asym, self.nele),
                       'ai': dir+'%s%02db.ai'%(asym, self.nele),
                       'ci': dir+'%s%02db.ci'%(asym, self.nele)}
        
        self.afiles = {'en': dir+'%s%02da.en'%(asym, self.nele),
                       'tr': dir+'%s%02da.tr'%(asym, self.nele),
                       'ce': dir+'%s%02da.ce'%(asym, self.nele),
                       'rr': dir+'%s%02da.rr'%(asym, self.nele),
                       'ai': dir+'%s%02da.ai'%(asym, self.nele),
                       'ci': dir+'%s%02da.ci'%(asym, self.nele)}

        return


    def add_excited(self, n0, n1, ibase):
        if (ibase == -1):
            base = self.grd_complex
        else:
            base = self.exc_complex[0].cgroup[ibase]
            
        n1.sort()
        
        cg = CGROUP('exc')
        k = len(self.exc_complex)
        bn = 'exc.%d'%(k)
        for i1 in n1:
            ex = COMPLEX('%s.%d.'%(bn,i1))
            if (ex.set_excited(n0, i1, base) == 0):
                if (self.nele in self.nele_sim and
                    k > self.n_shells):
                    ex.terms = ex.terms[:1]
                    ex.name = ex.name[:1]
                elif (len(ex.terms) == 3 and self.nele > 3):
                    if (k > 0):
                        ex.terms = ex.terms[:2]
                        ex.name = ex.name[:2]
                cg.add_complex(ex)

        self.exc_complex.append(cg)
        
        return
    

    def add_ionized(self, n0, ibase):
        if (ibase == -1):
            base = self.grd_complex
        else:
            base = self.exc_complex[0].cgroup[ibase]

        k = len(self.ion_complex.cgroup)
        ion = COMPLEX('%s.%d.'%('ion', k))
        ion.set_ionized(n0, base)
        if (len(ion.terms) == 3 and self.nele > 4):
            if (ibase == -1 and n0 != self.n_shells):
                ion.terms = ion.terms[:2]
                ion.name = ion.name[:2]
            elif (ibase > 0 and self.nele > 6):
                ion.terms = ion.terms[:2]
                ion.name = ion.name[:2]
        self.ion_complex.add_complex(ion)

        return


    def add_recombined(self, n0, ibase):
        base = self.ion_complex.cgroup[ibase]
        n0.sort()
        cg = self.exc_complex[ibase]
        if (ibase > self.n_shells and self.nele in self.nele_sim):
            bname = base.name[:1]
        else:
            bname = base.name
            if (len(bname) > 2 and ibase != 0):
                    bname = bname[:2]
        for i0 in n0:
            rec = COMPLEX('rec')
            if (rec.set_recombined(i0, bname) == 0):
                cg.add_complex(rec)
                
        return


    def run_tr_ce(self, a, b, c, d, tr = [-1], ce = 1):
        if (type(a) == StringType or type(a) == IntType):
            low = [a]
        else:
            low = a
        if (type(b) == StringType or type(b) == IntType):
            up = [b]
        else:
            up = b

	if (self.process['tr'] != 0):
            for m in tr:
                s = 'TR: %s -> %s %d'%(str(a), str(b), m)
                Print(s)
                TransitionTable(self.bfiles['tr'], low, up, m)

        if (ce != 0 and self.process['ce'] != 0):
            if (type(c) == StringType or type(c) == IntType):
                low = [c]
            else:
                low = c
            if (type(d) == StringType or type(d) == IntType):
                up = [d]
            else:
                up = d
            s = 'CE: %s -> %s'%(str(c), str(d))
            Print(s)
            CETable(self.bfiles['ce'], low, up)

        Reinit(radial = 1, excitation = 1, structure = 1)
        
        return


    def run_ci_rr(self, a, b, c, d, ci = 1, rr = 1):
        if (rr != 0 and self.process['rr'] != 0):
            if (type(c) == StringType or type(c) == IntType):
                low = [c]
            else:
                low = c
            if (type(d) == StringType or type(d) == IntType):
                up = [d]
            else:
                up = d
            s = 'RR: %s -> %s'%(c, d)
            Print(s)
            
            RRTable(self.bfiles['rr'], low, up)
            
            SetPEGrid(0)
            SetUsrPEGrid(0)
            SetRRTEGrid(0)
            
        if (ci != 0 and self.process['ci'] != 0):
            if (type(a) == StringType or type(a) == IntType):
                low = [a]
            else:
                low = a
            if (type(b) == StringType or type(b) == IntType):
                up = [b]
            else:
                up = b
            s = 'CI: %s -> %s'%(a, b)
            Print(s)
            CITable(self.bfiles['ci'], low, up)


        if (ci != 0 or rr != 0):
            Reinit(radial = 1, recombination = 1,
                   ionization = 1, structure = 1)
        
        return


    def run_ai(self, a, b, k):
	if (self.process['ai'] == 0):
	    return

        s = 'AI: %s -> %s %d'%(a, b, k)
        Print(s)
        if (type(a) == StringType or type(a) == IntType):
            low = [a]
        else:
            low = a
        if (type(b) == StringType or type(b) == IntType):
            up = [b]
        else:
            up = b

        AITable(self.bfiles['ai'], low, up, k)

        Reinit(radial = 1, recombination = 1, structure = 1)
        
        return
    

    def run_en(self):        
        SetAtom(self.asym)
        SetAngZCut(self.angz_cut0)
        
        g = self.grd_complex.name
        c = self.exc_complex[0].cgroup[0].name
        if (self.nele > 2):
            ConfigEnergy(0)
        if (self.nele < 3):
            OptimizeRadial(g)
        else:
            OptimizeRadial(g+c)
        if (self.nele > 2):
            ConfigEnergy(1)

        for ec in self.ecorrections:
            CorrectEnergy(ec[0], ec[1])

        # ground and the first excited complex are interacting
        Print('Structure: ground complex')
        Structure(self.bfiles['en'], g+c)
        
        Print('Structure: ionized complexes')
        c = self.ion_complex.cgroup
        for a in c:
            if (len(a.name) == 0):
                continue
            s = '    %s'%str(a.name)
            Print(s)
            Structure(self.bfiles['en'], a.name)

        Print('Structure: excited complexes')
        for i in range(len(self.exc_complex)):
            c = self.exc_complex[i].cgroup
            SetRecPWOptions(self.rec_pw_max[i])
            SetRecSpectator(self.nexc_max[i]+1, self.nfrozen[i])
            for j in range(len(c)):
                a = c[j]
                if (len(a.name) == 0):
                    continue
                if (i == 0 and j == 0):
                    continue
                if (a.nrec > 0):
                    s = '    (%s, %d)'%(str(a.name), a.nrec)
                    Print(s)
                    RecStates(self.bfiles['en'], a.name, a.nrec)
                else:
                    s = '    %s'%str(a.name)
                    Print(s)
                    Structure(self.bfiles['en'], a.name)

        return


    def run_tr_ce_all(self):
        g = self.grd_complex.name
        SetAngZCut(self.angz_cut0)
        SetTransitionCut(self.tr_cut0)
        self.run_tr_ce(g, g, g, g, tr=[-1,1,-2,2], ce=1)
        for i in range(self.n_shells):
            c = self.exc_complex[i].cgroup
            for j in range(len(c)):
                if (len(c[j].name) == 0):
                    continue
                a1 = g
                b = c[j].name
                nb = len(b)
                if (self.nterms[i] > 0 and nb > self.nterms[i]):
                    nb = self.nterms[i]
                    b = b[:nb]
                if (c[j].nrec == 0):
                    a2 = g
                else:
                    b = (b, c[j].nrec)
                    a2 = g[:1]
                if (i == 0 and (self.nele <= 2 or j == 0)):
                    tr = [-1, 1, -2, 2]
                    SetAngZCut(self.angz_cut0)
                    SetTransitionCut(self.tr_cut0)
                else:
                    tr = [-1]
                    SetAngZCut(self.angz_cut1)
                    SetTransitionCut(self.tr_cut1)
                nmax = c[j].nrec
                if (nmax == 0):
                    nmax = c[j].complex[-1][0]
                if (nmax <= self.nexc_rec[i]):
                    ce = 1
                else:
                    ce = 0
                if (ce != 0 or len(tr) != 0):
                    self.run_tr_ce(a1, b, a2, b, tr=tr, ce=ce)

        if (self.nele < 11):
            ce0 = 1
        else:
            ce0 = 0
        for i in range(len(self.exc_complex)):
            c = self.exc_complex[i].cgroup
            for j in range(len(c)):
                n2 = c[j].nrec
                b = c[j].name
                nb = len(b)
                if (nb == 0):
                    continue
                if (self.nterms[i] > 0 and nb > self.nterms[i]):
                    nb = self.nterms[i]
                    b = b[:nb]
                if (n2 == 0):
                    n2 = c[j].complex[-1][0]
                    SetAngZCut(self.angz_cut1)
                    SetTransitionCut(self.tr_cut1)
                else:
                    b = (b, n2)
                    SetAngZCut(self.angz_cut2)
                    SetTransitionCut(self.tr_cut2)
                if (i == 0 and j == 0):
                    ce = ce0
                else:
                    ce = 0
                for m in range(j+1):
                    a = c[m].name
                    na = len(a)
                    if (na == 0):
                        continue
                    if (self.nterms[i] > 0 and na > self.nterms[i]):
                        na = self.nterms[i]
                        a = a[:na]
                    nmax = c[m].nrec
                    if (nmax == 0):
                        nmax = c[m].complex[-1][0]
                    if (nmax <= self.n_decay[i]):
                        if (i == 0):
                            if (j == 0):
                                tr = [-1, 1]
                            else:
                                tr = [-1]
                        else:
                            if (n2 > self.nexc_max[i]):
                                tr = []
                            else:
                                tr = [-1]
                    else:
                        tr = []
                    if (c[m].nrec > 0):
                        a = (a, c[m].nrec)
                    if (ce != 0 or tr != []):
                        self.run_tr_ce(a, b, a[:1], b, tr=tr, ce=ce)

        for i in range(len(self.exc_complex)):
            c = self.exc_complex[i].cgroup
            d = self.exc_complex[0].cgroup
            for j in range(len(c)):
                b = c[j].name
                nb = len(b)
                if (nb == 0):
                    continue
                if (self.nterms[i] > 0 and nb > self.nterms[i]):
                    nb = self.nterms[i]
                    b = b[:nb]
                n2 = c[j].nrec
                if (n2 == 0):
                    n2 = c[j].complex[-1][0]
                    SetAngZCut(self.angz_cut1)
                    SetTransitionCut(self.tr_cut1)
                else:
                    b = (b, n2)
                    SetAngZCut(self.angz_cut2)
                    SetTransitionCut(self.tr_cut2)
                if (i == 0 and n2 <= self.n_decay[i]):
                    continue
                for m in range(len(d)):
                    a = d[m].name
                    na = len(a)
                    if (na == 0):
                        continue
                    if (self.nterms[0] > 0 and na > self.nterms[0]):
                        na = self.nterms[0]
                        a = a[:na]
                    n1 = d[m].nrec
                    if (n1 == 0):
                        n1 = d[m].complex[-1][0]
                        SetAngZCut(self.angz_cut1)
                        SetTransitionCut(self.tr_cut1)
                    else:
                        a = (a, n1)
                        SetAngZCut(self.angz_cut2)
                        SetTransitionCut(self.tr_cut2)
                    if (n1 != n2):
                        try:
                            nq = c[j].complex[-1][1]
                        except:
                            nq = 1
                        if (nq > 1):
                            n3 = n2
                        else:
                            try:
                                n3 = c[j].complex[-2][0]
                            except:
                                n3 = 0
                        if (n2 > self.nexc_max[i] or n3 != n1):
                            continue
                    self.run_tr_ce(a, b, a, b, tr=[-1], ce=0)
                        
        return

    
    def run_ci_rr_all(self):
        SetAngZCut(self.angz_cut0)
        g = self.grd_complex.name
        for i in range(self.n_shells):
            b = self.ion_complex.cgroup[i].name
            self.run_ci_rr(g, b, g, b, ci=1, rr=1)

        b = self.ion_complex.cgroup[0].name
        for i in range(len(self.exc_complex[0].cgroup)):
            a = self.exc_complex[0].cgroup[i].name
            na = len(a)
            if (na == 0):
                continue
            if (self.nterms[0] > 0 and na > self.nterms[0]):
                na = self.nterms[0]
                a = a[:na]
            n = self.exc_complex[0].cgroup[i].nrec
            if (n > 0):
                a = (a, n)
            b = b[:na]
            self.run_ci_rr(a, b, a, b, ci=1, rr=1)
             
        return

    
    def run_ai_all(self):
        for i in range(len(self.exc_complex)):
            c = self.exc_complex[i].cgroup
            if (i == 0):
                if (self.nele == 2):
                    continue
                if (self.nele == self.nele_max[self.n_shells-1]+1):
                    continue
                p = 1
            elif (i == self.n_shells-1 and
                  self.nele > self.nele_max[self.n_shells-1]+1):
                p = len(self.ion_complex.cgroup)
            else:
                p = 2
            for k in range(p):
                if (self.nele == self.nele_max[self.n_shells-1]+1):
                    t = k
                else:
                    if (p == 2 and k == 1):
                        t = self.n_shells
                    else:
                        t = k
                b = self.ion_complex.cgroup[t].name
                if (len(b) == 0):
                    continue
                if (p > 2 and t > 0):
                    if (i == t):
                        continue
                    nmax = self.ion_complex.cgroup[t].complex[-1][0]
                else:
                    nmax = 0
                for j in range(len(c)):
                    a = c[j].name
                    na = len(a)
                    if (na == 0):
                        continue
                    if (self.nterms[i] > 0 and na > self.nterms[i]):
                        na = self.nterms[i]
                        a = a[:na]
                    n = c[j].nrec
                    if (n > 0):
                        if (nmax > 0):
                            continue
                        a = (a, n)
                        SetAngZCut(self.angz_cut1)
                        SetAICut(self.ai_cut1)
                    else:
                        if (nmax > 0):
                            nmax1 = c[j].complex[-1][0]
                            if (nmax1 != nmax):
                                continue
                        SetAngZCut(self.angz_cut0)
                        SetAICut(self.ai_cut0)
                        
                    if (n > self.nrec_max[i] and
                        c[j].nrec_ext == 0):
                        continue
                    m = k*10+i
                    if (c[j].nrec_ext > 0):
                        m = n + (self.nrec_max[i]+1)*1000
                        m = -m
                    if (self.nele > 5):
                        b1 = b[:na]
                    else:
                        b1 = b
                    if (i == 0):
                        if (n > 0):
                            for c1 in a[0]:
                                for c2 in b1:
                                    self.run_ai(([c1], n), [c2], m)
                        else:
                            for c1 in a:
                                for c2 in b1:
                                    self.run_ai([c1], [c2], m)
                    else:
                        self.run_ai(a, b1, m)
                    
        return
    
    
    def print_table(self, type = [], v = 0):
        if (v != 0):
            MemENTable(self.bfiles['en'])

        for t in type:
            PrintTable(self.bfiles[t], self.afiles[t], v)
            
        return
  

    def set_configs(self):
        n = range(1, self.n_shells+1)
        n.reverse()
        nexc = self.nexc_max
        nrec = self.nrec_max
        
        for i in n:
            j = self.n_shells - i
            self.add_excited(i, range(nexc[j]+1), -1)
            self.add_ionized(i, -1)
            r = range(nexc[j]+1, nrec[j]+1)
            if (self.nrec_ext > nrec[j] and len(r) > 0):
                r.append(self.nrec_ext)
            self.add_recombined(r, j)
            try:
                self.exc_complex[j].cgroup[-1].nrec_ext = nrec[j]+1
            except:
                pass

        j = self.n_shells
        if (self.nele == self.nele_max[self.n_shells-1]+1):
            i = self.n_shells-1
            mrec = nrec[j+1]
            r = range(nexc[j]+1, mrec+1)
            if (self.nrec_ext > mrec and len(r) > 0):
                r.append(self.nrec_ext)
        else:
            i = self.n_shells
            mrec = nrec[j]
            r = range(nexc[j]+1, mrec+1)
            if (self.nrec_ext > mrec and len(r) > 0):
                r.append(self.nrec_ext)
        if (i > 0):
            self.add_excited(i, range(j+1, nexc[j]+1), 0)
            self.add_ionized(i, 0)
            self.add_recombined(r, j)
            try:
                self.exc_complex[j].cgroup[-1].nrec_ext = mrec+1
            except:
                pass

        if (self.nele > self.nele_max[self.n_shells-1]+1):
            i = self.n_shells
            j = self.n_shells+1
            self.add_excited(i, range(j+1, nexc[j]+1), 1)
            self.add_ionized(i, 1)
            r = range(nexc[j]+1, nrec[j]+1)
            if (self.nrec_ext > nrec[j] and len(r) > 0):
                r.append(self.nrec_ext)
            self.add_recombined(r, j)
            try:
                self.exc_complex[j].cgroup[-1].nrec_ext = nrec[j]+1
            except:
                pass
            
        # set up configurations
        c = self.grd_complex.name
        t = self.grd_complex.terms
        s = 'adding complex: %s, %s'%(str(c), str(t))
        Print(s)
        for i in range(len(t)):
            Config(t[i], group = c[i])
            
        for cg in self.exc_complex:
            for i in range(len(cg.cgroup)):
                if (cg.cgroup[i].nrec > 0):
                    continue
                c = cg.cgroup[i].name
                t = cg.cgroup[i].terms
                if (len(c) == 0):
                    continue
                s = 'adding complex: %s, %s'%(str(c), str(t))
                Print(s)
                for j in range(len(t)):
                    Config(t[j], group = c[j])

        cg = self.ion_complex.cgroup
        for i in range(len(cg)):
            c = cg[i].name
            if (len(c) == 0):
                continue
            t = cg[i].terms
            s = 'adding complex: %s, %s'%(str(c), str(t))
            Print(s)
            for j in range(len(t)):
                Config(t[j], group = c[j])
 
        return


    def run_fac(self):
        start = time.time()
        start0 = start

        # solve the structure
        s = 'EN...'
        Print(s)
        self.run_en()        
        stop = time.time()
        s = 'Done %10.3E s'%(stop-start)
        print s

        # transition rates and collisional excitation
        start = stop  
        s = 'TR and/or CE...'
        Print(s)  
        self.run_tr_ce_all() 
        stop = time.time()
        s = 'Done %10.3E s'%(stop-start)
        print s

        # CI and/or RR
        start = stop
        s = 'CI and/or RR...'
        Print(s)
        self.run_ci_rr_all()
        stop = time.time()
        s = 'Done %10.3E s'%(stop-start)
        print s 

        start = stop 
        s = 'AI...'
        Print(s)
        self.run_ai_all()
        stop = time.time()
        s = 'Done %10.3E s'%(stop-start)
        print s

        s = 'Total Time: %10.3E s'%(stop-start0)
        print s

        return


def atomic_data(nele, asym, iprint=1, dir='', **kw):
    if (type(nele) == IntType):
        n = [nele]
        if (kw.has_key('ecorrections')):
            kw['ecorrections'] = [kw['ecorrections']]
    else:
        n = nele

    for im in range(len(n)):
        m = n[im]
        atom = ATOM(m)
	if (kw.has_key('no_ce')):
	    atom.process['ce'] = not kw['no_ce']
	if (kw.has_key('no_tr')):
	    atom.process['tr'] = not kw['no_tr']
	if (kw.has_key('no_rr')):
	    atom.process['rr'] = not kw['no_rr']
	if (kw.has_key('no_ci')):
	    atom.process['ci'] = not kw['no_ci']
	if (kw.has_key('no_ai')):
	    atom.process['ai'] = not kw['no_ai']
        if (kw.has_key('nexc_max')):
            atom.nexc_max = kw['nexc_max']
        if (kw.has_key('nrec_max')):
            atom.nrec_max = kw['nrec_max']
        if (kw.has_key('nrec_ext')):
            atom.nrec_ext = kw['nrec_ext']
        if (kw.has_key('nexc_rec')):
            atom.nexc_rec = kw['nexc_rec']
        if (kw.has_key('n_decay')):
            atom.n_decay = kw['n_decay']
        if (kw.has_key('rec_pw_max')):
            atom.rec_pw_max = kw['rec_pw_max']
        if (kw.has_key('ecorrections')):
            atom.ecorrections = kw['ecorrections'][im]

        atom.set_configs()
        s = 'NELE = %d'%m
        Print(s)
        s = 'ATOM = %s'%asym
        z = ATOMICSYMBOL.index(asym.lower())
        s = '%s, Z = %d'%(s, z)
        Print(s)
        atom.set_atom(asym, dir=dir)
        atom.run_fac()
        if (iprint >= 0):
            atom.print_table(['en', 'tr', 'ce', 'rr', 'ai', 'ci'], iprint)
                
        Reinit(0)

    return
