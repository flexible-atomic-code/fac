from pfac.crm import *
from pfac import const
from math import *
import sys
import time
import copy
import string
import cPickle
import biggles
import pprint

class TABLE:
    def __init__(self,
                 fname = '',
                 title='',
                 authors=[],
                 date='',
                 separator0 = '',
                 separator = ''):
        self.fname = fname
        self.title = title
        self.authors = authors
        if (time):
            self.date = date
        else:
            self.date = time.localtime()
        if (separator0):
            self.separator = separator
        else:
            self.separator0 = '='*72+'\n'
        if (separator):
            self.separator = separator
        else:
            self.separator = '-'*72+'\n'
            
        self.columns = []
        self.notes = []

    def add_column(self, **c):
        if (not c.has_key('label')):
            raise 'column must have a label'
        if (not c.has_key('format')):
            raise 'column must have a format'
        fmt = c['format']
        c['width'] = int(float(fmt[1:]))
        if (fmt[0] == 'A'):
            fmt = '%' + fmt[1:] + 's'
        elif (fmt[0] == 'I'):
            fmt = '%' + fmt[1:] + 'd'
        elif (fmt[0] == 'F'):
            fmt = '%' + fmt[1:] + 'f'
        elif (fmt[0] == 'E'):
            fmt = '%' + fmt[1:] + 'E'
        else:
            raise('Format unsupported')
        c['format0'] = fmt
        if (not c.has_key('units')):
            c['units'] = 'None'
        if (not c.has_key('description')):
            c['description'] = 'Column %d'%len(self.columns)
        if (not c.has_key('padding')):
            c['padding'] = ' '
        if (c.has_key('notes')):
            notes = c['notes']
            for a in notes:
                if (type(a) == type(1)):
                    k = a-1
                    if (k >= len(self.notes)):
                        raise 'Notes %d has not been defined'%(k+1)
                elif (type(a) != type(())):
                    raise 'New Notes must be in a tuple'
                k = a[0]-1
                t = a[1]
                if (k < len(self.notes)):
                    raise 'Notes %d already exists'%(k+1)
                elif (k > len(self.notes)):
                    raise 'Next notes must be %d'%(len(self.notes)+1)
                else:
                    self.notes.append(t)
                    c['notes'] = a[0]
                    
        self.columns.append(c)
            
    def open(self, mode, fname=''):
        if (fname):
            self.fname = fname
        self.file = open(self.fname, mode)

    def close(self):
        self.file.close()

    def write_header(self):
        f = self.file
        s = 'Title:   %s\n'%(self.title)
        f.write(s)
        a = ', '.join(self.authors)
        s = 'Authors: %s\n'%(a)
        f.write(s)
        f.write(self.separator0)
        s = 'byte-by-byte description of file: %s\n'%(self.fname)
        f.write(s)
        f.write(self.separator)
        s = '  %8s %6s %20s %10s %-s\n'%('Bytes', 'Format', 'Units',
                                         'Label', 'Explanation')
        f.write(s)
        f.write(self.separator)
        
        b0 = 1
        for i in range(len(self.columns)):
            c = self.columns[i]
            label = c['label']
            d = c['description']
            units = c['units']
            fmt = c['format']
            w = c['width']
            p = c['padding']
            b1 = b0 + w-1                
            s = '  %3d-%4d %6s %20s %10s %-s'%(b0, b1, fmt,
                                               units, label, d)
            if (c.has_key('notes')):
                s = s+'; (%d)'%(c['notes'])
            s = s + '\n'
            f.write(s)
            b0 = 1 + b1 + len(p)
        f.write(self.separator)
        if (len(self.notes) > 0):
            for i in range(len(self.notes)):
                s = 'Note (%d): '%(i+1)
                p = ' '*len(s)
                t = pad_text(self.notes[i], p)
                s = s+t+'\n'
                f.write(s)
            f.write(self.separator)
            
        self.dstart = f.tell()
        
    def write_raw(self, *data):
        if (len(data) != len(self.columns)):
            raise 'num. of data items does not match columns'
        s = ''
        for i in range(len(data)):
            fmt = self.columns[i]['format0']
            s = s + fmt%(data[i]) + self.columns[i]['padding']
        s = s + '\n'
        self.file.write(s)

    def read_header(self):
        f = self.file
        a = f.readline()
        a = a.split(':')
        self.title = a[1].strip()
        a = f.readline()
        a = a.split(':')
        self.author = a[1].split(',')
        a = f.readline()
        self.separator0 = a
        a = f.readline()
        a = f.readline()
        self.separator = a
        a = f.readline()
        a = f.readline()
        has_notes = 0
        while (1):
            a = f.readline()
            if (a[0] == self.separator[0]):
                break
            fmt = a[11:17].strip()
            units = a[18:38].strip()
            label = a[39:49].strip()
            d = a[50:-1].strip()
            d = d.split(';')
            if (len(d) > 1):
                has_notes = 1
            self.add_column(label=label,
                            units=units,
                            format=fmt,
                            description=d[0])
        if (has_notes):
            while (1):
                a = f.readline()
                if (a[0] == self.separator[0]):
                    break
        self.dstart = f.tell()
        
    def read_columns(self, k, filter='', start=0, stop=-1):
        nk = len(k)
        r = []
        for i in range(nk):
            r.append([])
        f = self.file
        nc = len(self.columns)
        n = range(nc)
        i = 0
        while (i < start):
            c = f.readline()
        
        while (1):
            if (stop >= 0):
                if (i >= stop):
                    break
            c = f.readline()
            if (not c):
                break
            c = c.split()
            for i in n:
                fmt = self.columns[i]['format'][0]
                if (fmt == 'I'):
                    c[i] = int(c[i])
                elif (fmt == 'F' or fmt == 'E'):
                    c[i] = float(c[i])
            
            if (filter):
                if (not eval(filter)):
                    continue
            for i in range(nk):
                r[i].append(c[k[i]])
        return r

    def convert2tex(self, tfn, k, filter='', start=0, stop=-1):
        f = open(tfn, 'w')
        self.open('r')
        self.read_header()
        d = self.read_columns(k, filter=filter, start=start, stop=stop)
        nd = len(d)
        for i in range(len(d[0])):
            s = ''
            for j in range(nd):
                fmt = self.columns[k[j]]['format0']
                s = s + fmt%(d[j][i])
                if (j == nd-1):
                    s = s + '\\\\'
                else:
                    s = s + ' & '
            s = s + '\n'
            f.write(s)
        self.close()
        f.close()
        
    def rewind(self):
        self.file.seek(self.dstart)
        
    
def pad_text(t, p):
    s = t.split('\n')
    a = ''
    for i in range(len(s)):
        a = a+s[i]
        if (i < len(s)-1):
            a = a + '\n' + p
    return a
    
def tabulate_states(dfile, neles, z = 26, dir='', pref='Fe', suffix='b'):
    tbl = TABLE(fname=dfile,
                title='Energy Levels for Z=%d'%z,
                authors=['M. F. Gu'],
                date=time.localtime())
    d = 'Num. of Electrons'
    tbl.add_column(label='NELE',
                   units='None',
                   description=d,
                   format='I2')
    d = 'Level Index'
    tbl.add_column(label='Index',
                   units='None',
                   description=d,
                   format='I6')
    d = 'Level Energy'
    notetext = 'For each ion, the true ground state energy is given, \n'
    notetext = notetext+'while the energies of excited states are given '
    notetext = notetext+'relative to the ground state.'
    tbl.add_column(label='Energy',
                   units='eV',
                   description=d,
                   format='E11.4',
                   notes=[(1, notetext)])
    d = 'Parity'
    tbl.add_column(label='P',
                   units='None',
                   description=d,
                   format='I1')
    d = 'Twice of Total Angular Momentum'
    tbl.add_column(label='2J',
                   units='None',
                   description=d,
                   format='I2')
    d = 'N Complex'
    tbl.add_column(label='NComplex',
                   units='None',
                   description=d,
                   format='A12')
    d = 'Non-relativistic Configuration'
    tbl.add_column(label='ConfigNR',
                   units='None',
                   description=d,
                   format='A12')
    d = 'Relativistic Configuration'
    tbl.add_column(label='ConfigR',
                   description=d,
                   format='A48')
    tbl.open('w')
    tbl.write_header()
    for k in neles:
        efile = dir+'%s%02d%s.en'%(pref, k, suffix)
        c = get_complexes(k)
        print '%d %s %s'%(k, efile, str(c))
        i = 0
        while (1):
            a = LevelInfor(efile, i)
            if (not a[3].strip() in c):
                break
            e = a[0]*const.Ryd_eV*2.0
            if (i == 0):
                e0 = e
            else:
                e = e-e0
            col = (k, i, e, a[1], a[2],
                   a[3].strip(), a[4].strip(), a[5].strip())
            tbl.write_raw(*col)
            i = i + 1
    tbl.close()
            
def tabulate_trates(dfile, neles, z=26, pref='Fe'):
    tbl = TABLE(fname=dfile,
                title='Total Ionization and Recombination Rate Coefficients',
                authors=['M. F. Gu'],
                date=time.localtime())
    d = 'Num. of Electrons'
    tbl.add_column(label='NELE', units='None',
                   description=d, format='I2')
    d = 'Temperature'
    tbl.add_column(label='Temp', units='[K]',
                   description=d, format='F4.2')
    d = 'Total DR rate coefficients'
    tbl.add_column(label='DR', units='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    d = 'Total DR Arnaud & Raymond'
    tbl.add_column(label='DR_AR', units='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    d = 'Total RR rate coefficients'
    tbl.add_column(label='RR', units='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    d = 'Total RR Arnaud & Raymond'
    tbl.add_column(label='RR_AR', units='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    d = 'Total DCI rate coefficients'
    tbl.add_column(label='CI', units='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    d = 'Total DCI Arnaud & Raymond'
    tbl.add_column(label='DCI_AR', units='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    d = 'Total EA rate coefficients'
    tbl.add_column(label='EA', units='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    d = 'Total EA Arnaud & Raymond'
    tbl.add_column(label='EA_AR', units='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    tbl.open('w')
    tbl.write_header()
    for k in neles:
        dir0 = '%s%02d/'%(pref, k)
        rates2 = cPickle.load(open(dir0+'rates2.sav', 'r'))
        rates3 = cPickle.load(open(dir0+'rates3.sav', 'r'))
        logt = rates2['logt']
        nt = len(logt)
        if (k == neles[0]):
            tdr = rates2['tdr'][1][1]
            trr = rates2['trr'][1][1]
            for i in range(nt):
                b = (10**logt[i])*const.kb
                (a0, a1, a2) = Recomb(z, k-1, b)
                (b0, b1, b2) = Ionis(z, k-1, b)
                tbl.write_raw(k-1, logt[i],
                              tdr[i], a2,
                              trr[i], a1,
                              0.0, b2,
                              0.0, b1)
        tdr = rates3['tdr'][2][1]
        trr = rates3['trr'][2][1]
        tci = rates2['tci'][0][1]
        tea = rates2['tea'][0][1]     
        for i in range(nt):
            b = (10**logt[i])*const.kb
            (a0, a1, a2) = Recomb(z, k, b)
            (b0, b1, b2) = Ionis(z, k, b)
            tbl.write_raw(k, logt[i],
                          tdr[i], a2,
                          trr[i], a1,
                          tci[i], b2,
                          tea[i], b1)
        if (k == neles[-1]):
            tci = rates3['tci'][1][1]
            tea = rates3['tea'][1][1]
            for i in range(nt):
                b = (10**logt[i])*const.kb
                (a0, a1, a2) = Recomb(z, k+1, b)
                (b0, b1, b2) = Ionis(z, k+1, b)
                tbl.write_raw(k+1, logt[i],
                              0.0, a2,
                              0.0, a1,
                              tci[i], b2,
                              tea[i], b1)
                
    tbl.close()    
    
def tabulate_rates(dfile, neles, z=26, pref='Fe'):
    tbl = TABLE(fname=dfile,
                title='Line Formation Rate Coefficients',
                authors=['M. F. Gu'],
                date=time.localtime())
    d = 'Num. of electrons'
    tbl.add_column(label='NELE',
                   units='None',
                   description=d,
                   format='I2')
    d = 'Level Index'
    tbl.add_column(label='Index',
                   units='None',
                   description=d,
                   format='I3')
    d = 'Temperature'
    tbl.add_column(label='Temp',
                   units='[K]',
                   description=d,
                   format='F4.2')
    units = 's^-1^'
    d = 'Total Depletion Rate'
    tbl.add_column(label = 'RT',
                   units=units,
                   description = d,
                   format = 'E8.2')
    units = '10^-10^cm^3^/s'
    d = 'Collisional Excitation' 
    tbl.add_column(label='CE',
                   units=units,
                   description=d,
                   format='E8.2')
    d = 'Resonance Excitation'
    tbl.add_column(label='RE',
                   units=units,
                   description=d,
                   format='E8.2')
    d = 'Radiative Recombination'
    tbl.add_column(label='RR',
                   units=units,
                   description=d,
                   format='E8.2')
    d = 'CE + n=3 Cascades'
    tbl.add_column(label='CE+CS3',
                   units=units,
                   description=d,
                   format='E8.2')
    d = 'CE + All Cascades'
    tbl.add_column(label='CE+CS',
                   units=units,
                   description=d,
                   format='E8.2')
    d = 'RE + All Cascades'
    tbl.add_column(label='RE+CS',
                   units=units,
                   description=d,
                   format='E8.2')
    d = 'DR + RR + n=3 Cascades'
    tbl.add_column(label='DR+RR+CS3',
                   units=units,
                   description=d,
                   format='E8.2')
    d = 'DR + RR + All Cascades'
    tbl.add_column(label='DR+RR+CS',
                   units=units,
                   description=d,
                   format='E8.2')
    d = 'Collisional Ionization'
    tbl.add_column(label='CI',
                   units=units,
                   description=d,
                   format='E8.2')

    tbl.open('w')
    tbl.write_header()
    for k in neles:
        dir0 = '%s%02d/'%(pref, k)
        rates3 = cPickle.load(open(dir0+'rates3.sav', 'r'))
        rates2 = cPickle.load(open(dir0+'rates2.sav', 'r'))
        rates1 = cPickle.load(open(dir0+'rates1.sav', 'r'))
        logt = rates2['logt']
        nt = len(logt)
        rt1 = rates1['rt']
        ce1 = rates1['ce']
        cs1 = rates1['cs']
        cs2 = rates2['cs']
        rr2 = rates2['rr']
        cs3 = rates3['cs']
        ci3 = rates3['ci']
        re3 = rates3['re']
        a1 = rates2['abund']
        b = []
        for i in range(len(rt1)):
            b.append([0.0]*nt)
        c = []
        for i in range(10):
            c.append(copy.deepcopy(b))
        for i in range(nt):
            p1 = a1[i][k-1]
            p2 = a1[i][k]
            p = p1/p2
            for m in range(len(rt1)):
                c[0][m][i] = rt1[m][2][i]
            for m in range(len(ce1)):
                n = ce1[m][1]
                c[1][n][i] = ce1[m][2][i]
            for m in range(len(cs1)):
                n = cs1[m][1]
                c[5][n][i] = cs1[m][2][i]
                c[4][n][i] = cs1[m][3][i]
            for m in range(len(cs2)):
                if (cs2[m][0] != k):
                    continue
                n = cs2[m][1]
                c[8][n][i] = (cs2[m][2][i]-c[5][n][i])/p
                c[7][n][i] = (cs2[m][3][i]-c[4][n][i])/p
                if (c[7][n][i] < 0):
                    c[7][n][i] = 0.0
                if (c[8][n][i] < c[7][n][i]):
                    c[8][n][i] = c[7][n][i]
            for m in range(len(rr2)):
                if (rr2[m][0] != k):
                    continue
                n = rr2[m][1]
                c[3][n][i] = rr2[m][2][i]
            for m in range(len(re3)):
                if (re3[m][0] != k):
                    continue
                n = re3[m][1]
                c[2][n][i] = re3[m][2][i]
            for m in range(len(ci3)):
                if (ci3[m][0] != k):
                    continue
                n = ci3[m][1]
                c[9][n][i] = ci3[m][2][i]
            for m in range(len(cs3)):
                if (cs3[m][0] != k):
                    continue
                n = cs3[m][1]
                c[6][n][i] = (cs3[m][2][i]-c[5][n][i]-p*c[8][n][i])
                if (c[6][n][i] < 0.0):
                    c[6][n][i] = 0.0
        for m in range(len(ce1)):
            for i in range(nt):
                col = [ce1[m][0], ce1[m][1]]
                col.append(logt[i])
                for j in range(10):
                    if (j == 4 or j == 5):
                        c[j][m][i] = c[j][m][i]+c[1][m][i]
                    elif (j == 6):
                        c[j][m][i] = c[j][m][i]+c[2][m][i]
                    elif (j == 7 or j == 8):
                        c[j][m][i] = c[j][m][i]+c[3][m][i]
                    col.append(c[j][m][i])
                col = tuple(col)
                tbl.write_raw(*col)
        
    tbl.close()
    
        
def write_trates(f, r, header, nele):
    s = '# %s\n'%(header)
    f.write(s)
    for i in range(len(r)):
        if (r[i][0] == nele):
            continue
        s = '%2d '%(r[i][0])
        for a in r[i][1]:
            s = s + '%9.3E '%(a)
        s = s[:-1] + '\n'
        f.write(s)
    f.write('\n')

def write_rates(f, r, header, nele):
    s = '# %s\n'%(header)
    f.write(s)
    for i in range(len(r)):
        if (r[i][0] == nele):
            continue
        for j in range(2, len(r[i])):
            s = '%2d %4d  '%(r[i][0], r[i][1])
            for a in r[i][j]:
                s = s + '%9.3E '%(a)
            s = s[:-1] + '\n'
            f.write(s)
    f.write('\n')
    

def save_rates(rates, sfile, dfile, **kwd):
    rates.update(kwd)
    rates['tdr'] = copy.deepcopy(rates['tdc'])
    for i in range(1, len(rates['tdc'])):
        rates['tdr'][i][1] = map(lambda x,y: x-y,
                                 rates['tdc'][i][1],
                                 rates['tre'][i-1][1])
    f = open(sfile, 'w')
    cPickle.dump(rates, f)
    f.close()

    f = open(dfile, 'w')
    if (rates.has_key('temp')):
        temp = rates['temp']
    else:
        temp = []
    if (rates.has_key('logt')):
        logt = rates['logt']
    else:
        logt = []
    if (rates.has_key('abund')):
        abund = rates['abund']
    else:
        abund = []
        
    neles = map(lambda x:x[0], rates['tdc'])
    nt = len(temp)
    s = '#   LogT(K)  Temp(eV)   '
    for k in neles:
        s = s + 'NELE=%-5d '%(k)
    s = s[:-1] + '\n'
    f.write(s)
    for i in range(nt):
        s = '%2d '%i
        if (logt):
            s = s + '%7.3f  '%(logt[i])
        s = s + '%10.4E '%(temp[i])
        if (abund):
            for k in neles:
                s = s + '%10.4E '%(abund[i][k])
        s = s[:-1]+'\n'
        f.write(s)
    f.write('\n')
    
    if (rates.has_key('tdc')):
        write_trates(f, rates['tdc'],
                     'Total Dielectronic Capture', neles[0])
    if (rates.has_key('tdr')):
        write_trates(f, rates['tdr'],
                     'Total Dielectronic Recombination', neles[0])
    if (rates.has_key('trr')):
        write_trates(f, rates['trr'],
                     'Total Radiative Recombination', neles[0])
    if (rates.has_key('tea')):
        write_trates(f, rates['tea'],
                     'Total Excitation Autoionization', neles[-1])
    if (rates.has_key('tci')):
        write_trates(f, rates['tci'],
                     'Total Direct Ionization', neles[-1])
    if (rates.has_key('rt')):
        write_rates(f, rates['rt'],
                    'Total Depletion Rate', -1)
    if (rates.has_key('cs')):
        write_rates(f, rates['cs'],
                    'Radiative Cascades', -1)
    if (rates.has_key('ce')):
        write_rates(f, rates['ce'],
                    'Direct Excitation', -1)
    if (rates.has_key('re')):
        write_rates(f, rates['re'],
                    'Resonance Excitation', neles[-1])
    if (rates.has_key('ci')):
        write_rates(f, rates['ci'],
                    'Direct Ionization', neles[-1])
    if (rates.has_key('rr')):
        write_rates(f, rates['rr'],
                    'Radiative Recombination', neles[0])

    f.close()
    
    
def read_rates(nt, nd, nele, pref='Fe', dir='', nion=2):
    c = get_complexes(nele)
    complexes = [c[1]]
    if (nion > 1):
        c = get_complexes(nele-1)
        complexes.append(c[1])
        if (nion > 2):
            c = get_complexes(nele+1)
            complexes.append(c[1])
    re = []
    ci = []
    rr = []
    cs = []
    ce = []
    rt = []
    tdc = []
    tre = []
    trr = []
    tpi = []
    tci = []
    tea = []
    for t in range(nt):
        for d in range(nd):
            rt_file = dir+'%s%02da_t%02dd%di%d.rt'%(pref, nele, t, d, nion)
            f = open(rt_file, 'r')
            ilev = -1
            den = [-1]*3
            ire = 0
            ici = 0
            irr = 0
            ics = 0
            ice = 0
            irt = 0
            itot = 0
            while (1):
                a = f.readline()
                if (not a):
                    break
                if (len(a) < 4):
                    continue
                if (a[:4] == 'NELE'):
                    a = string.split(a)
                    nel = int(a[2])
                elif (a[:4] == 'ILEV'):
                    a = string.split(a)
                    ilev = int(a[2])
                    nilev = 1
                    r3 = 0.0
                elif (a[:4] == ' Sum'):
                    a = string.split(a)
                    b = []
                    for c in a[1:-1]:
                        b.append(float(c))
                    if (den[0] > 0):
                        b[0] = b[0]/den[0]
                        b[1] = b[1]/den[0]
                    if (den[1] > 0):
                        b[2] = b[2]/den[1]
                    if (den[2] > 0):
                        b[3] = b[3]/den[2]
                        b[4] = b[4]/den[2]
                        b[5] = b[5]/den[2]
                    if (t == 0 and d == 0):
                        trr.append([nel, [b[0]]])
                        tdc.append([nel, [b[1]]])
                        tre.append([nel, [b[2]]])
                        tpi.append([nel, [b[3]]])
                        tea.append([nel, [b[4]]])
                        tci.append([nel, [b[5]]])
                    else:
                        trr[itot][1].append(b[0])
                        tdc[itot][1].append(b[1])
                        tre[itot][1].append(b[2])
                        tpi[itot][1].append(b[3])
                        tea[itot][1].append(b[4])
                        tci[itot][1].append(b[5])
                    itot = itot + 1
                elif (ilev >= 0):
                    if (a[:4] == 'Dens'):
                        a = string.split(a)
                        d0 = float(a[2])
                        tp = 0.0
                    elif (a[2:4] == '-1'):
                        a = string.split(a)
                        den[0] = float(a[1])
                        if (not den[0]):
                            continue
                        b = float(a[4])
                        if (b > 0):
                            tp = tp + b
                            b = b/den[0]
                            if (t == 0 and d == 0):
                                rr.append([nel, ilev, [b]])
                            else:
                                rr[irr][2].append(b)
                            irr = irr + 1
                    elif (a[2:4] == '-2'):
                        a = string.split(a)
                        den[1] = float(a[1])
                        if (not den[1]):
                            continue
                        b = float(a[2])
                        if (b > 0):
                            tp = tp + b
                            b = b/den[1]
                            r3 = r3/den[1]
                            if (t == 0 and d == 0):
                                cs.append([nel, ilev, [b], [r3]])
                            else:
                                cs[ics][2].append(b)
                                cs[ics][3].append(r3)
                            ics = ics + 1
                        b = float(a[3])
                        if (b > 0):
                            tp = tp + b
                            b = b/den[1]
                            if (t == 0 and d == 0):
                                ce.append([nel, ilev, [b]])
                            else:
                                ce[ice][2].append(b)
                            ice = ice + 1
                        b = float(a[5])
                        if (b > 0):
                            tp = tp + b
                            b = b/den[1]
                            if (t == 0 and d == 0):
                                re.append([nel, ilev, [b]])
                            else:
                                re[ire][2].append(b)
                            ire = ire + 1
                    elif (a[2:4] == '-3'):
                        a = string.split(a)
                        den[2] = float(a[1])
                        if (den[2] > 0):
                            b = float(a[6])
                            if (b > 0):
                                tp = tp + b
                                b = b/den[2]
                                if (t == 0 and d == 0):
                                    ci.append([nel, ilev, [b]])
                                else:
                                    ci[ici][2].append(b)
                                ici = ici + 1
                        tp = tp/d0
                        if (t == 0 and d == 0):
                            rt.append([nel, ilev, [tp]])
                        else:
                            rt[irt][2].append(tp)
                        irt = irt + 1
                    else:
                        c = a[72:-1]
                        if (len(c) > 1):
                            c = string.strip(c)
                            if c in complexes:
                                a = string.split(a[:72])
                                b = float(a[2])
                                if (b > 0):
                                    r3 = r3 + b
                                
                            
    def compare(x, y):
        if (x[0] < y[0]):
            return -1
        elif (x[0] > y[0]):
            return 1
        else:
            if (x[1] < y[1]):
                return -1
            elif (x[1] > y[1]):
                return 1
            else:
                return 0
    cs.sort(compare)
    ce.sort(compare)
    re.sort(compare)
    ci.sort(compare)
    rr.sort(compare)
    rt.sort(compare)
    return {'tdc': tdc,
            'tre': tre,
            'trr': trr,
            'tpi': tpi,
            'tea': tea,
            'tci': tci,
            'rt':  rt,
            'cs':  cs,
            'ce':  ce,
            'rr':  rr,
            're':  re,
            'ci':  ci}

                    
def get_tgrid(z, nele, dt = 0.15, amin = 1E-2, limits=[]):
    if (len(limits) == 0):
        limits = [5.0, 8.0]
        
    (tmax,a) = MaxAbund(z, nele)
    amax = a[nele-1:nele+1]
    logtm = log10(tmax/const.kb)
    logtm = limits[0]+int((logtm-limits[0])/dt)*dt
    a0 = 1.0
    logt = [logtm]
    t = [const.kb*10**(logtm)]
    a = FracAbund(z, t[0])
    ab = [a]
    while (a0 > amin and logt[-1] < limits[1]):
        logt0 = logt[-1] + dt
        t0 = const.kb*10**(logt0)
        t.append(t0)
        logt.append(logt0)
        a = FracAbund(z, t0)
        ab.append(a)
        a0 = max(map(lambda x,y:x/y, a[nele-1:nele+1], amax))
        
    a0 = 1.0
    while (a0 > amin and logt[0] > limits[0]):
        logt0 = logt[0] - dt
        t0 = const.kb*10**(logt0)
        t.insert(0, t0)
        logt.insert(0, logt0)
        a = FracAbund(z, t0)
        ab.insert(0, a)
        a0 = max(map(lambda x,y:x/y, a[nele-1:nele+1], amax))
        
    return (t, logt, ab)

    
def get_complexes(nelectrons):
    n = 1
    nele = nelectrons
    g = []
    while (nele > 0):
        nqm = 2*n*n
        if (nele < nqm):
            g.append((n,nele))
        else:
            g.append((n,nqm))
        nele = nele - nqm
        n = n + 1
    c0 = ''
    c1 = ''
    for a in g:
        if (a != g[-1]):
            c0 = c0 + '%d*%d '%a
            c1 = c1 + '%d*%d '%a
        else:
            c0 = c0 + '%d*%d'%a
            if (a[1] == 1):
                c1 = c1 + '%d*%d'%(a[0]+1,1)
            else:
                c1 = c1 + '%d*%d %d*%d'%(a[0], a[1]-1, a[0]+1, 1)
    return (c0, c1)

def spectrum(neles, temp, den, population, pref,
             suf = 'b', dir0 = '', dir1= '', nion = 3,
             dist = 0, cascade = 0, rrc = 0):
    for k in neles:
        rate = get_complexes(k)
        if (nion > 1):
            rate = rate + get_complexes(k-1)
            if (nion > 2):
                rate = rate + get_complexes(k+1)
                
        print 'NELE = %d'%k
        f1 = '%s%02d%s'%(pref, k-1, suf)
        f2 = '%s%02d%s'%(pref, k, suf)
        f3 = '%s%02d%s'%(pref, k+1, suf)
        AddIon(k, 0.0, dir0+f2) 
        if (nion == 3):
            AddIon(k+1, 0.0, dir0+f3)
        if (nion > 1):
            if (k > 1):
                SetBlocks(0.0, dir0+f1)
            else:
                SetBlocks(0.0)
        else:
            SetBlocks(-1.0)

        for i in range(len(temp)):
            p1 = population[i][k-1]
            p2 = population[i][k]
            p3 = population[i][k+1]
            p1 = p1/p2
            p3 = p3/p2
            p2 = 1.0
            print 'Temp = %10.3E'%(temp[i])
            print 'Abund: %10.3E %10.3E %10.3E'%(p1, p2, p3)

            SetEleDist(dist, temp[i], -1.0, -1.0)
            print 'CE rates...'
            SetCERates(1)
            print 'TR rates...'
            SetTRRates(0)
            if (nion > 1):
                print 'RR rates...'
                SetRRRates(0)
                print 'CI rates...'
                SetCIRates(0)
                print 'AI rates...'
                SetAIRates(1)
                SetAbund(k-1, p1)
            SetAbund(k, p2)
            if (nion == 3):
                SetAbund(k+1, p3)
                
            for d in range(len(den)):
                print 'Density = %10.3E'%den[d]
                SetEleDensity(den[d])
                SetCascade(cascade, 1E-4)
                print 'Init blocks...'
                InitBlocks()
                s = 't%02dd%di%d'%(i, d, nion)
                rt_file = dir1+'%s_%s.rt'%(f2,s)
                sp_file = dir1+'%s_%s.sp'%(f2,s)
                rt_afile = dir1+'%sa_%s.rt'%(f2[:-1],s)
                sp_afile = dir1+'%sa_%s.sp'%(f2[:-1],s)

                LevelPopulation()
                Cascade()

                rt = (rt_file,)+rate
                RateTable(*rt)
                SpecTable(sp_file, rrc)
                PrintTable(rt_file, rt_afile, 1)
                PrintTable(sp_file, sp_afile, 1)
                sys.stdout.flush()
                ReinitCRM(2)
            ReinitCRM(1)
        ReinitCRM()

def read_lines(file, nele):
    k = -1
    s = []
    for line in open(file, 'r').readlines():
        a = string.split(line)
        if (len(a) != 5):
            continue
        if (int(a[0]) != nele):
            continue
        e = float(a[3])
        i = int(a[1])
        j = int(a[2])
        s0 = float(a[4])
        s.append([i, j, e, s0])
                
    return s

def id_lines(w, s, name, eps=''):
    n = len(w)
    biggles.configure('screen', 'width', 800)
    biggles.configure('screen', 'height', 500)
    p = biggles.FramedPlot()
    p.aspect_ratio = 0.7
    xmin = min(w)
    xmax = max(w)
    for i in range(n):
        q1 = (w[i], 0)
        q2 = (w[i], s[i])
        l = biggles.Line(q1, q2)
        p.add(l)
        l = biggles.Label(q2[0], q2[1], '%d: %s'%(i,name[i]),
                          color='red', textangle=90)
        p.add(l)
    p.add(biggles.Line((xmin,0), (xmax,0), color='red'))
    if (eps):
        p.write_eps(eps)
    else:
        p.show()
