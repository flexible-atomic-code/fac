"""
This module implements an ASCII TABLE class, conforming to the format
of machine readable tables in ApJ.
"""

import time

def pad_text(t, p):
    s = t.split('\n')
    a = ''
    for i in range(len(s)):
        a = a+s[i]
        if (i < len(s)-1):
            a = a + '\n' + p
    return p+a
    
class TABLE:
    def __init__(self,
                 fname='',
                 title='',
                 authors=[],
                 date='',
                 separator0='',
                 separator=''):
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
        if (not c.has_key('unit')):
            c['unit'] = 'None'
        if (not c.has_key('description')):
            c['description'] = 'Column %d'%len(self.columns)
        if (not c.has_key('padding')):
            c['padding'] = ' '
                    
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
        s = '  %8s %6s %20s %10s  %-s\n'%('Bytes', 'Format', 'Units',
                                          'Label', 'Explanation')
        f.write(s)
        f.write(self.separator)
        
        b0 = 1
        for i in range(len(self.columns)):
            c = self.columns[i]
            label = c['label']
            d = c['description']
            if (c.has_key('note')):
                d = '*'+d
            else:
                d = ' '+d
            unit = c['unit']
            fmt = c['format']
            w = c['width']
            p = c['padding']
            b1 = b0 + w-1                
            s = '  %3d-%4d %6s %20s %10s %-s'%(b0, b1, fmt,
                                               unit, label, d)
            s = s + '\n'
            f.write(s)
            b0 = 1 + b1 + len(p)
        f.write(self.separator)
        s = ''
        for i in range(len(self.columns)):
            c = self.columns[i]
            if (c.has_key('note')):
                s = 'Note on %s:\n'%c['label']
                p = ' '*4
                t = pad_text(c['note'], p)
                s = s + t + '\n'
                f.write(s)
        if (s):
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
            unit = a[18:38].strip()
            label = a[39:49].strip()
            d = a[50:-1]
            if (d[0] == '*'):
                has_notes = 1
                self.add_column(label=label,
                                unit=unit,
                                format=fmt,
                                description=d[1:],
                                note=' ')
            else:
                self.add_column(label=label,
                                unit=unit,
                                format=fmt,
                                description=d[1:])
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
