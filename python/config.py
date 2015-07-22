#
#   FAC - Flexible Atomic Code
#   Copyright (C) 2001-2015 Ming Feng Gu
# 
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
# 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
# 
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import string
import fac

# spectroscopic notation of orbital angular momentums.
_orbital_symbols = {}
for i in range(len(fac.SPECSYMBOL)):
    _orbital_symbols[fac.SPECSYMBOL[i]] = i
_group_name = ''
_closed_shells = []

def avgconfig(configs):
    shells = string.split(configs)
    acfg = []
    for s in shells:
        n = len(s)
        i = n-1
        nq = 1.0
        while (i > 0):
            if (s[i] == '+' or s[i] == '-'):
                break
            try:
                a = float(s[i:n])
                nq = a
                i = i - 1
            except:
                break
        s = s[0:i+1]
        c = distribute(s)
        n = len(c)
        nq = nq/n
        for t in c:
            p = t[0]
            acfg.append((p[0], p[1], p[2], nq))
    fac.SetAvgConfig(acfg)

    
def closed(*configs):
    global _closed_shells

    if len(configs) == 0:
        _closed_shells = []
    for con in configs:
        shells = string.split(con)
        for s in shells:
            c = distribute(s)
            for t in c:
                p = t[0]
                q = 2 * p[1] + p[2] + 1
                _closed_shells.append((p[0], p[1], p[2], q))

    
# specify configurations.
def config(*configs, **group):
    global _group_name, _closed_shells

# If the group name is not specified, the one last specified is used.
    if group.has_key('group'):
        _group_name = group['group']

    if len(configs) == 0:
        if len(_closed_shells) == 0:
            fac.AddConfig(_group_name, [(1, 0, 1, 0)])
        else:
            fac.AddConfig(_group_name, _closed_shells)
        return
    
    for c in configs:
        shells = string.split(c)
        configurations = distribute(shells[0])
            
        if len(shells) > 1:
            for s in shells[1:]:
                sconf = distribute(s)
                new_conf = []
                for t in configurations:
                    for next_shell in sconf:
                        new_conf.append(t+next_shell)
                configurations = new_conf
        for t in configurations:
            fac.AddConfig(_group_name, _closed_shells + t)


def expand(c):
    shells = string.split(c)
    configurations = distribute(shells[0])
            
    if len(shells) > 1:
        for s in shells[1:]:
            sconf = distribute(s)
            new_conf = []
            for t in configurations:
                for next_shell in sconf:
                    new_conf.append(t+next_shell)
            configurations = new_conf

    return configurations


def cname(c):
    name = ''
    for a in c:
        if (a[2] > 0):
            j = '+'
        else:
            j = '-'
        k = fac.SPECSYMBOL[a[1]]
        s = '%d%s%s%d '%(a[0],k,j,a[3])
        name = name+s
        
    return name[:-1]


# For a given specified shell with wild casts, generate all possible
# distribution of electrons. eg., '2p2' -> ['2p-2', '2p-1 2p+1', '2p+2']
def distribute(s):
    
    # the first char must be the principle quantum number.
    len_s = len(s)
    for i0 in range(len_s):
        if not s[i0:i0+1].isdigit():
            break
    
    n = string.atoi(s[0:i0])

    # 2nd char or those in the [] is the orbital angular momentum.
    # * means any allowed value.
    if s[i0] == '*':
        l = range(n)
    elif s[i0] == '[':
        i0 = i0+1
        i1 = i0
        while (s[i1] != ']'):
            i1 = i1+1
        b = string.split(s[i0:i1],',')
        l = []
        for a in b:
            if a.isdigit():
                l.append(string.atoi(a))
            else:
                l.append(_orbital_symbols[a])
        i0 = i1
    else:
        l = [_orbital_symbols[s[i0]]]

    # prepare to calculate the half of maximum occupation number.
    half_q = 0
    len_l = len(l)
    for m in l:
        half_q = half_q + m

    # 3rd char may indicate whether j = l+/-0.5. If it's not +/-, j may be
    # either value. The remaining chars (4th on if the 3rd is +/-, 3rd on
    # otherwise) indicate the occupation number.

    i0 = i0+1
    if len_s == i0:
        j = [-1, 1]
        q = 1
        half_q = 2*half_q + len_l
    elif len_s > i0:
        i1 = i0 + 1
        if s[i0] == '+':
            j = [1]
            if len_s > i1:
                q = string.atoi(s[i1:])
            else:
                q = 1
            half_q = half_q + len_l
        elif s[i0] == '-':
            j = [-1]
            if len_s > i1:
                q = string.atoi(s[i1:])
            else:
                q = 1
        else:
            j = [-1,1]
            q = string.atoi(s[i0:])
            half_q = 2*half_q + len_l
            
    # if the occupation number is more than half of the maximum,
    # convert it to the hole states.
    hole = 0 
    if q > half_q:
        q = 2*half_q - q
        hole = 1
        
    shells = []
    for m in l:
        if m == 0:
            mj = 1
            maxq = 2*m + mj + 1
            s = []
            qq = range(min(maxq,q)+1)
            if not hole:
                qq.reverse()
            for mq in qq:
                if hole:
                    mq = maxq - mq
                s.append([(n, m, mj, mq)])
            shells.append(s)
        else:
            for mj in j:
                maxq = 2*m + mj +1
                s = []
                qq = range(min(maxq,q)+1)
                if not hole:
                    qq.reverse()
                for mq in qq:
                    if hole:
                        mq = maxq - mq
                    s.append([(n, m, mj, mq)])
                shells.append(s)

    if hole:
        q = 2*half_q - q
    all_configs = []
    configurations = shells[0]
    tq = []
    for t in configurations:
        tq.append(t[0][3])
    if len(shells) > 1:
        for s in shells[1:]:
            new_conf = []
            ntq = []
            for i in range(len(configurations)):
                t = configurations[i]
                for next_shell in s:
                    x = tq[i] + next_shell[0][3]
                    if x > q:
                        continue
                    new_conf.append(t+next_shell)
                    ntq.append(x)
            configurations = new_conf
            tq = ntq

    for t in configurations:
        qq = 0
        tt = []
        for s in t:
            if s[3] > 0:
                qq = qq + s[3]
                tt.append(s)
        if qq == q:            
            all_configs.append(tt)
    return all_configs
            
