from pfac.spm import *
from pfac import const

nt = 9
logt = range(nt)
logt = map(lambda x:x*0.18+6.3, logt)
temp = map(lambda x:(10**x)*const.kb, logt)

z = 26
population = []
for i in range(len(temp)):
    a = FracAbund(z, temp[i])
    population.append(a)

neles = range(2, 12)
den = [1.0]

trans = [(3,2), (4,2)]
spectrum(neles[:-1], temp, den, population,
         'Fe', dir='../', nion=3)
spectrum(neles[-1:], temp, den, population,
         'Fe', dir='../', nion=2)
spectrum(neles, temp, den, population,
         'Fe', dir='../', nion=1)
