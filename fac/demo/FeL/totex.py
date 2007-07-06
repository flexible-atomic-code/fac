from pfac.table import *

tbl = TABLE('rates.tbl')
tbl.convert2tex('rates.tex', [1,2,4,5,8,9,11,12], filter='c[0]==10')

tbl = TABLE('trates.tbl')
tbl.convert2tex('trates.tex', range(1,8), filter='c[0]==10 or c[0]==11')
