# Sample DR calculation
import config
import fac
import dr

fac.SetAtom('Fe',26)
fac.SetAngZOptions(0, 1E-2, 1E-3)
config.closed('1s')

cfgs = ['2*1', '3*1']
gr   = ['n2', 'n3']
correlations = [['n2'],['n3']]

ngr = len(gr)

for i in range(ngr):
    config.config(cfgs[i], group = gr[i])   

fac.OptimizeRadial(gr)
print 'done optimize'
for c in correlations:
    fac.Structure(c)
print 'done structure'

prefix = 'Li'
fac.LevelTable('%s.lev'%prefix)
fac.TransitionTable(gr[0:1], gr, '%s.tr'%prefix, -1)

fac.SetTransitionCut(1E-3)
fac.SetAICut(1E-3)

rec_level = [0]
ground = [gr[0:1], gr[0:1]]
ai_group = [gr[0:1], gr[1:2]]
decay_group = [gr[0:1], gr[0:1]]

nopen = fac.DROpen(ai_group[0])
print nopen
ng1 = dr.getngrid(nopen, 60)
ng2 = filter(lambda x: x <= 20, ng1)
p = dr.getngrid([3, 4, 5], nopen[0]-1)
ng2[0:0] = p
ngrid = [ng1, ng2]
f = open('%s.ngrid'%prefix, 'w')
f.write(str(nopen))
f.write('\n')
for i in ngrid:
    print i
    f.write(str(i))
    f.write('\n')
f.close()

channels = [0, 1]
nmin = [2, 3]
nj0 = [11, 3]
max_kl = [12, 6]
nspec = [6, 6]

dr.dr(prefix, gr, rec_level, ground, ai_group, decay_group,
      channels, ngrid, nmin, nj0, max_kl, nspec, correlations)


for i in range(len(channels)):
    dr.drall(prefix, ngrid[i], channels[i], rate=[0.1, 1E4, 101],
             cross = [3.0, 0], strength=1)
    rt_file = dr.constructfn(prefix, 'rt', -1, channels[i])
    dr.sumrate(rt_file, prefix, ngrid[i], nopen, channels[i], 400)





