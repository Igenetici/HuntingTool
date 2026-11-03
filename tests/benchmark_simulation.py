"""Accuracy on simulated heterozygous genotypes at 2%, 5% and 8% per-base error.
Run: python tests/benchmark_simulation.py"""
import random, sys
import os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from collections import Counter
from huntingtool.structure import parse_read
from huntingtool.caller import call_alleles
from tests.simulate import reads
rng=random.Random(3)
cfgs=['canonical']*6+['CAA-loss','CAACAG-dup','CCA-loss','CAACAGCCGCCA-loss']
res=Counter()
fails=[]
for err in (0.02,0.05,0.08):
  for t in range(60):
    a1=(rng.randint(12,26),rng.choice(cfgs),rng.choice([7,7,7,10]),1.0)
    a2=(rng.choice([rng.randint(12,30),rng.randint(36,90)]),'canonical',rng.choice([7,7,10]),rng.uniform(0.3,1))
    if abs(a1[0]-a2[0])<2 and a1[2]==a2[2]: continue
    hap=Counter(); adj={}
    for s,q in reads([a1,a2],1500,err,seed=t):
        st=parse_read(s)
        if st and st.status=='complete':
            hap[(st.config,st.ccg,st.cag)]+=1; adj.setdefault(st.config,st.adj)
    al,notes=call_alleles(hap,adj)
    got=sorted((a.cag,a.config,a.ccg) for a in al)
    exp=sorted((x[0],x[1],x[2]) for x in (a1,a2))
    ok=got==exp
    okcag=sorted(a.cag for a in al)==sorted(x[0] for x in (a1,a2))
    res[(err,'exact')]+=ok; res[(err,'cag')]+=okcag; res[(err,'n')]+=1
    if not ok and len(fails)<12: fails.append((err,exp,got,notes[:2]))
for e in (0.02,0.05,0.08): print(e, res[(e,'exact')],'/',res[(e,'n')],'exact; cag ok',res[(e,'cag')])
for f in fails: print(f)
