pA=0.05
pBA=0.99
pNotBNotA=0.95

pBNotA=1-pNotBNotA


pB=pBA*pA+pBNotA*(1-pA)

pAB=pBA*pA/pB
pANotB=(1-pBA)*pA/(1-pB)


print "pB=",pB
print "pAB=",pAB
print "pANotB=",pANotB

