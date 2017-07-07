import os
text=open('text.tm').read()
pbs=open('BEM.tpbs').read()
zr=[2.095,2.105,100]
gap=(zr[1]-zr[0])/zr[2]
div=10
for i in range(zr[2]/div):
    open('text{0}.m'.format(i),'w').write(text.format(zr[0]+gap*i*div,zr[0]+gap*(i+1)*div,div))
    open('BEM{0}.pbs'.format(i),'w').write(pbs.format(i))

os.system('for((i=0;i<{0};i++));do qsub BEM$i.pbs;done'.format(zr[2]/div))
