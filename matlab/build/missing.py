#!/bin/python
from modeller import *
from modeller.automodel import *

log.verbose();
e=environ()

code0='G120'
code1='G120-ALL'

m0=model(e,file=code0)
m1=model(e,file=code1)

aln=alignment(e)
aln.append_model(m0, align_codes=code0)
aln.append_model(m1, align_codes=code1)

# rough align with default parameters; this is easy because the sequences are identical where coordinates are defined
aln.salign();

aln.write(file=code0+'.ali')

#model=loopmodel(e, alnfile=code0+'.ali', knowns=code0, sequence=code1)
model=automodel(e, alnfile=code0+'.ali', knowns=code0, sequence=code1)
model.starting_model=1
model.ending_model=10

#model.loop.starting_model=1
#model.loop.ending_model  =1

model.make()
