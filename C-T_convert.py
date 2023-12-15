#!/usr/bin/env python
# coding: utf-8

# In[1]:


# read reference,
reference = "FMR1_BisTemplate.fa"
from Bio import SeqIO
ref_data = SeqIO.parse(open(reference),'fasta')
for fasta in ref_data:
    name, ref = fasta.id, str(fasta.seq)
    ref_split = ref.strip()
    print(ref)


# In[2]:


outfile = open("C-T_convert.fa", 'w')
cCount=0
for i in range(len(ref_split)):
    if ref_split[i] == "C":
        cCount += 1
        seq = list(ref_split)
        seq[i] = "T"
        ref_mod = ''.join(seq)
        print(">" + str(cCount), file=outfile)
        print(ref_mod, file=outfile)
    else:
        pass
outfile.close()


# In[ ]:




