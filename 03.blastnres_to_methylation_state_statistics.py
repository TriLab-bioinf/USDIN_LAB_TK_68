#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re
import pandas as pd
filename = "ref_to_all_inserts_q0.txt"
data = open(filename, 'r')


# In[2]:


# collect the alignments
isC = []
seqs = ''
alignDict = {}
scoreDict = {}
for line in data:
    if line.isspace():
        continue
    elif line[:1] == ">":
        isC = []
        seqs = ''
        id = re.sub(">","",line.strip().split()[0])
    elif line[1:6] == "Score":
        score = line.strip().split()[2]
        scoreDict[id] = score
    elif line[:6] == "Query ":
        refs = line.strip().split()[2]
        refs = refs.upper()
        for i in range(len(refs)):
            if refs[i] == "C":
                isC.append(i)
            else: pass
    #elif re.search(id,line):
    elif line[:6] == "Sbjct ":
        maps = line.strip().split()[2]
        for m in isC:
            seqs = seqs + maps[m]
        if id in alignDict:
            alignDict[id] = alignDict[id] + seqs
        else:
            alignDict[id] = seqs
    else: pass
data.close()


# In[3]:


## quality control
df1 = pd.DataFrame(scoreDict.items(), columns=['ID', 'score'])
df2 = pd.DataFrame(alignDict.items(), columns=['ID', 'sequence'])
df_join = pd.merge(df1,df2,on='ID')


# In[4]:


df_join['numC'] = ''
df_join['numT'] = ''
df_join['numdash'] = ''
df_join['Cpcnt'] = ''
df_join['class'] = ''
df_join['Phred'] = ''


# In[5]:


# disable chained assignments
pd.options.mode.chained_assignment = None 
for i in range(df_join.shape[0]):
    df_join['numC'].iloc[i] = len(re.findall("[C]", df_join["sequence"].iloc[i]))
    df_join['numT'].iloc[i] = len(re.findall("T", df_join["sequence"].iloc[i]))
    df_join['numdash'].iloc[i] = len(re.findall("A|G|-", df_join["sequence"].iloc[i]))
    df_join['Cpcnt'].iloc[i] = df_join['numC'].iloc[i]/(df_join['numC'].iloc[i]+df_join['numT'].iloc[i])
    df_join['Phred'].iloc[i] = float(re.sub('^.*_q', '', df_join['ID'].iloc[i]))
    if df_join['Phred'].iloc[i] < 20:
        df_join['class'].iloc[i] = '<20'
    elif df_join['Phred'].iloc[i] >= 20 and df_join['Phred'].iloc[i] <25:
        df_join['class'].iloc[i] = '20~25'
    elif df_join['Phred'].iloc[i] >= 25 and df_join['Phred'].iloc[i] <30:
        df_join['class'].iloc[i] = '25~30'
    elif df_join['Phred'].iloc[i] >= 30 and df_join['Phred'].iloc[i] <35:
        df_join['class'].iloc[i] = '30~35'
    elif df_join['Phred'].iloc[i] >= 35 and df_join['Phred'].iloc[i] <40:
        df_join['class'].iloc[i] = '35~40'
    else:
        df_join['class'].iloc[i] = '>40'


# In[6]:


df_join.score = df_join.score.astype(float)
df_join.Phred =df_join.Phred.astype(float)


# In[7]:


df_join


# In[8]:


# phred score vs number of dash
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="darkgrid")

sns.boxplot(x='class', y='numdash', data=df_join, order=["<20", "20~25", "25~30","30~35", "35~40", ">40"]).set(title='Phred score vs number of dash')
plt.show()


# In[9]:


# number of dash vs alignment score
sns.set(style="darkgrid")

sns.boxplot(x='numdash', y='score', data=df_join)
plt.show()


# In[10]:


# phred score vs alignment score
sns.boxplot(x='class', y='score', data=df_join, order=["<20", "20~25", "25~30","30~35", "35~40", ">40"]).set(title='Phred score vs alignment socre')
plt.show()


# In[11]:


## correlation
sns.lmplot(x='Phred', y='score', data=df_join)
from scipy import stats
stats.pearsonr(df_join['Phred'], df_join['score'])


# In[12]:


df_join.to_csv('blastn_q0.res.csv')


# In[13]:


# filter sequences with number of dash<=1 and Phred>30
df_filter = df_join.loc[(df_join['numdash']<=1) & (df_join['Phred']>30)]
alignDict2 = {}
alignDict2 = {k: v for k, v in alignDict.items() if k in df_filter['ID'].values}


# In[14]:


lenQuery = 38
# set up lists to later calculate % methylation at each position
Cpcnt = [0] * lenQuery
tpcnt = [0] * lenQuery

# count the number of alleles with a given % of methylation in bins of 10
allPcnt = [0] * 10

# count the number of unique alleles
allType = {} 


# In[15]:


# then see if those positions are C or T in the test sequences
for id in sorted(alignDict2):
    checkC = alignDict2[id]
    if len(checkC) != lenQuery:
        checkC = checkC + ((lenQuery - len(checkC)) * "X")
    else: pass
    # collect the C or t value and update % methylation at each position
    cCount = ""
    for i in range(38):
        if checkC[i] == "C":
            cCount += "C"
            Cpcnt[i] +=1
        elif checkC[i] == "T":
            cCount += "t"
            tpcnt[i] +=1
        else:
            cCount += "-"
    # get the pcent methylation
    numC = cCount.count("C")
    numt = cCount.count("t")
    if numC + numt > 0:
        pcent = numC/(numC + numt)
        pcentIndex = int(10 * pcent)
        if pcentIndex == 10:
            pcentIndex = 9
        else: pass
        allPcnt[pcentIndex] +=1
        # check for allele uniqueness
        if cCount in allType:
            allType[cCount] +=1
        else:
            allType[cCount] =1
    else: pass


# In[16]:


# calculate % methylation at each position and prepare output
cpgNum = ""
cpgVal = ""
for i in range(len(Cpcnt)):
    cpgNum += str(i + 1) + "\t"
    cpgVal += str(round(Cpcnt[i]/(Cpcnt[i] + tpcnt[i]), 2)) + "\t"


# In[17]:


#outfile = open("Methylation_statistics_nhmmer_q"+str(phred_cutoff)+".txt", 'w')
outfile = open("Methylation_statistics_blastn_q30_dash1.txt", 'w')
print("Number of alleles: ", len(alignDict2), file=outfile)
print("", file=outfile)
print("Percent methylation of each CpG", file=outfile)
print(cpgNum, file=outfile)
print(cpgVal, file=outfile)
print("", file=outfile)
print("Number of alleles with methylation percent", file=outfile)
print("<10\t<20\t<30\t<40\t<50\t<60\t<70\t<80\t<90\t90+", file=outfile)
methpcent = ""


# In[18]:


for value in allPcnt:
    methpcent += str(value) + "\t"
print(methpcent, file=outfile)
print("", file=outfile)
print("Unique alleles found", file=outfile)
print("Frequency\tNumt\tNumdash\tNumC\tSequence", file=outfile)
sortAllele = sorted(allType.items(), key=lambda x:x[1], reverse=True)
for allele in sortAllele:
    outLine = str(allele[1]) + "\t" + str(allele[0].count("t")) + "\t" + str(allele[0].count("-")) + "\t" + str(allele[0].count("C"))
    for char in allele[0]:
        outLine += "\t" + char
    print(outLine, file=outfile)
outfile.close()


# In[ ]:




