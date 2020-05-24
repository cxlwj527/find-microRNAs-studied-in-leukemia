# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 09:32:31 2020
This script is aimed to find which miRNAs are reported in PubMed related to\
one designated disease or gene. 
@author: Weijun
"""
import textwrap
import re
import os
from Bio import Entrez

search_term="leukemia and microRNA"

Entrez.email = "cxlwj527@gmail.com"
handle = Entrez.egquery(term=search_term)
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
    if row['DbName']=="pubmed":
        print(row["Count"])
        paper_numbers=row["Count"]

handle = Entrez.esearch(db="pubmed", \
                        term=search_term, \
                        retmax=paper_numbers)
record = Entrez.read(handle)
handle.close()
idlist=record["IdList"]
#print(idlist)

from Bio import Medline
handle=Entrez.efetch(db="pubmed", id=idlist, rettype='medline', retmode='text')
records=Medline.parse(handle)
records=list(records)
#print(len(records))

miRNA_regEx1=re.compile(r'mir-?\s?\d+[a-z]?-[35]p', re.IGNORECASE)
miRNA_regEx2=re.compile(r'mir-?\s?\d+[a-z]?', re.IGNORECASE)
miRNA_regEx3=re.compile(r'mirna-?\s?\d+[a-z]?-[35]p', re.IGNORECASE)
miRNA_regEx4=re.compile(r'mirna-?\s?\d+[a-z]?', re.IGNORECASE)
miRNA_regEx5=re.compile(r'microrna-?\s?\d+[a-z]?-[35]p', re.IGNORECASE)
miRNA_regEx6=re.compile(r'microrna-?\s?\d+[a-z]?', re.IGNORECASE)
miRNA_regEx7=re.compile(r'let-?\s?\d+[a-z]?-[35]p', re.IGNORECASE)
miRNA_regEx8=re.compile(r'let-?\s?\d+[a-z]?', re.IGNORECASE)
miRNA_regExs=[miRNA_regEx1,miRNA_regEx2,miRNA_regEx3,miRNA_regEx4,\
              miRNA_regEx5,miRNA_regEx6,miRNA_regEx7,miRNA_regEx8]
 

paper_find = []
for record in records:
    try:
        for miRNA_regEx in miRNA_regExs:
            match=re.search(miRNA_regEx, record['AB'])
            if match != None:
                paper_find.append(record['PMID'])
            else:
                continue
    except KeyError:
        continue

paper_find=list(set(paper_find))
print(len(paper_find))



miRNA_counts={}
with open("processing.txt","a") as handle:
    for record in records:
        if record["PMID"] in paper_find:
            try:
                handle.write(textwrap.fill(record["SO"],width=100)+'\n')
                handle.write(textwrap.fill(record["TI"],width=100)+'\n')
                handle.write(textwrap.fill(record["AB"],width=100)+'\n')
                handle.write("PMID"+record["PMID"]+'\n')
                find_miRNAs=[]
                for miRNA_regEx in miRNA_regExs:
                    miRNAs=re.findall(miRNA_regEx, record["AB"])
                    find_miRNAs = find_miRNAs + miRNAs
                    unique_miRNAs=list(set(find_miRNAs))
                    #handle.write("miRNAs in paper PMID"+record["PMID"]+'\n')
                    print("miRNA in paper PMID", record["PMID"])
                    print(unique_miRNAs)
                    
                miRNA_list=[]
                for miRNA in unique_miRNAs:
                    handle.write('miRNA found: '+miRNA+'\n')
                    print('miRNA found: ',miRNA)
                    if re.match(miRNA_regEx1, miRNA):
                        miRNA=re.sub('miR-', 'miR', miRNA)
                        miRNA=re.sub('miR ', 'miR', miRNA)
                        miRNA=re.sub('mir-', 'miR', miRNA)
                        miRNA=re.sub('mir ', 'miR', miRNA)
                        miRNA=re.sub('mir', 'miR', miRNA)
                        miRNA=re.sub('MIR-', 'miR', miRNA)
                        miRNA=re.sub('MIR ', 'miR', miRNA)
                        miRNA=re.sub('MIR', 'miR', miRNA)
                        miRNA=re.sub('Mir-', 'miR', miRNA)
                        miRNA=re.sub('Mir ', 'miR', miRNA)
                        miRNA=re.sub('Mir', 'miR', miRNA)
                        miRNA=re.sub('MiR-', 'miR', miRNA)
                        miRNA=re.sub('MiR ', 'miR', miRNA)
                        miRNA=re.sub('MiR', 'miR', miRNA)
                        #miRNA=re.sub('3p', '', miRNA)
                        #miRNA=re.sub('5p', '', miRNA)
                        miRNA=re.sub('-3p', '', miRNA)
                        miRNA=re.sub('-5p', '', miRNA)
                        handle.write('Formatted miRNA: '+ miRNA+'\n')
                        print('Formatted miRNA: ', miRNA)
                        miRNA_list.append(miRNA)
                          
                    elif re.match(miRNA_regEx2, miRNA):
                        miRNA=re.sub('miR-', 'miR', miRNA)
                        miRNA=re.sub('miR ', 'miR', miRNA)
                        miRNA=re.sub('mir-', 'miR', miRNA)
                        miRNA=re.sub('mir ', 'miR', miRNA)
                        miRNA=re.sub('mir', 'miR', miRNA)
                        miRNA=re.sub('MIR-', 'miR', miRNA)
                        miRNA=re.sub('MIR ', 'miR', miRNA)
                        miRNA=re.sub('MIR', 'miR', miRNA)
                        miRNA=re.sub('Mir-', 'miR', miRNA)
                        miRNA=re.sub('Mir ', 'miR', miRNA)
                        miRNA=re.sub('Mir', 'miR', miRNA)
                        miRNA=re.sub('MiR-', 'miR', miRNA)
                        miRNA=re.sub('MiR ', 'miR', miRNA)
                        miRNA=re.sub('MiR', 'miR', miRNA)
                        handle.write('Formatted miRNA: '+ miRNA+'\n')
                        print('Formatted miRNA: ', miRNA)
                        miRNA_list.append(miRNA)
                    elif re.match(miRNA_regEx3, miRNA):
                        miRNA=re.sub('miRNA-', 'miR', miRNA)
                        miRNA=re.sub('miRNA ', 'miR', miRNA)
                        miRNA=re.sub('miRNA', 'miR', miRNA)
                        miRNA=re.sub('MiRNA-', 'miR', miRNA)
                        miRNA=re.sub('MiRNA ', 'miR', miRNA)
                        miRNA=re.sub('MiRNA', 'miR', miRNA)
                        miRNA=re.sub('mirna-', 'miR', miRNA)
                        miRNA=re.sub('mirna ', 'miR', miRNA)
                        miRNA=re.sub('mirna', 'miR', miRNA)
                        miRNA=re.sub('MIRNA-', 'miR', miRNA)
                        miRNA=re.sub('MIRNA ', 'miR', miRNA)
                        miRNA=re.sub('MIRNA', 'miR', miRNA)
                        #miRNA=re.sub('3p', '', miRNA)
                        #miRNA=re.sub('5p', '', miRNA)
                        miRNA=re.sub('-3p', '', miRNA)
                        miRNA=re.sub('-5p', '', miRNA)
                        handle.write('Formatted miRNA: '+ miRNA+'\n')
                        print('Formatted miRNA: ', miRNA)
                        miRNA_list.append(miRNA)
                    elif re.match(miRNA_regEx4, miRNA):
                        miRNA=re.sub('miRNA-', 'miR', miRNA)
                        miRNA=re.sub('miRNA ', 'miR', miRNA)
                        miRNA=re.sub('miRNA', 'miR', miRNA)
                        miRNA=re.sub('MiRNA-', 'miR', miRNA)
                        miRNA=re.sub('MiRNA ', 'miR', miRNA)
                        miRNA=re.sub('MiRNA', 'miR', miRNA)
                        miRNA=re.sub('mirna-', 'miR', miRNA)
                        miRNA=re.sub('mirna ', 'miR', miRNA)
                        miRNA=re.sub('mirna', 'miR', miRNA)
                        miRNA=re.sub('MIRNA-', 'miR', miRNA)
                        miRNA=re.sub('MIRNA ', 'miR', miRNA)
                        miRNA=re.sub('MIRNA', 'miR', miRNA)
                        handle.write('Formatted miRNA: '+ miRNA+'\n')
                        print('Formatted miRNA: ', miRNA)
                        miRNA_list.append(miRNA)
                    elif re.match(miRNA_regEx5, miRNA):
                        miRNA=re.sub('microRNA-', 'miR', miRNA)
                        miRNA=re.sub('microRNA ', 'miR', miRNA)
                        miRNA=re.sub('microRNA', 'miR', miRNA)
                        miRNA=re.sub('microrna-', 'miR', miRNA)
                        miRNA=re.sub('microrna ', 'miR', miRNA)
                        miRNA=re.sub('microrna', 'miR', miRNA)
                        miRNA=re.sub('MICRORNA-', 'miR', miRNA)
                        miRNA=re.sub('MICRORNA', 'miR', miRNA)
                        miRNA=re.sub('MICRORNA ', 'miR', miRNA)
                        miRNA=re.sub('MicroRNA-', 'miR', miRNA)
                        miRNA=re.sub('MicroRNA ', 'miR', miRNA)
                        miRNA=re.sub('MicroRNA', 'miR', miRNA)
                        miRNA=re.sub('Microrna-', 'miR', miRNA)
                        miRNA=re.sub('Microrna ', 'miR', miRNA)
                        miRNA=re.sub('Microrna', 'miR', miRNA)
                        #miRNA=re.sub('3p', '', miRNA)
                        #miRNA=re.sub('5p', '', miRNA)
                        miRNA=re.sub('-3p', '', miRNA)
                        miRNA=re.sub('-5p', '', miRNA)
                        handle.write('Formatted miRNA: '+ miRNA+'\n')
                        print('Formatted miRNA: ', miRNA)
                        miRNA_list.append(miRNA)
                    elif re.match(miRNA_regEx6, miRNA):
                        miRNA=re.sub('microRNA-', 'miR', miRNA)
                        miRNA=re.sub('microRNA ', 'miR', miRNA)
                        miRNA=re.sub('microRNA', 'miR', miRNA)
                        miRNA=re.sub('microrna-', 'miR', miRNA)
                        miRNA=re.sub('microrna ', 'miR', miRNA)
                        miRNA=re.sub('microrna', 'miR', miRNA)
                        miRNA=re.sub('MICRORNA-', 'miR', miRNA)
                        miRNA=re.sub('MICRORNA', 'miR', miRNA)
                        miRNA=re.sub('MICRORNA ', 'miR', miRNA)
                        miRNA=re.sub('MicroRNA-', 'miR', miRNA)
                        miRNA=re.sub('MicroRNA ', 'miR', miRNA)
                        miRNA=re.sub('MicroRNA', 'miR', miRNA)
                        miRNA=re.sub('Microrna-', 'miR', miRNA)
                        miRNA=re.sub('Microrna ', 'miR', miRNA)
                        miRNA=re.sub('Microrna', 'miR', miRNA)
                        handle.write('Formatted miRNA: '+ miRNA+'\n')
                        print('Formatted miRNA: ', miRNA)
                        miRNA_list.append(miRNA)
                    elif re.match(miRNA_regEx7, miRNA):
                        miRNA=re.sub('Let-', 'Let', miRNA)
                        miRNA=re.sub('Let ', 'Let', miRNA)
                        miRNA=re.sub('Let', 'Let', miRNA)
                        miRNA=re.sub('let-', 'Let', miRNA)
                        miRNA=re.sub('let ', 'Let', miRNA)
                        miRNA=re.sub('let', 'Let', miRNA)
                        miRNA=re.sub('LET-', 'Let', miRNA)
                        miRNA=re.sub('LET ', 'Let', miRNA)
                        miRNA=re.sub('LET', 'Let', miRNA)
                        #miRNA=re.sub('3p', '', miRNA)
                        #miRNA=re.sub('5p', '', miRNA)
                        miRNA=re.sub('-3p', '', miRNA)
                        miRNA=re.sub('-5p', '', miRNA)
                        handle.write('Formatted miRNA: '+ miRNA+'\n')
                        print('Formatted miRNA: ', miRNA)
                        miRNA_list.append(miRNA)
                    elif re.match(miRNA_regEx8, miRNA):
                        miRNA=re.sub('Let-', 'Let', miRNA)
                        miRNA=re.sub('Let ', 'Let', miRNA)
                        miRNA=re.sub('Let', 'Let', miRNA)
                        miRNA=re.sub('let-', 'Let', miRNA)
                        miRNA=re.sub('let ', 'Let', miRNA)
                        miRNA=re.sub('let', 'Let', miRNA)
                        miRNA=re.sub('LET-', 'Let', miRNA)
                        miRNA=re.sub('LET ', 'Let', miRNA)
                        miRNA=re.sub('LET', 'Let', miRNA)
                        handle.write('Formatted miRNA: '+ miRNA+'\n')
                        print('Formatted miRNA: ', miRNA)
                        miRNA_list.append(miRNA)
                        
                   
                print(len(miRNA_list))
                print(len(set(miRNA_list)))
                for miRNA in list(set(miRNA_list)):
                    miRNA_counts[miRNA]=miRNA_counts.get(miRNA,0)+1
                handle.write('\n')
            except KeyError:
                continue
                
outfilename=format('_'.join(search_term.split(' ')))+'.txt'
if os.path.exists(outfilename):
    os.remove(outfilename)
    os.rename('Processing.txt',outfilename)
else:
    os.rename('Processing.txt',outfilename)

#for miRNA, count in sorted(miRNA_counts.items(), \
                           #key=lambda kv: kv[1], reverse=True):
    #print("%s: %s" % (miRNA, count))
    
import operator    
sorted_miRNA_counts=dict(sorted(miRNA_counts.items(), \
                                key=operator.itemgetter(1),reverse=True))

from itertools import islice
#def take(n, iterable):
    #return list(islice(iterable, n))

miRNAs = list(sorted_miRNA_counts.keys())
if len(miRNAs) >20:
    #miRNAs_20=take(20, sorted_miRNA_counts.items())
    miRNAs_20=dict(islice(sorted_miRNA_counts.items(), 20))
else:
    miRNAs_20=sorted_miRNA_counts
    
import csv
with open ('processing.csv', 'w') as fh:
    fh.write('miRNA,counts\n\n')
    for miRNA in sorted_miRNA_counts.keys():
        fh.write('%s, %s\n'%(miRNA,miRNA_counts[miRNA]))

outfilename=format('_'.join(search_term.split(' '))) +'.csv'
if os.path.exists(outfilename):
    os.remove(outfilename)
    os.rename('Processing.csv',outfilename)
else:
    os.rename('Processing.csv',outfilename)


import matplotlib.pyplot as plt
import numpy as np
#plt.bar(range(len(miRNA_counts)), miRNA_counts.values(), align='center')
#plt.xticks(range(len(miRNA_counts)), list(miRNA_counts.keys()))
#plt.show()

fig, ax = plt.subplots()

# A little data preparation
miRNAs = list(miRNAs_20.keys())
counts= miRNAs_20.values()
x = np.arange(len(miRNAs))

# Plot each bar plot. Note: manually calculating the 'dodges' of the bars
ax.bar(x, counts, color='#0343df')


# Customise some display properties
ax.set_ylabel('Papers')
ax.set_title('Most 20 studied miRNAs in '+ search_term.split(' ')[0])
ax.set_xticks(x)    # This ensures we have one tick per miRNA, otherwise we get fewer
ax.set_xticklabels(miRNAs, rotation=75)
#ax.legend()

# Ask Matplotlib to show the plot
plt.show()