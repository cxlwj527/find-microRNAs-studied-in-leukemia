# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 09:32:31 2020
This script is aimed to find which miRNAs are reported in PubMed related to\
one designated disease or gene. 
@author: Weijun Liu
"""
from Bio import Entrez
search_term="leukemia and microRNA"
Entrez.email = "your email"
handle = Entrez.egquery(term=search_term)
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
    if row['DbName']=="pubmed":
        paper_numbers=row["Count"]

handle = Entrez.esearch(db="pubmed", \
                        term=search_term, \
                        retmax=paper_numbers)
record = Entrez.read(handle)
handle.close()
idlist=record["IdList"]

from Bio import Medline
handle=Entrez.efetch(db="pubmed", id=idlist, rettype='medline', retmode='text')
records=Medline.parse(handle)
records=list(records)


import re
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


miRNA_counts={}
for record in records:
    if record["PMID"] in paper_find:
        try:
            find_miRNAs=[]
            for miRNA_regEx in miRNA_regExs:
                miRNAs=re.findall(miRNA_regEx, record["AB"])
                find_miRNAs = find_miRNAs + miRNAs
                unique_miRNAs=list(set(find_miRNAs))
                    
            miRNA_list=[]
            for miRNA in unique_miRNAs:
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
                    miRNA=re.sub('-3p', '', miRNA)
                    miRNA=re.sub('-5p', '', miRNA)
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
                    miRNA=re.sub('-3p', '', miRNA)
                    miRNA=re.sub('-5p', '', miRNA)
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
                    miRNA=re.sub('-3p', '', miRNA)
                    miRNA=re.sub('-5p', '', miRNA)
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
                    miRNA=re.sub('-3p', '', miRNA)
                    miRNA=re.sub('-5p', '', miRNA)
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
                    miRNA_list.append(miRNA)
            for miRNA in list(set(miRNA_list)):
                miRNA_counts[miRNA]=miRNA_counts.get(miRNA,0)+1
        except KeyError:
            continue
                
import operator    
sorted_miRNA_counts=dict(sorted(miRNA_counts.items(), \
                                key=operator.itemgetter(1),reverse=True))
from itertools import islice
miRNAs = list(sorted_miRNA_counts.keys())
if len(miRNAs) >20:
    miRNAs_20=dict(islice(sorted_miRNA_counts.items(), 20))
else:
    miRNAs_20=sorted_miRNA_counts
    
import matplotlib.pyplot as plt
import numpy as np
fig, ax = plt.subplots()
miRNAs = list(miRNAs_20.keys())
counts= miRNAs_20.values()
x = np.arange(len(miRNAs))
ax.bar(x, counts, color='#0343df')
ax.set_ylabel('Papers')
ax.set_title('Most studied 20 miRNAs in '+ search_term.split(' ')[0])
ax.set_xticks(x)    
ax.set_xticklabels(miRNAs, rotation=75)
plt.show()
