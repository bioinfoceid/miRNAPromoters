"""
On the computational prediction of miRNA promoters
Charalampos Michail1, Aigli Korfiati1,2, Konstantinos Theofilatos2, Spiros Likothanassis1, Seferina Mavroudi3
1 Department of Computer Engineering and Informatics, University of Patras, Patra, Greece
{cmichail,korfiati,likothan}@ceid.upatras.gr
2InSyBio Ltd, London, UK
k.theofilatos@insybio.com
3 Department of Social Work, Technological Institute of Western Greece, Patra, Greece
mavroudi@teiwest.gr
"""


#!/usr/bin/python
import csv
from __future__ import print_function
import random 
import re
import csv
from collections import defaultdict

vim = open("tixaia_wc.txt", "r")
v = open("data.txt", "r")
f5 = open("characteristics.txt","w")
f2 = open("cry.txt","w")
f1 = open("chromatin_results.txt","w")

data = list(csv.reader(open('wgEncodeBroadHmmGm12878HMM.csv')))
data2 = list(csv.reader(open('wgEncodeBroadHmmH1hescHMM.csv')))
data3 = list(csv.reader(open('wgEncodeBroadHmmHmecHMM.csv')))
data4 = list(csv.reader(open('wgEncodeBroadHmmHsmmHMM_1.csv')))
data5 = list(csv.reader(open('wgEncodeBroadHmmHuvecHMM.csv')))
data6 = list(csv.reader(open('wgEncodeBroadHmmNhekHMM.csv')))
data7 = list(csv.reader(open('wgEncodeBroadHmmNhlfHMM_1.csv')))

der = list(csv.reader(open(chromatin_data.csv')))

sum1=0
sum2=0
sum3=0
sum4=0
sum5=0
sum6=0
sum7=0
sum8=0
sum9=0
sum10=0
sum11=0
sum12=0
sum13=0
sum14=0
sum15=0

print("1_Active_Promoter(Gm12878) \t 2_Weak_Promoter(Gm12878) \t 3_Poised_Promoter(Gm12878) \t 4_Strong_Enhancer(Gm12878) \t 5_Strong_Enhancer(Gm12878) \t 6_Weak_Enhancer(Gm12878) \t 7(Gm12878) \t 8_Insulator(Gm12878) \t 9_Txn_Transition(Gm12878) \t 10_Txn_Elongation(Gm12878) \t 11_Weak_Txn(Gm12878) \t 12_Repressed(Gm12878) \t 13_Heterochrom/lo(Gm12878) \t 14_Repetitive/CNV(Gm12878) \t 15_Repetitive/CNV(Gm12878) \t 1_Active_Promoter(H1hesc) \t 2_Weak_Promoter(H1hesc) \t 3_Poised_Promoter(H1hesc) \t 4_Strong_Enhancer(H1hesc) \t 5_Strong_Enhancer(H1hesc) \t 6_Weak_Enhancer(H1hesc) \t 7_Weak_Enhancer(H1hesc) \t 8_Insulator(H1hesc) \t 9_Txn_Transition(H1hesc) \t 10_Txn_Elongation(H1hesc) \t 11_Weak_Txn(H1hesc) \t 12_Repressed(H1hesc) \t 13_Heterochrom/lo(H1hesc) \t 14_Repetitive/CNV(H1hesc) \t 15_Repetitive/CNV(H1hesc) \t 1_Active_Promoter(Hmec) \t 2_Weak_Promoter(Hmec) \t 3_Poised_Promoter(Hmec) \t 4_Strong_Enhancer(Hmec) \t 5_Strong_Enhancer(Hmec) \t 6_Weak_Enhancer(Hmec) \t 7_Weak_Enhancer(Hmec) \t 8_Insulator(Hmec) \t 9_Txn_Transition(Hmec) \t 10_Txn_Elongation(Hmec) \t 11_Weak_Txn(Hmec) \t 12_Repressed(Hmec) \t 13_Heterochrom/lo(Hmec) \t 14_Repetitive/CNV(Hmec) \t 15_Repetitive/CNV(Hmec) \t 1_Active_Promoter(Hsmm) \t 2_Weak_Promoter(Hsmm) \t 3_Poised_Promoter(Hsmm) \t 4_Strong_Enhancer(Hsmm) \t 5_Strong_Enhancer(Hsmm) \t 6_Weak_Enhancer(Hsmm) \t 7_Weak_Enhancer(Hsmm) \t 8_Insulator(Hsmm) \t 9_Txn_Transition(Hsmm) \t 10_Txn_Elongation(Hsmm) \t 11_Weak_Txn(Hsmm) \t 12_Repressed(Hsmm) \t 13_Heterochrom/lo(Hsmm) \t 14_Repetitive/CNV(Hsmm) \t 15_Repetitive/CNV(Hsmm) \t 1_Active_Promoter(Huvec) \t 2_Weak_Promoter(Huvec) \t 3_Poised_Promoter(Huvec) \t 4_Strong_Enhancer(Huvec) \t 5_Strong_Enhancer(Huvec) \t 6_Weak_Enhancer(Huvec) \t 7_Weak_Enhancer(Huvec) \t 8_Insulator(Huvec) \t 9_Txn_Transition(Huvec) \t 10_Txn_Elongation(Huvec) \t 11_Weak_Txn(Huvec) \t 12_Repressed(Huvec) \t 13_Heterochrom/lo(Huvec) \t 14_Repetitive/CNV(Huvec) \t 15_Repetitive/CNV(Huvec) \t 1_Active_Promoter(Nhek) \t 2_Weak_Promoter(Nhek) \t 3_Poised_Promoter(Nhek) \t 4_Strong_Enhancer(Nhek) \t 5_Strong_Enhancer(Nhek) \t 6_Weak_Enhancer(Nhek) \t 7_Weak_Enhancer(Nhek) \t 8_Insulator(Nhek) \t 9_Txn_Transition(Nhek) \t 10_Txn_Elongation(Nhek) \t 11_Weak_Txn(Nhek) \t 12_Repressed(Nhek) \t 13_Heterochrom/lo(Nhek) \t 14_Repetitive/CNV(Nhek) \t 15_Repetitive/CNV(Nhek) \t 1_Active_Promoter(Nhlf) \t 2_Weak_Promoter(Nhlf) \t 3_Poised_Promoter(Nhlf) \t 4_Strong_Enhancer(Nhlf) \t 5_Strong_Enhancer(Nhlf) \t 6_Weak_Enhancer(Nhlf) \t 7_Weak_Enhancer(Nhlf) \t 8_Insulator(Nhlf) \t 9_Txn_Transition(Nhlf) \t 10_Txn_Elongation(Nhlf) \t 11_Weak_Txn(Nhlf) \t 12_Repressed(Nhlf) \t 13_Heterochrom/lo(Nhlf) \t 14_Repetitive/CNV(Nhlf) \t 15_Repetitive/CNV(Nhlf)", file = f1)

x=len(data) 
x2=len(data2)
x3=len(data3)
x4=len(data4)
x5=len(data5)
x6=len(data6)
x7=len(data7)
a_count = 0
g_count = 0
c_count = 0
t_count = 0
atskew=0.0
gcskew=0.0
pith=0.0
at_sum=0
cg_sum=0
cg_count=0
elaxisto=0
megisto=0
dna=''
upstream=1
diafora_rep=0
sum_rep=0
sum_palindromes = 0
sum_wc= 0.0
all_counts = []
all_pith = []
all_nucleotides = [] 
organisms = []
myseq=''

y=0
duk=0

def sliding_window(sequence, winSize, step):
        numOfChunks = ((len(sequence)-winSize)/step)+1
        numOfChunks=int(numOfChunks)
        for i in range(0,numOfChunks*step,step):
                yield sequence[i:i+winSize]

				
def reverse_complement(s):
    complements = {'a':'t', 't':'a', 'g':'c', 'c':'g' }
    return ''.join([complements[c] for c in reversed(s)])


def reverse_palindromes(s):
    results = []

    l = len(s)

    for i in range(l):
        for j in range(6, 501):

            if i + j > l:
                continue

            s1 = s[i:i+j]
            s2 = reverse_complement(s1)

            if s1 == s2:
                results.append((i + 1, j))

    return results	

#-----------------------------word commonality--------------------------------------------

for line in vim:
     if '>' in line :
         myseq= myseq + 'X'
     if not line.startswith(">"):
         myseq=myseq+line
     if  '>end' in line :
         for base1 in ['a', 't', 'g', 'c']:
                 for base2 in ['a', 't', 'g', 'c']:
                         for base3 in ['a', 't', 'g', 'c']:
                             trinucleotide = base1 + base2 + base3
                             all_nucleotides.append(trinucleotide)
                             counta = myseq.count(trinucleotide)
                             #print("count is " + str(count) + " for " + trinucleotide)
                             all_counts.append(counta)
         myseq=''
        # print (all_counts)
         elaxisto=min(all_counts)
         megisto=max(all_counts)
         print(all_counts,file=f2)
         print(all_nucleotides,file=f2)
         #print ("min", x)
         #print ("max", y)

j=0
i=0
for j in xrange (64):
     if megisto != 0 :
         print(all_counts[j])
         pith=float(all_counts[j] - elaxisto)/(megisto-elaxisto)
         all_pith.append(pith)
print(all_pith,file=f2)
i=0
j=0				

#------------------------end word commonality---------------------------------------------
				
#------------------------ first line------------------------------------------------------				
print ("\t obs/exp(AA) \t obs/exp(AAA) \t obs/exp(AAT) \t obs/exp(AAG) \t obs/exp(AAC) \t obs/exp(AT) \t obs/exp(ATA) \t obs/exp(ATT) \t obs/exp(ATG) \t obs/exp(ATC) \t obs/exp(AG) \t obs/exp(AGA) \t obs/exp(AGT) \t obs/exp(AGG) \t obs/exp(AGC) \t obs/exp(AC) \t obs/exp(ACA) \t obs/exp(ACT) \t obs/exp(ACG) \t obs/exp(ACC) \t obs/exp(TA) \t obs/exp(TAA) \t obs/exp(TAT) \t obs/exp(TAG) \t obs/exp(TAC) \t obs/exp(TT) \t obs/exp(TTA) \t obs/exp(TTT) \t obs/exp(TTG) \t obs/exp(TTC) \t obs/exp(TG) \t obs/exp(TGA) \t obs/exp(TGT) \t obs/exp(TGG) \t obs/exp(TGC) \t obs/exp(TC) \t obs/exp(TCA) \t obs/exp(TCT) \t obs/exp(TCG) \t obs/exp(TCC) \t obs/exp(GA) \t obs/exp(GAA) \t obs/exp(GAT) \t obs/exp(GAG) \t obs/exp(GAC) \t obs/exp(GT) \t obs/exp(GTA) \t obs/exp(GTT) \t obs/exp(GTG) \t obs/exp(GTC) \t obs/exp(GG) \t obs/exp(GGA) \t obs/exp(GGT) \t obs/exp(GGG) \t obs/exp(GGC) \t obs/exp(GC) \t obs/exp(GCA) \t obs/exp(GCT) \t obs/exp(GCG) \t obs/exp(GCC) \t obs/exp(CA) \t obs/exp(CAA) \t obs/exp(CAT) \t obs/exp(CAG) \t obs/exp(CAC) \t obs/exp(CT) \t obs/exp(CTA) \t obs/exp(CTT) \t obs/exp(CTG) \t obs/exp(CTC) \t obs/exp(CG) \t obs/exp(CGA) \t obs/exp(CGT) \t obs/exp(CGG) \t obs/exp(CGC) \t obs/exp(CC) \t obs/exp(CCA) \t obs/exp(CCT) \t obs/exp(CCG) \t obs/exp(CCC) \t upstream(AA) \t upstream(AAA) \t upstream(AAT) \t upstream(AAG) \t upstream(AAC) \t upstream(AT) \t upstream(ATA) \t upstream(ATT) \t upstream(ATG) \t upstream(ATC) \t upstream(AG) \t upstream(AGA) \t upstream(AGT) \t upstream(AGG) \t upstream(AGC) \t upstream(AC) \t upstream(ACA) \t upstream(ACT) \t upstream(ACG) \t upstream(ACC) \t upstream(TA) \t upstream(TAA) \t upstream(TAT) \t upstream(TAG) \t upstream(TAC) \t upstream(TT) \t upstream(TTA) \t upstream(TTT) \t upstream(TTG) \t upstream(TTC) \t upstream(TG) \t upstream(TGA) \t upstream(TGT) \t upstream(TGG) \t upstream(TGC) \t upstream(TC) \t upstream(TCA) \t upstream(TCT) \t upstream(TCG) \t upstream(TCC) \t upstream(GA) \t upstream(GAA) \t upstream(GAT) \t upstream(GAG) \t upstream(GAC) \t upstream(GT) \t upstream(GTA) \t upstream(GTT) \t upstream(GTG) \t upstream(GTC) \t upstream(GG) \t upstream(GGA) \t upstream(GGT) \t upstream(GGG) \t upstream(GGC) \t upstream(GC) \t upstream(GCA) \t upstream(GCT) \t upstream(GCG) \t upstream(GCC) \t upstream(CA) \t upstream(CAA) \t upstream(CAT) \t upstream(CAG) \t upstream(CAC) \t upstream(CT) \t upstream(CTA) \t upstream(CTT) \t upstream(CTG) \t upstream(CTC) \t upstream(CG) \t upstream(CGA) \t upstream(CGT) \t upstream(CGG) \t upstream(CGC) \t upstream(CC) \t upstream(CCA) \t upstream(CCT) \t upstream(CCG) \t upstream(CCC) \t upstream(G) \t upstream(A) \t upstream(T) \t upstream(C) \t downstream(AA) \t downstream(AAA) \t downstream(AAT) \t downstream(AAG) \t downstream(AAC) \t downstream(AT) \t downstream(ATA) \t downstream(ATT) \t downstream(ATG) \t downstream(ATC) \t downstream(AG) \t downstream(AGA) \t downstream(AGT) \t downstream(AGG) \t downstream(AGC) \t downstream(AC) \t downstream(ACA) \t downstream(ACT) \t downstream(ACG) \t downstream(ACC) \t downstream(TA) \t downstream(TAA) \t downstream(TAT) \t downstream(TAG) \t downstream(TAC) \t downstream(TT) \t downstream(TTA) \t downstream(TTT) \t downstream(TTG) \t downstream(TTC) \t downstream(TG) \t downstream(TGA) \t downstream(TGT) \t downstream(TGG) \t downstream(TGC) \t downstream(TC) \t downstream(TCA) \t downstream(TCT) \t downstream(TCG) \t downstream(TCC) \t downstream(GA) \t downstream(GAA) \t downstream(GAT) \t downstream(GAG) \t downstream(GAC) \t downstream(GT) \t downstream(GTA) \t downstream(GTT) \t downstream(GTG) \t downstream(GTC) \t downstream(GG) \t downstream(GGA) \t downstream(GGT) \t downstream(GGG) \t downstream(GGC) \t downstream(GC) \t downstream(GCA) \t downstream(GCT) \t downstream(GCG) \t downstream(GCC) \t downstream(CA) \t downstream(CAA) \t downstream(CAT) \t downstream(CAG) \t downstream(CAC) \t downstream(CT) \t downstream(CTA) \t downstream(CTT) \t downstream(CTG) \t downstream(CTC) \t downstream(CG) \t downstream(CGA) \t downstream(CGT) \t downstream(CGG) \t downstream(CGC) \t downstream(CC) \t downstream(CCA) \t downstream(CCT) \t downstream(CCG) \t downstream(CCC) \t downstream(G) \t downstream(A) \t downstream(T) \t downstream(C) \t GC-skew(1) \t AT-skew(1) \t repetitives(1) \t palindromes(1) \t CpG island(1) \t CpGs(1) \t WC(1) \t GC-skew(2) \t AT-skew(2) \t repetitives(2) \t palindromes(2) \t CpG island(2) \t CpGs(2) \t WC(2) \t GC-skew(3) \t AT-skew(3) \t repetitives(3) \t palindromes(3) \t CpG island(3) \t CpGs(3) \t WC(3) \t GC-skew(4) \t AT-skew(4) \t repetitives(4) \t palindromes(4) \t CpG island(4) \t CpGs(4) \t WC(4) \t GC-skew(5) \t AT-skew(5) \t repetitives(5) \t palindromes(5) \t CpG island(5) \t CpGs(5) \t WC(5) \t GC-skew(6) \t AT-skew(6) \t repetitives(6) \t palindromes(6) \t CpG island(6) \t CpGs(6) \t WC(6) \t GC-skew(7) \t AT-skew(7) \t repetitives(7) \t palindromes(7) \t CpG island(7) \t CpGs(7) \t WC(7) \t GC-skew(8) \t AT-skew(8) \t repetitives(8) \t palindromes(8) \t CpG island(8) \t CpGs(8) \t WC(8) \t GC-skew(9) \t AT-skew(9) \t repetitives(9) \t palindromes(9) \t CpG island(9) \t CpGs(9) \t WC(9) \t GC-skew(10) \t AT-skew(10) \t repetitives(10) \t palindromes(10) \t CpG island(10) \t CpGs(10) \t WC(10) \t GC-skew(11) \t AT-skew(11) \t repetitives(11) \t palindromes(11) \t CpG island(11) \t CpGs(11) \t WC(11) \t GC-skew(12) \t AT-skew(12) \t repetitives(12) \t palindromes(12) \t CpG island(12) \t CpGs(12) \t WC(12) \t GC-skew(13) \t AT-skew(13) \t repetitives(13) \t palindromes(13) \t CpG island(13) \t CpGs(13) \t WC(13) \t GC-skew(14) \t AT-skew(14) \t repetitives(14) \t palindromes(14) \t CpG island(14) \t CpGs(14) \t WC(14) \t GC-skew(15) \t AT-skew(15) \t repetitives(15) \t palindromes(15) \t CpG island(15) \t CpGs(15) \t WC(15) \t GC-skew(16) \t AT-skew(16) \t repetitives(16) \t palindromes(16) \t CpG island(16) \t CpGs(16) \t WC(16) \t GC-skew(17) \t AT-skew(17) \t repetitives(17) \t palindromes(17) \t CpG island(17) \t CpGs(17) \t WC(17) \t GC-skew(18) \t AT-skew(18) \t repetitives(18) \t palindromes(18) \t CpG island(18) \t CpGs(18) \t WC(18) \t GC-skew(19) \t AT-skew(19) \t repetitives(19) \t palindromes(19) \t CpG island(19) \t CpGs(19) \t WC(19) \t GC-skew(20) \t AT-skew(20) \t repetitives(20) \t palindromes(20) \t CpG island(20) \t CpGs(20) \t WC(20) \t GC-skew(21) \t AT-skew(21) \t repetitives(21) \t palindromes(21) \t CpG island(21) \t CpGs(21) \t WC(21) \t GC-skew(22) \t AT-skew(22) \t repetitives(22) \t palindromes(22) \t CpG island(22) \t CpGs(22) \t WC(22) \t GC-skew(23) \t AT-skew(23) \t repetitives(23) \t palindromes(23) \t CpG island(23) \t CpGs(23) \t WC(23) \t GC-skew(24) \t AT-skew(24) \t repetitives(24) \t palindromes(24) \t CpG island(24) \t CpGs(24) \t WC(24) \t GC-skew(25) \t AT-skew(25) \t repetitives(25) \t palindromes(25) \t CpG island(25) \t CpGs(25) \t WC(25) \t GC-skew(26) \t AT-skew(26) \t repetitives(26) \t palindromes(26) \t CpG island(26) \t CpGs(26) \t WC(26) \t GC-skew(27) \t AT-skew(27) \t repetitives(27) \t palindromes(27) \t CpG island(27) \t CpGs(27) \t WC(27) \t GC-skew(28) \t AT-skew(28) \t repetitives(28) \t palindromes(28) \t CpG island(28) \t CpGs(28) \t WC(28) \t GC-skew(29) \t AT-skew(29) \t repetitives(29) \t palindromes(29) \t CpG island(29) \t CpGs(29) \t WC(29) \t GC-skew(30) \t AT-skew(30) \t repetitives(30) \t palindromes(30) \t CpG island(30) \t CpGs(30) \t WC(30) \t GC-skew(31) \t AT-skew(31) \t repetitives(31) \t palindromes(31) \t CpG island(31) \t CpGs(31) \t WC(31) \t GC-skew(32) \t AT-skew(32) \t repetitives(32) \t palindromes(32) \t CpG island(32) \t CpGs(32) \t WC(32) \t GC-skew(33) \t AT-skew(33) \t repetitives(33) \t palindromes(33) \t CpG island(33) \t CpGs(33) \t WC(33) \t GC-skew(34) \t AT-skew(34) \t repetitives(34) \t palindromes(34) \t CpG island(34) \t CpGs(34) \t WC(34) \t GC-skew(35) \t AT-skew(35) \t repetitives(35) \t palindromes(35) \t CpG island(35) \t CpGs(35) \t WC(35) \t GC-skew(36) \t AT-skew(36) \t repetitives(36) \t palindromes(36) \t CpG island(36) \t CpGs(36) \t WC(36) \t GC-skew(37) \t AT-skew(37) \t repetitives(37) \t palindromes(37) \t CpG island(37) \t CpGs(37) \t WC(37) \t GC-skew(38) \t AT-skew(38) \t repetitives(38) \t palindromes(38) \t CpG island(38) \t CpGs(38) \t WC(38) \t GC-skew(39) \t AT-skew(39) \t repetitives(39) \t palindromes(39) \t CpG island(39) \t CpGs(39) \t WC(39) \t organism",file=f5)				
#print("\n", file=f5)


for line in v:
     if not line.startswith(">"):
         dna=dna+line
         dna=dna.strip(' \t\n\r')
         
         
            					 
     if  '>' in line :
         if '>d' in line:
             organisms.append(line)
             len_dna=len(dna)

         if(duk > 1):
                 print("\t",organisms[y], end="", file=f5)
                 y= y + 1
	#----------------k-mers---------------------------
         for base11 in ['a', 't', 'g', 'c']:
                 monocleotide=base11
                 g1_count=dna.count('g')
                 a1_count=dna.count('a')
                 t1_count=dna.count('t')
                 c1_count=dna.count('c')
                
                 for base21 in ['a', 't', 'g', 'c']:
                         dinucleotide=base11+base21
                         di1_count=dna.count(dinucleotide)
                         #print(dinucleotide, "for: ",di_count)
                         #print(float((di_count)/(r.count(base1)+r.count(base2))))
                         if ((dna.count(base11)+dna.count(base21)) != 0): 
                                 ob=float(di1_count)*len_dna/(dna.count(base11)*dna.count(base21))
                                 #print (dinucleotide, "has obs: ", ob, file=f5)
                                 #print(dinucleotide, "for",di_count)
                                 print ("\t",ob, end="", file=f5)
                         for base31 in ['a', 't', 'g', 'c']:
                                 trinucleotide = base11 + base21 + base31
                                 tri1_count = dna.count(trinucleotide)
                                 if ((dna.count(base11)+dna.count(base21)+dna.count(base31)) != 0): 
                                         tri_ob=float(tri1_count)/(dna.count(base11)+dna.count(base21)+dna.count(base31))
                                         #print (trinucleotide, "has obs: ", tri_ob, file=f5)
                                         print ("\t",tri_ob, end="", file=f5)

			 
			 
         myvect = sliding_window(dna, 1000, 1001)
         for r in myvect:
             #if (upstream==1):
            #     print("\t upstream from TSS")
             #    upstream=0
            # else:
              #   print("\t downstream from TSS")
               #  upstream=1
               
             for base1 in ['a', 't', 'g', 'c']:
                     monocleotide=base1
                     g2_count=r.count('g')
                     a2_count=r.count('a')
                     t2_count=r.count('t')
                     c2_count=r.count('c')

                     for base2 in ['a', 't', 'g', 'c']:
                             dinucleotide=base1+base2
                             di_count=r.count(dinucleotide)
                             print("\t",di_count, end="", file=f5)
                                 #print(float((di_count)/(r.count(base1)+r.count(base2))))
                                 #if ((r.count(base1)+r.count(base2)) != 0): 
                                    # ob=float(di_count)/(r.count(base1)+r.count(base2))
                                     #print (dinucleotide, "for: ", ob)
                            # print(dinucleotide, "is :",di_count)
                             for base3 in ['a', 't', 'g', 'c']:
                                 trinucleotide = base1 + base2 + base3
                                 tri_count = r.count(trinucleotide)
                                 print ("\t", tri_count,end="",file=f5)
             print("\t",g2_count, end="", file=f5)
             print("\t",a2_count, end="", file=f5)
             print("\t",t2_count, end="",file=f5)
             print("\t",c2_count, end="",file=f5)
                 
	 
	 
     #---------------AT and GC skew------------------------    
         myvect = sliding_window(dna, 100, 50)
         for r in myvect:
                # print("\n",r,"\n", file=f5)
                 a_count=r.count('a')
                 g_count=r.count('g')
                 t_count=r.count('t')
                 c_count=r.count('c')
                 at_sum=a_count+t_count
                 cg_sum=c_count+g_count
                 at_sub=a_count-t_count
                 cg_sub=c_count-g_count
                 if (at_sum==0) :
                     atskew=0
                     #print(atskew)
                 else :
                     atskew=float(at_sub)/(at_sum)
                 #print(atskew)
                 if (cg_sum==0) :
                     gcskew = 0
                     #print(gcskew)
                 else :
                     gcskew=float(cg_sub)/(cg_sum)
                 print("\t",gcskew,"\t",atskew,end="", file=f5)


#------------------------repetitive---------------------------------------------

                 runs = re.finditer(r"[aa]{2,100}", r )
                 for match in runs:
                    run_start = match.start()
                    run_end = match.end()
                    #print("AA region from " + str(run_start) + " to " + str(run_end),file=f5)
                    diafora_rep = run_end - run_start
                    sum_rep     = sum_rep + diafora_rep
                 #print("athroisma rep gia AA",sum_rep, file=f5)



                 runs1 = re.finditer(r"[tt]{2,100}", r)
                 for match in runs1:
                    run_start = match.start()
                    run_end = match.end()
                  #  print("TT region from " + str(run_start) + " to " + str(run_end),file=f5)
                    diafora_rep = run_end - run_start
                    sum_rep     = sum_rep + diafora_rep
                 #print("athroisma rep gia TT k panw",sum_rep, file=f5)



	
	
                 runs2 = re.finditer(r"[cc]{2,100}", r)
                 for match in runs2:
                    run_start = match.start()
                    run_end = match.end()
                   # print("CC region from " + str(run_start) + " to " + str(run_end),file=f5)
                    diafora_rep = run_end - run_start
                    sum_rep     = sum_rep + diafora_rep
                # print("athroisma rep gia CC kai panw",sum_rep, file=f5)
	
	
                 runs3 = re.finditer(r"[gg]{2,100}", r)
                 for match in runs3:
                    run_start = match.start()
                    run_end = match.end()
                    #print("GG region from " + str(run_start) + " to " + str(run_end),file=f5)
                    diafora_rep = run_end - run_start
                    sum_rep     = sum_rep + diafora_rep


                 print("\t", sum_rep, end="", file=f15)

#--------------------------------Palindromes-------------------------------------------------------
				 
                 #print ("\n Palindromes are at position and have length:",file=f5)
                 results = reverse_palindromes(r)
                 for v in results:
                     #print ("\n".join([' '.join(map(str, v[1:])) ]),file=f5)
                     weak =' '.join(map(str, v[1:]))
                     xak = int(weak)
                     sum_palindromes = sum_palindromes + xak
                 print ("\t",sum_palindromes, end="", file=f5)
                 
#--------------------------------CpG-------------------------------------------------------------				 
				 
                 cg_count=r.count("cg")
                 x=len(r)
                 sum_count=c_count+g_count
                 if (sum_count!=0):
                    obs_cg=float(cg_count*x)/sum_count
                    print("obs",obs_cg)
                    cg_content=float(sum_count*100)/(a_count+t_count+sum_count)
                    print("content",cg_content)
                 else:
					obs_cg=0
					cg_content=0
                 if (a_count!=0) or (c_count!=0) or (g_count!=0) or (t_count!=0):
					#print ("obs",obs_cg, file=f5)
					#print "cg content", cg_content
                     if (obs_cg>8.5)and (cg_content>21):
                         print ("\t 1",end="", file=f5)
                         print ("\t", cg_count,end="",file=f5)
                         
                     else:
                     	 print ("\t 0",end="", file=f5)
                         print ("\t 0",end="",file=f5)
					     
                         #print ("CpG numbers:",0)

#------------------------------Word Commonality--------------------------------------------------		 
                 for k in xrange(64):
                         if (all_nucleotides[k] in r): 
                             sum_wc= sum_wc +all_pith[k]
                 print ("\t",sum_wc,end="",file=f5)
                 sum_wc = 0.0





                 sum_rep = 0
                 sum_palindromes = 0		 
				 		 
		

         cg_count = 0
         duk = duk + 1
	 
         dna=''
        # print(line,end="",file=f5)

print("\t",organisms[y], file = f5)		
		 
#------------------------------- Chromatin States -------------------------------------
y_13=len(der)
for j in xrange(y_13):
    for i in xrange(x):
        if (data[i][3] == "1_Active_Promoter") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum1 = sum1 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum1 = sum1 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum1 = sum1 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum1 = sum1 + apot
                    else:
                        #print ("0", file=f1)
                        sum1 = sum1 + 0
            else:
                #print ("0", file=f1)
                sum1 = sum1 + 0

    


        elif (data[i][3] == "2_Weak_Promoter") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum2 = sum2 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum2 = sum2 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum2 = sum2 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum2 = sum2 + apot
                    else:
                        #print ("0", file=f1)
                        sum2 = sum2 + 0
            else:
                #print ("0", file=f1)
                sum2 = sum2 + 0


    
        elif (data[i][3] == "3_Poised_Promoter") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum3 = sum3 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum3 = sum3 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum3 = sum3 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum3 = sum3 + apot
                    else:
                        #print ("0", file=f1)
                        sum3 = sum3 + 0
            else:
                #print ("0", file=f1)
                sum3 = sum3 + 0


        elif (data[i][3] == "4_Strong_Enhancer") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum4 = sum4 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum4 = sum4 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum4 = sum4 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum4 = sum4 + apot
                    else:
                        #print ("0", file=f1)
                        sum4 = sum4 + 0
            else:
                #print ("0", file=f1)
                sum4 = sum4 + 0


        elif (data[i][3] == "5_Strong_Enhancer") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum5 = sum5 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum5 = sum5 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum5 = sum5 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum5 = sum5 + apot
                    else:
                        #print ("0", file=f1)
                        sum5 = sum5 + 0
            else:
                #print ("0", file=f1)
                sum5 = sum5 + 0


        elif (data[i][3] == "6_Weak_Enhancer") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum6 = sum6 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum6 = sum6 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum6 = sum6 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum6 = sum6 + apot
                    else:
                        #print ("0", file=f1)
                        sum6 = sum6 + 0
            else:
                #print ("0", file=f1)
                sum6 = sum6 + 0



        elif (data[i][3] == "7_Weak_Enhancer") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum7 = sum7 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum7 = sum7 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum7 = sum7 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum7 = sum7 + apot
                    else:
                        #print ("0", file=f1)
                        sum7 = sum7 + 0
            else:
                #print ("0", file=f1)
                sum7 = sum7 + 0



        elif (data[i][3] == "8_Insulator") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum8 = sum8 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum8 = sum8 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum8 = sum8 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum8 = sum8 + apot
                    else:
                        #print ("0", file=f1)
                        sum8 = sum8 + 0
            else:
                #print ("0", file=f1)
                sum8 = sum8 + 0





        elif (data[i][3] == "9_Txn_Transition") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum9 = sum9 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum9 = sum9 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum9 = sum9 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum9 = sum9 + apot
                    else:
                        #print ("0", file=f1)
                        sum9 = sum9 + 0
            else:
                #print ("0", file=f1)
                sum9 = sum9 + 0




        elif (data[i][3] == "10_Txn_Elongation") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum10 = sum10 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum10 = sum10 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum10 = sum10 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum10 = sum10 + apot
                    else:
                        #print ("0", file=f1)
                        sum10 = sum10 + 0
            else:
                #print ("0", file=f1)
                sum10 = sum10 + 0


        elif (data[i][3] == "11_Weak_Txn") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum11 = sum11 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum11 = sum11 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum11 = sum11 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum11 = sum11 + apot
                    else:
                        #print ("0", file=f1)
                        sum11 = sum11 + 0
            else:
                #print ("0", file=f1)
                sum11 = sum11 + 0


        elif (data[i][3] == "12_Repressed") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum12 = sum12 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum12 = sum12 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum12 = sum12 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum12 = sum12 + apot
                    else:
                        #print ("0", file=f1)
                        sum12 = sum12 + 0
            else:
                #print ("0", file=f1)
                sum12 = sum12 + 0


        elif (data[i][3] == "13_Heterochrom/lo") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum13 = sum13 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum13 = sum13 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum13 = sum13 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum13 = sum13 + apot
                    else:
                        #print ("0", file=f1)
                        sum13 = sum13 + 0
            else:
                #print ("0", file=f1)
                sum13 = sum13 + 0


        elif (data[i][3] == "14_Repetitive/CNV") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum14 = sum14 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum14 = sum14 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum14 = sum14 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum14 = sum14 + apot
                    else:
                        #print ("0", file=f1)
                        sum14 = sum14 + 0
            else:
                #print ("0", file=f1)
                sum14 = sum14 + 0



        elif (data[i][3] == "15_Repetitive/CNV") :
            if (data[i][0]==der[j][0]):
                    if (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][2])) :
                        sum15 = sum15 + 1

                    elif (int(data[i][1]) < int(der[j][1])) and (int(data[i][2]) > int(der[j][1])) and (int(data[i][2]) <= int(der[j][2])) :
                        apot = float(int(data[i][2])- int(der[j][1]))/2000
                        sum15 = sum15 + apot
                          # print(apot, file=f1)
                    elif (int(data[i][1]) >= int(der[j][1])) and (int(data[i][1]) < int(der[j][2])) :
                        if (int(data[i][2]) < int(der[j][2])) :
                            apot = float(int(data[i][2])- int(data[i][1]))/2000
                            sum15 = sum15 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data[i][1]))/2000
                            #print(apot, file=f1)
                            sum15 = sum15 + apot
                    else:
                        #print ("0", file=f1)
                        sum15 = sum15 + 0
            else:
                #print ("0", file=f1)
                sum15 = sum15 + 0




    print (sum1,"\t",end="",file = f1)
    print (sum2,"\t",end="",file = f1)
    print (sum3,"\t",end="",file = f1)
    print (sum4,"\t",end="",file = f1)
    print (sum5,"\t",end="",file = f1)
    print (sum6,"\t",end="",file = f1)
    print (sum7,"\t",end="",file = f1)
    print (sum8,"\t",end="",file = f1)
    print (sum9,"\t",end="",file = f1)
    print (sum10,"\t",end="",file = f1)
    print (sum11,"\t",end="",file = f1)
    print (sum12,"\t",end="",file = f1)
    print (sum13,"\t",end="",file = f1)
    print (sum14,"\t",end="",file = f1)
    print (sum15,"\t",end="",file = f1)

    sum1=0
    sum2=0
    sum3=0
    sum4=0
    sum5=0
    sum6=0
    sum7=0
    sum8=0
    sum9=0
    sum10=0
    sum11=0
    sum12=0
    sum13=0
    sum14=0
    sum15=0
    i=0



    for i in xrange(x2):
        if (data2[i][3] == "1_Active_Promoter") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum1 = sum1 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum1 = sum1 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum1 = sum1 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum1 = sum1 + apot
                    else:
                        #print ("0", file=f1)
                        sum1 = sum1 + 0
            else:
                #print ("0", file=f1)
                sum1 = sum1 + 0


        elif (data2[i][3] == "2_Weak_Promoter") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum2 = sum2 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum2 = sum2 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum2 = sum2 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum2 = sum2 + apot
                    else:
                        #print ("0", file=f1)
                        sum2 = sum2 + 0
            else:
                #print ("0", file=f1)
                sum2 = sum2 + 0


    
        elif (data2[i][3] == "3_Poised_Promoter") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum3 = sum3 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum3 = sum3 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum3 = sum3 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum3 = sum3 + apot
                    else:
                        #print ("0", file=f1)
                        sum3 = sum3 + 0
            else:
                #print ("0", file=f1)
                sum3 = sum3 + 0


        elif (data2[i][3] == "4_Strong_Enhancer") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum4 = sum4 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum4 = sum4 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum4 = sum4 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum4 = sum4 + apot
                    else:
                        #print ("0", file=f1)
                        sum4 = sum4 + 0
            else:
                #print ("0", file=f1)
                sum4 = sum4 + 0


        elif (data2[i][3] == "5_Strong_Enhancer") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum5 = sum5 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum5 = sum5 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum5 = sum5 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum5 = sum5 + apot
                    else:
                        #print ("0", file=f1)
                        sum5 = sum5 + 0
            else:
                #print ("0", file=f1)
                sum5 = sum5 + 0


        elif (data2[i][3] == "6_Weak_Enhancer") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum6 = sum6 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum6 = sum6 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum6 = sum6 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum6 = sum6 + apot
                    else:
                        #print ("0", file=f1)
                        sum6 = sum6 + 0
            else:
                #print ("0", file=f1)
                sum6 = sum6 + 0



        elif (data2[i][3] == "7_Weak_Enhancer") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum7 = sum7 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum7 = sum7 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum7 = sum7 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum7 = sum7 + apot
                    else:
                        #print ("0", file=f1)
                        sum7 = sum7 + 0
            else:
                #print ("0", file=f1)
                sum7 = sum7 + 0



        elif (data2[i][3] == "8_Insulator") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum8 = sum8 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum8 = sum8 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum8 = sum8 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum8 = sum8 + apot
                    else:
                        #print ("0", file=f1)
                        sum8 = sum8 + 0
            else:
                #print ("0", file=f1)
                sum8 = sum8 + 0



        elif (data2[i][3] == "9_Txn_Transition") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum9 = sum9 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum9 = sum9 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum9 = sum9 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum9 = sum9 + apot
                    else:
                        #print ("0", file=f1)
                        sum9 = sum9 + 0
            else:
                #print ("0", file=f1)
                sum9 = sum9 + 0


        elif (data2[i][3] == "10_Txn_Elongation") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum10 = sum10 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum10 = sum10 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum10 = sum10 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum10 = sum10 + apot
                    else:
                        #print ("0", file=f1)
                        sum10 = sum10 + 0
            else:
                #print ("0", file=f1)
                sum10 = sum10 + 0


        elif (data2[i][3] == "11_Weak_Txn") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum11 = sum11 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum11 = sum11 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum11 = sum11 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum11 = sum11 + apot
                    else:
                        #print ("0", file=f1)
                        sum11 = sum11 + 0
            else:
                #print ("0", file=f1)
                sum11 = sum11 + 0


        elif (data2[i][3] == "12_Repressed") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum12 = sum12 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum12 = sum12 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum12 = sum12 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum12 = sum12 + apot
                    else:
                        #print ("0", file=f1)
                        sum12 = sum12 + 0
            else:
                #print ("0", file=f1)
                sum12 = sum12 + 0


        elif (data2[i][3] == "13_Heterochrom/lo") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum13 = sum13 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum13 = sum13 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum13 = sum13 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum13 = sum13 + apot
                    else:
                        #print ("0", file=f1)
                        sum13 = sum13 + 0
            else:
                #print ("0", file=f1)
                sum13 = sum13 + 0


        elif (data2[i][3] == "14_Repetitive/CNV") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum14 = sum14 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum14 = sum14 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum14 = sum14 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum14 = sum14 + apot
                    else:
                        #print ("0", file=f1)
                        sum14 = sum14 + 0
            else:
                #print ("0", file=f1)
                sum14 = sum14 + 0



        elif (data2[i][3] == "15_Repetitive/CNV") :
            if (data2[i][0]==der[j][0]):
                    if (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][2])) :
                        sum15 = sum15 + 1

                    elif (int(data2[i][1]) < int(der[j][1])) and (int(data2[i][2]) > int(der[j][1])) and (int(data2[i][2]) <= int(der[j][2])) :
                        apot = float(int(data2[i][2])- int(der[j][1]))/2000
                        sum15 = sum15 + apot
                          # print(apot, file=f1)
                    elif (int(data2[i][1]) >= int(der[j][1])) and (int(data2[i][1]) < int(der[j][2])) :
                        if (int(data2[i][2]) < int(der[j][2])) :
                            apot = float(int(data2[i][2])- int(data2[i][1]))/2000
                            sum15 = sum15 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data2[i][1]))/2000
                            #print(apot, file=f1)
                            sum15 = sum15 + apot
                    else:
                        #print ("0", file=f1)
                        sum15 = sum15 + 0
            else:
                #print ("0", file=f1)
                sum15 = sum15 + 0




    print (sum1,"\t",end="",file = f1)
    print (sum2,"\t",end="",file = f1)
    print (sum3,"\t",end="",file = f1)
    print (sum4,"\t",end="",file = f1)
    print (sum5,"\t",end="",file = f1)
    print (sum6,"\t",end="",file = f1)
    print (sum7,"\t",end="",file = f1)
    print (sum8,"\t",end="",file = f1)
    print (sum9,"\t",end="",file = f1)
    print (sum10,"\t",end="",file = f1)
    print (sum11,"\t",end="",file = f1)
    print (sum12,"\t",end="",file = f1)
    print (sum13,"\t",end="",file = f1)
    print (sum14,"\t",end="",file = f1)
    print (sum15,"\t",end="",file = f1)

    sum1=0
    sum2=0
    sum3=0
    sum4=0
    sum5=0
    sum6=0
    sum7=0
    sum8=0
    sum9=0
    sum10=0
    sum11=0
    sum12=0
    sum13=0
    sum14=0
    sum15=0
    i=0


    for i in xrange(x3):
        if (data3[i][3] == "1_Active_Promoter") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum1 = sum1 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum1 = sum1 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum1 = sum1 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum1 = sum1 + apot
                    else:
                        #print ("0", file=f1)
                        sum1 = sum1 + 0
            else:
                #print ("0", file=f1)
                sum1 = sum1 + 0


        elif (data3[i][3] == "2_Weak_Promoter") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum2 = sum2 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum2 = sum2 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum2 = sum2 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum2 = sum2 + apot
                    else:
                        #print ("0", file=f1)
                        sum2 = sum2 + 0
            else:
                #print ("0", file=f1)
                sum2 = sum2 + 0


    
        elif (data3[i][3] == "3_Poised_Promoter") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum3 = sum3 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum3 = sum3 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum3 = sum3 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum3 = sum3 + apot
                    else:
                        #print ("0", file=f1)
                        sum3 = sum3 + 0
            else:
                #print ("0", file=f1)
                sum3 = sum3 + 0


        elif (data3[i][3] == "4_Strong_Enhancer") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum4 = sum4 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum4 = sum4 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum4 = sum4 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum4 = sum4 + apot
                    else:
                        #print ("0", file=f1)
                        sum4 = sum4 + 0
            else:
                #print ("0", file=f1)
                sum4 = sum4 + 0


        elif (data3[i][3] == "5_Strong_Enhancer") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum5 = sum5 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum5 = sum5 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum5 = sum5 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum5 = sum5 + apot
                    else:
                        #print ("0", file=f1)
                        sum5 = sum5 + 0
            else:
                #print ("0", file=f1)
                sum5 = sum5 + 0


        elif (data3[i][3] == "6_Weak_Enhancer") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum6 = sum6 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum6 = sum6 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum6 = sum6 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum6 = sum6 + apot
                    else:
                        #print ("0", file=f1)
                        sum6 = sum6 + 0
            else:
                #print ("0", file=f1)
                sum6 = sum6 + 0



        elif (data3[i][3] == "7_Weak_Enhancer") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum7 = sum7 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum7 = sum7 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum7 = sum7 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum7 = sum7 + apot
                    else:
                        #print ("0", file=f1)
                        sum7 = sum7 + 0
            else:
                #print ("0", file=f1)
                sum7 = sum7 + 0



        elif (data3[i][3] == "8_Insulator") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum8 = sum8 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum8 = sum8 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum8 = sum8 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum8 = sum8 + apot
                    else:
                        #print ("0", file=f1)
                        sum8 = sum8 + 0
            else:
                #print ("0", file=f1)
                sum8 = sum8 + 0



        elif (data3[i][3] == "9_Txn_Transition") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum9 = sum9 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum9 = sum9 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum9 = sum9 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum9 = sum9 + apot
                    else:
                        #print ("0", file=f1)
                        sum9 = sum9 + 0
            else:
                #print ("0", file=f1)
                sum9 = sum9 + 0


        elif (data3[i][3] == "10_Txn_Elongation") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum10 = sum10 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum10 = sum10 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum10 = sum10 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum10 = sum10 + apot
                    else:
                        #print ("0", file=f1)
                        sum10 = sum10 + 0
            else:
                #print ("0", file=f1)
                sum10 = sum10 + 0


        elif (data3[i][3] == "11_Weak_Txn") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum11 = sum11 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum11 = sum11 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum11 = sum11 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum11 = sum11 + apot
                    else:
                        #print ("0", file=f1)
                        sum11 = sum11 + 0
            else:
                #print ("0", file=f1)
                sum11 = sum11 + 0


        elif (data3[i][3] == "12_Repressed") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum12 = sum12 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum12 = sum12 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum12 = sum12 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum12 = sum12 + apot
                    else:
                        #print ("0", file=f1)
                        sum12 = sum12 + 0
            else:
                #print ("0", file=f1)
                sum12 = sum12 + 0


        elif (data3[i][3] == "13_Heterochrom/lo") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum13 = sum13 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum13 = sum13 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum13 = sum13 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum13 = sum13 + apot
                    else:
                        #print ("0", file=f1)
                        sum13 = sum13 + 0
            else:
                #print ("0", file=f1)
                sum13 = sum13 + 0


        elif (data3[i][3] == "14_Repetitive/CNV") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum14 = sum14 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum14 = sum14 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum14 = sum14 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum14 = sum14 + apot
                    else:
                        #print ("0", file=f1)
                        sum14 = sum14 + 0
            else:
                #print ("0", file=f1)
                sum14 = sum14 + 0



        elif (data3[i][3] == "15_Repetitive/CNV") :
            if (data3[i][0]==der[j][0]):
                    if (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][2])) :
                        sum15 = sum15 + 1

                    elif (int(data3[i][1]) < int(der[j][1])) and (int(data3[i][2]) > int(der[j][1])) and (int(data3[i][2]) <= int(der[j][2])) :
                        apot = float(int(data3[i][2])- int(der[j][1]))/2000
                        sum15 = sum15 + apot
                          # print(apot, file=f1)
                    elif (int(data3[i][1]) >= int(der[j][1])) and (int(data3[i][1]) < int(der[j][2])) :
                        if (int(data3[i][2]) < int(der[j][2])) :
                            apot = float(int(data3[i][2])- int(data3[i][1]))/2000
                            sum15 = sum15 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data3[i][1]))/2000
                            #print(apot, file=f1)
                            sum15 = sum15 + apot
                    else:
                        #print ("0", file=f1)
                        sum15 = sum15 + 0
            else:
                #print ("0", file=f1)
                sum15 = sum15 + 0




    print (sum1,"\t",end="",file = f1)
    print (sum2,"\t",end="",file = f1)
    print (sum3,"\t",end="",file = f1)
    print (sum4,"\t",end="",file = f1)
    print (sum5,"\t",end="",file = f1)
    print (sum6,"\t",end="",file = f1)
    print (sum7,"\t",end="",file = f1)
    print (sum8,"\t",end="",file = f1)
    print (sum9,"\t",end="",file = f1)
    print (sum10,"\t",end="",file = f1)
    print (sum11,"\t",end="",file = f1)
    print (sum12,"\t",end="",file = f1)
    print (sum13,"\t",end="",file = f1)
    print (sum14,"\t",end="",file = f1)
    print (sum15,"\t",end="",file = f1)

    sum1=0
    sum2=0
    sum3=0
    sum4=0
    sum5=0
    sum6=0
    sum7=0
    sum8=0
    sum9=0
    sum10=0
    sum11=0
    sum12=0
    sum13=0
    sum14=0
    sum15=0
    i=0



    for i in xrange(x4):
        if (data4[i][3] == "1_Active_Promoter") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum1 = sum1 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum1 = sum1 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum1 = sum1 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum1 = sum1 + apot
                    else:
                        #print ("0", file=f1)
                        sum1 = sum1 + 0
            else:
                #print ("0", file=f1)
                sum1 = sum1 + 0

    


        elif (data4[i][3] == "2_Weak_Promoter") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum2 = sum2 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum2 = sum2 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum2 = sum2 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum2 = sum2 + apot
                    else:
                        #print ("0", file=f1)
                        sum2 = sum2 + 0
            else:
                #print ("0", file=f1)
                sum2 = sum2 + 0


    
        elif (data4[i][3] == "3_Poised_Promoter") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum3 = sum3 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum3 = sum3 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum3 = sum3 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum3 = sum3 + apot
                    else:
                        #print ("0", file=f1)
                        sum3 = sum3 + 0
            else:
                #print ("0", file=f1)
                sum3 = sum3 + 0


        elif (data4[i][3] == "4_Strong_Enhancer") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum4 = sum4 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum4 = sum4 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum4 = sum4 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum4 = sum4 + apot
                    else:
                        #print ("0", file=f1)
                        sum4 = sum4 + 0
            else:
                #print ("0", file=f1)
                sum4 = sum4 + 0


        elif (data4[i][3] == "5_Strong_Enhancer") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum5 = sum5 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum5 = sum5 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum5 = sum5 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum5 = sum5 + apot
                    else:
                        #print ("0", file=f1)
                        sum5 = sum5 + 0
            else:
                #print ("0", file=f1)
                sum5 = sum5 + 0


        elif (data4[i][3] == "6_Weak_Enhancer") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum6 = sum6 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum6 = sum6 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum6 = sum6 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum6 = sum6 + apot
                    else:
                        #print ("0", file=f1)
                        sum6 = sum6 + 0
            else:
                #print ("0", file=f1)
                sum6 = sum6 + 0



        elif (data4[i][3] == "7_Weak_Enhancer") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum7 = sum7 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum7 = sum7 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum7 = sum7 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum7 = sum7 + apot
                    else:
                        #print ("0", file=f1)
                        sum7 = sum7 + 0
            else:
                #print ("0", file=f1)
                sum7 = sum7 + 0



        elif (data4[i][3] == "8_Insulator") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum8 = sum8 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum8 = sum8 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum8 = sum8 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum8 = sum8 + apot
                    else:
                        #print ("0", file=f1)
                        sum8 = sum8 + 0
            else:
                #print ("0", file=f1)
                sum8 = sum8 + 0





        elif (data4[i][3] == "9_Txn_Transition") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum9 = sum9 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum9 = sum9 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum9 = sum9 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum9 = sum9 + apot
                    else:
                        #print ("0", file=f1)
                        sum9 = sum9 + 0
            else:
                #print ("0", file=f1)
                sum9 = sum9 + 0




        elif (data4[i][3] == "10_Txn_Elongation") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum10 = sum10 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum10 = sum10 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum10 = sum10 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum10 = sum10 + apot
                    else:
                        #print ("0", file=f1)
                        sum10 = sum10 + 0
            else:
                #print ("0", file=f1)
                sum10 = sum10 + 0


        elif (data4[i][3] == "11_Weak_Txn") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum11 = sum11 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum11 = sum11 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum11 = sum11 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum11 = sum11 + apot
                    else:
                        #print ("0", file=f1)
                        sum11 = sum11 + 0
            else:
                #print ("0", file=f1)
                sum11 = sum11 + 0


        elif (data4[i][3] == "12_Repressed") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum12 = sum12 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum12 = sum12 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum12 = sum12 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum12 = sum12 + apot
                    else:
                        #print ("0", file=f1)
                        sum12 = sum12 + 0
            else:
                #print ("0", file=f1)
                sum12 = sum12 + 0


        elif (data4[i][3] == "13_Heterochrom/lo") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum13 = sum13 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum13 = sum13 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum13 = sum13 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum13 = sum13 + apot
                    else:
                        #print ("0", file=f1)
                        sum13 = sum13 + 0
            else:
                #print ("0", file=f1)
                sum13 = sum13 + 0


        elif (data4[i][3] == "14_Repetitive/CNV") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum14 = sum14 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum14 = sum14 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum14 = sum14 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum14 = sum14 + apot
                    else:
                        #print ("0", file=f1)
                        sum14 = sum14 + 0
            else:
                #print ("0", file=f1)
                sum14 = sum14 + 0



        elif (data4[i][3] == "15_Repetitive/CNV") :
            if (data4[i][0]==der[j][0]):
                    if (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][2])) :
                        sum15 = sum15 + 1

                    elif (int(data4[i][1]) < int(der[j][1])) and (int(data4[i][2]) > int(der[j][1])) and (int(data4[i][2]) <= int(der[j][2])) :
                        apot = float(int(data4[i][2])- int(der[j][1]))/2000
                        sum15 = sum15 + apot
                          # print(apot, file=f1)
                    elif (int(data4[i][1]) >= int(der[j][1])) and (int(data4[i][1]) < int(der[j][2])) :
                        if (int(data4[i][2]) < int(der[j][2])) :
                            apot = float(int(data4[i][2])- int(data4[i][1]))/2000
                            sum15 = sum15 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data4[i][1]))/2000
                            #print(apot, file=f1)
                            sum15 = sum15 + apot
                    else:
                        #print ("0", file=f1)
                        sum15 = sum15 + 0
            else:
                #print ("0", file=f1)
                sum15 = sum15 + 0




    print (sum1,"\t",end="",file = f1)
    print (sum2,"\t",end="",file = f1)
    print (sum3,"\t",end="",file = f1)
    print (sum4,"\t",end="",file = f1)
    print (sum5,"\t",end="",file = f1)
    print (sum6,"\t",end="",file = f1)
    print (sum7,"\t",end="",file = f1)
    print (sum8,"\t",end="",file = f1)
    print (sum9,"\t",end="",file = f1)
    print (sum10,"\t",end="",file = f1)
    print (sum11,"\t",end="",file = f1)
    print (sum12,"\t",end="",file = f1)
    print (sum13,"\t",end="",file = f1)
    print (sum14,"\t",end="",file = f1)
    print (sum15,"\t",end="",file = f1)

    sum1=0
    sum2=0
    sum3=0
    sum4=0
    sum5=0
    sum6=0
    sum7=0
    sum8=0
    sum9=0
    sum10=0
    sum11=0
    sum12=0
    sum13=0
    sum14=0
    sum15=0
    i=0



    for i in xrange(x5):
        if (data5[i][3] == "1_Active_Promoter") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum1 = sum1 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum1 = sum1 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum1 = sum1 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum1 = sum1 + apot
                    else:
                        #print ("0", file=f1)
                        sum1 = sum1 + 0
            else:
                #print ("0", file=f1)
                sum1 = sum1 + 0

    


        elif (data5[i][3] == "2_Weak_Promoter") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum2 = sum2 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum2 = sum2 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum2 = sum2 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum2 = sum2 + apot
                    else:
                        #print ("0", file=f1)
                        sum2 = sum2 + 0
            else:
                #print ("0", file=f1)
                sum2 = sum2 + 0


    
        elif (data5[i][3] == "3_Poised_Promoter") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum3 = sum3 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum3 = sum3 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum3 = sum3 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum3 = sum3 + apot
                    else:
                        #print ("0", file=f1)
                        sum3 = sum3 + 0
            else:
                #print ("0", file=f1)
                sum3 = sum3 + 0


        elif (data5[i][3] == "4_Strong_Enhancer") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum4 = sum4 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum4 = sum4 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum4 = sum4 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum4 = sum4 + apot
                    else:
                        #print ("0", file=f1)
                        sum4 = sum4 + 0
            else:
                #print ("0", file=f1)
                sum4 = sum4 + 0


        elif (data5[i][3] == "5_Strong_Enhancer") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum5 = sum5 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum5 = sum5 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum5 = sum5 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum5 = sum5 + apot
                    else:
                        #print ("0", file=f1)
                        sum5 = sum5 + 0
            else:
                #print ("0", file=f1)
                sum5 = sum5 + 0


        elif (data5[i][3] == "6_Weak_Enhancer") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum6 = sum6 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum6 = sum6 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum6 = sum6 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum6 = sum6 + apot
                    else:
                        #print ("0", file=f1)
                        sum6 = sum6 + 0
            else:
                #print ("0", file=f1)
                sum6 = sum6 + 0



        elif (data5[i][3] == "7_Weak_Enhancer") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum7 = sum7 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum7 = sum7 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum7 = sum7 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum7 = sum7 + apot
                    else:
                        #print ("0", file=f1)
                        sum7 = sum7 + 0
            else:
                #print ("0", file=f1)
                sum7 = sum7 + 0



        elif (data5[i][3] == "8_Insulator") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum8 = sum8 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum8 = sum8 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum8 = sum8 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum8 = sum8 + apot
                    else:
                        #print ("0", file=f1)
                        sum8 = sum8 + 0
            else:
                #print ("0", file=f1)
                sum8 = sum8 + 0





        elif (data5[i][3] == "9_Txn_Transition") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum9 = sum9 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum9 = sum9 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum9 = sum9 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum9 = sum9 + apot
                    else:
                        #print ("0", file=f1)
                        sum9 = sum9 + 0
            else:
                #print ("0", file=f1)
                sum9 = sum9 + 0




        elif (data5[i][3] == "10_Txn_Elongation") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum10 = sum10 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum10 = sum10 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum10 = sum10 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum10 = sum10 + apot
                    else:
                        #print ("0", file=f1)
                        sum10 = sum10 + 0
            else:
                #print ("0", file=f1)
                sum10 = sum10 + 0


        elif (data5[i][3] == "11_Weak_Txn") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum11 = sum11 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum11 = sum11 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum11 = sum11 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum11 = sum11 + apot
                    else:
                        #print ("0", file=f1)
                        sum11 = sum11 + 0
            else:
                #print ("0", file=f1)
                sum11 = sum11 + 0


        elif (data5[i][3] == "12_Repressed") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum12 = sum12 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum12 = sum12 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum12 = sum12 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum12 = sum12 + apot
                    else:
                        #print ("0", file=f1)
                        sum12 = sum12 + 0
            else:
                #print ("0", file=f1)
                sum12 = sum12 + 0


        elif (data5[i][3] == "13_Heterochrom/lo") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum13 = sum13 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum13 = sum13 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum13 = sum13 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum13 = sum13 + apot
                    else:
                        #print ("0", file=f1)
                        sum13 = sum13 + 0
            else:
                #print ("0", file=f1)
                sum13 = sum13 + 0


        elif (data5[i][3] == "14_Repetitive/CNV") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum14 = sum14 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum14 = sum14 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum14 = sum14 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum14 = sum14 + apot
                    else:
                        #print ("0", file=f1)
                        sum14 = sum14 + 0
            else:
                #print ("0", file=f1)
                sum14 = sum14 + 0



        elif (data5[i][3] == "15_Repetitive/CNV") :
            if (data5[i][0]==der[j][0]):
                    if (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][2])) :
                        sum15 = sum15 + 1

                    elif (int(data5[i][1]) < int(der[j][1])) and (int(data5[i][2]) > int(der[j][1])) and (int(data5[i][2]) <= int(der[j][2])) :
                        apot = float(int(data5[i][2])- int(der[j][1]))/2000
                        sum15 = sum15 + apot
                          # print(apot, file=f1)
                    elif (int(data5[i][1]) >= int(der[j][1])) and (int(data5[i][1]) < int(der[j][2])) :
                        if (int(data5[i][2]) < int(der[j][2])) :
                            apot = float(int(data5[i][2])- int(data5[i][1]))/2000
                            sum15 = sum15 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data5[i][1]))/2000
                            #print(apot, file=f1)
                            sum15 = sum15 + apot
                    else:
                        #print ("0", file=f1)
                        sum15 = sum15 + 0
            else:
                #print ("0", file=f1)
                sum15 = sum15 + 0




    print (sum1,"\t",end="",file = f1)
    print (sum2,"\t",end="",file = f1)
    print (sum3,"\t",end="",file = f1)
    print (sum4,"\t",end="",file = f1)
    print (sum5,"\t",end="",file = f1)
    print (sum6,"\t",end="",file = f1)
    print (sum7,"\t",end="",file = f1)
    print (sum8,"\t",end="",file = f1)
    print (sum9,"\t",end="",file = f1)
    print (sum10,"\t",end="",file = f1)
    print (sum11,"\t",end="",file = f1)
    print (sum12,"\t",end="",file = f1)
    print (sum13,"\t",end="",file = f1)
    print (sum14,"\t",end="",file = f1)
    print (sum15,"\t",end="",file = f1)

    sum1=0
    sum2=0
    sum3=0
    sum4=0
    sum5=0
    sum6=0
    sum7=0
    sum8=0
    sum9=0
    sum10=0
    sum11=0
    sum12=0
    sum13=0
    sum14=0
    sum15=0
    i=0


    for i in xrange(x6):
        if (data6[i][3] == "1_Active_Promoter") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum1 = sum1 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum1 = sum1 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum1 = sum1 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum1 = sum1 + apot
                    else:
                        #print ("0", file=f1)
                        sum1 = sum1 + 0
            else:
                #print ("0", file=f1)
                sum1 = sum1 + 0

    


        elif (data6[i][3] == "2_Weak_Promoter") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum2 = sum2 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum2 = sum2 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum2 = sum2 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum2 = sum2 + apot
                    else:
                        #print ("0", file=f1)
                        sum2 = sum2 + 0
            else:
                #print ("0", file=f1)
                sum2 = sum2 + 0


    
        elif (data6[i][3] == "3_Poised_Promoter") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum3 = sum3 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum3 = sum3 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum3 = sum3 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum3 = sum3 + apot
                    else:
                        #print ("0", file=f1)
                        sum3 = sum3 + 0
            else:
                #print ("0", file=f1)
                sum3 = sum3 + 0


        elif (data6[i][3] == "4_Strong_Enhancer") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum4 = sum4 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum4 = sum4 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum4 = sum4 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum4 = sum4 + apot
                    else:
                        #print ("0", file=f1)
                        sum4 = sum4 + 0
            else:
                #print ("0", file=f1)
                sum4 = sum4 + 0


        elif (data6[i][3] == "5_Strong_Enhancer") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum5 = sum5 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum5 = sum5 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum5 = sum5 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum5 = sum5 + apot
                    else:
                        #print ("0", file=f1)
                        sum5 = sum5 + 0
            else:
                #print ("0", file=f1)
                sum5 = sum5 + 0


        elif (data6[i][3] == "6_Weak_Enhancer") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum6 = sum6 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum6 = sum6 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum6 = sum6 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum6 = sum6 + apot
                    else:
                        #print ("0", file=f1)
                        sum6 = sum6 + 0
            else:
                #print ("0", file=f1)
                sum6 = sum6 + 0



        elif (data6[i][3] == "7_Weak_Enhancer") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum7 = sum7 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum7 = sum7 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum7 = sum7 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum7 = sum7 + apot
                    else:
                        #print ("0", file=f1)
                        sum7 = sum7 + 0
            else:
                #print ("0", file=f1)
                sum7 = sum7 + 0



        elif (data6[i][3] == "8_Insulator") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum8 = sum8 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum8 = sum8 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum8 = sum8 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum8 = sum8 + apot
                    else:
                        #print ("0", file=f1)
                        sum8 = sum8 + 0
            else:
                #print ("0", file=f1)
                sum8 = sum8 + 0





        elif (data6[i][3] == "9_Txn_Transition") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum9 = sum9 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum9 = sum9 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum9 = sum9 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum9 = sum9 + apot
                    else:
                        #print ("0", file=f1)
                        sum9 = sum9 + 0
            else:
                #print ("0", file=f1)
                sum9 = sum9 + 0




        elif (data6[i][3] == "10_Txn_Elongation") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum10 = sum10 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum10 = sum10 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum10 = sum10 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum10 = sum10 + apot
                    else:
                        #print ("0", file=f1)
                        sum10 = sum10 + 0
            else:
                #print ("0", file=f1)
                sum10 = sum10 + 0


        elif (data6[i][3] == "11_Weak_Txn") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum11 = sum11 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum11 = sum11 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum11 = sum11 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum11 = sum11 + apot
                    else:
                        #print ("0", file=f1)
                        sum11 = sum11 + 0
            else:
                #print ("0", file=f1)
                sum11 = sum11 + 0


        elif (data6[i][3] == "12_Repressed") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum12 = sum12 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum12 = sum12 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum12 = sum12 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum12 = sum12 + apot
                    else:
                        #print ("0", file=f1)
                        sum12 = sum12 + 0
            else:
                #print ("0", file=f1)
                sum12 = sum12 + 0


        elif (data6[i][3] == "13_Heterochrom/lo") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum13 = sum13 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum13 = sum13 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum13 = sum13 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum13 = sum13 + apot
                    else:
                        #print ("0", file=f1)
                        sum13 = sum13 + 0
            else:
                #print ("0", file=f1)
                sum13 = sum13 + 0


        elif (data6[i][3] == "14_Repetitive/CNV") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum14 = sum14 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum14 = sum14 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum14 = sum14 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum14 = sum14 + apot
                    else:
                        #print ("0", file=f1)
                        sum14 = sum14 + 0
            else:
                #print ("0", file=f1)
                sum14 = sum14 + 0



        elif (data6[i][3] == "15_Repetitive/CNV") :
            if (data6[i][0]==der[j][0]):
                    if (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][2])) :
                        sum15 = sum15 + 1

                    elif (int(data6[i][1]) < int(der[j][1])) and (int(data6[i][2]) > int(der[j][1])) and (int(data6[i][2]) <= int(der[j][2])) :
                        apot = float(int(data6[i][2])- int(der[j][1]))/2000
                        sum15 = sum15 + apot
                          # print(apot, file=f1)
                    elif (int(data6[i][1]) >= int(der[j][1])) and (int(data6[i][1]) < int(der[j][2])) :
                        if (int(data6[i][2]) < int(der[j][2])) :
                            apot = float(int(data6[i][2])- int(data6[i][1]))/2000
                            sum15 = sum15 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data6[i][1]))/2000
                            #print(apot, file=f1)
                            sum15 = sum15 + apot
                    else:
                        #print ("0", file=f1)
                        sum15 = sum15 + 0
            else:
                #print ("0", file=f1)
                sum15 = sum15 + 0




    print (sum1,"\t",end="",file = f1)
    print (sum2,"\t",end="",file = f1)
    print (sum3,"\t",end="",file = f1)
    print (sum4,"\t",end="",file = f1)
    print (sum5,"\t",end="",file = f1)
    print (sum6,"\t",end="",file = f1)
    print (sum7,"\t",end="",file = f1)
    print (sum8,"\t",end="",file = f1)
    print (sum9,"\t",end="",file = f1)
    print (sum10,"\t",end="",file = f1)
    print (sum11,"\t",end="",file = f1)
    print (sum12,"\t",end="",file = f1)
    print (sum13,"\t",end="",file = f1)
    print (sum14,"\t",end="",file = f1)
    print (sum15,"\t",end="",file = f1)

    sum1=0
    sum2=0
    sum3=0
    sum4=0
    sum5=0
    sum6=0
    sum7=0
    sum8=0
    sum9=0
    sum10=0
    sum11=0
    sum12=0
    sum13=0
    sum14=0
    sum15=0
    i=0


    for i in xrange(x4):
        if (data7[i][3] == "1_Active_Promoter") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum1 = sum1 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum1 = sum1 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum1 = sum1 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum1 = sum1 + apot
                    else:
                        #print ("0", file=f1)
                        sum1 = sum1 + 0
            else:
                #print ("0", file=f1)
                sum1 = sum1 + 0

    


        elif (data7[i][3] == "2_Weak_Promoter") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum2 = sum2 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum2 = sum2 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum2 = sum2 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum2 = sum2 + apot
                    else:
                        #print ("0", file=f1)
                        sum2 = sum2 + 0
            else:
                #print ("0", file=f1)
                sum2 = sum2 + 0


    
        elif (data7[i][3] == "3_Poised_Promoter") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum3 = sum3 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum3 = sum3 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum3 = sum3 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum3 = sum3 + apot
                    else:
                        #print ("0", file=f1)
                        sum3 = sum3 + 0
            else:
                #print ("0", file=f1)
                sum3 = sum3 + 0


        elif (data7[i][3] == "4_Strong_Enhancer") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum4 = sum4 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum4 = sum4 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum4 = sum4 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum4 = sum4 + apot
                    else:
                        #print ("0", file=f1)
                        sum4 = sum4 + 0
            else:
                #print ("0", file=f1)
                sum4 = sum4 + 0


        elif (data7[i][3] == "5_Strong_Enhancer") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum5 = sum5 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum5 = sum5 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum5 = sum5 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum5 = sum5 + apot
                    else:
                        #print ("0", file=f1)
                        sum5 = sum5 + 0
            else:
                #print ("0", file=f1)
                sum5 = sum5 + 0


        elif (data7[i][3] == "6_Weak_Enhancer") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum6 = sum6 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum6 = sum6 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum6 = sum6 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum6 = sum6 + apot
                    else:
                        #print ("0", file=f1)
                        sum6 = sum6 + 0
            else:
                #print ("0", file=f1)
                sum6 = sum6 + 0



        elif (data7[i][3] == "7_Weak_Enhancer") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum7 = sum7 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum7 = sum7 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum7 = sum7 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum7 = sum7 + apot
                    else:
                        #print ("0", file=f1)
                        sum7 = sum7 + 0
            else:
                #print ("0", file=f1)
                sum7 = sum7 + 0



        elif (data7[i][3] == "8_Insulator") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum8 = sum8 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum8 = sum8 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum8 = sum8 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum8 = sum8 + apot
                    else:
                        #print ("0", file=f1)
                        sum8 = sum8 + 0
            else:
                #print ("0", file=f1)
                sum8 = sum8 + 0





        elif (data7[i][3] == "9_Txn_Transition") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum9 = sum9 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum9 = sum9 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum9 = sum9 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum9 = sum9 + apot
                    else:
                        #print ("0", file=f1)
                        sum9 = sum9 + 0
            else:
                #print ("0", file=f1)
                sum9 = sum9 + 0




        elif (data7[i][3] == "10_Txn_Elongation") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum10 = sum10 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum10 = sum10 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum10 = sum10 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum10 = sum10 + apot
                    else:
                        #print ("0", file=f1)
                        sum10 = sum10 + 0
            else:
                #print ("0", file=f1)
                sum10 = sum10 + 0


        elif (data7[i][3] == "11_Weak_Txn") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum11 = sum11 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum11 = sum11 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum11 = sum11 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum11 = sum11 + apot
                    else:
                        #print ("0", file=f1)
                        sum11 = sum11 + 0
            else:
                #print ("0", file=f1)
                sum11 = sum11 + 0


        elif (data7[i][3] == "12_Repressed") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum12 = sum12 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum12 = sum12 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum12 = sum12 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum12 = sum12 + apot
                    else:
                        #print ("0", file=f1)
                        sum12 = sum12 + 0
            else:
                #print ("0", file=f1)
                sum12 = sum12 + 0


        elif (data7[i][3] == "13_Heterochrom/lo") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum13 = sum13 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum13 = sum13 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum13 = sum13 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum13 = sum13 + apot
                    else:
                        #print ("0", file=f1)
                        sum13 = sum13 + 0
            else:
                #print ("0", file=f1)
                sum13 = sum13 + 0


        elif (data7[i][3] == "14_Repetitive/CNV") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum14 = sum14 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum14 = sum14 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum14 = sum14 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum14 = sum14 + apot
                    else:
                        #print ("0", file=f1)
                        sum14 = sum14 + 0
            else:
                #print ("0", file=f1)
                sum14 = sum14 + 0



        elif (data7[i][3] == "15_Repetitive/CNV") :
            if (data7[i][0]==der[j][0]):
                    if (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][2])) :
                        sum15 = sum15 + 1

                    elif (int(data7[i][1]) < int(der[j][1])) and (int(data7[i][2]) > int(der[j][1])) and (int(data7[i][2]) <= int(der[j][2])) :
                        apot = float(int(data7[i][2])- int(der[j][1]))/2000
                        sum15 = sum15 + apot
                          # print(apot, file=f1)
                    elif (int(data7[i][1]) >= int(der[j][1])) and (int(data7[i][1]) < int(der[j][2])) :
                        if (int(data7[i][2]) < int(der[j][2])) :
                            apot = float(int(data7[i][2])- int(data7[i][1]))/2000
                            sum15 = sum15 + apot
                            #print(apot, file=f1)
                        else:
                            apot = float(int(der[j][2])- int(data7[i][1]))/2000
                            #print(apot, file=f1)
                            sum15 = sum15 + apot
                    else:
                        #print ("0", file=f1)
                        sum15 = sum15 + 0
            else:
                #print ("0", file=f1)
                sum15 = sum15 + 0




    print (sum1,"\t",end="",file = f1)
    print (sum2,"\t",end="",file = f1)
    print (sum3,"\t",end="",file = f1)
    print (sum4,"\t",end="",file = f1)
    print (sum5,"\t",end="",file = f1)
    print (sum6,"\t",end="",file = f1)
    print (sum7,"\t",end="",file = f1)
    print (sum8,"\t",end="",file = f1)
    print (sum9,"\t",end="",file = f1)
    print (sum10,"\t",end="",file = f1)
    print (sum11,"\t",end="",file = f1)
    print (sum12,"\t",end="",file = f1)
    print (sum13,"\t",end="",file = f1)
    print (sum14,"\t",end="",file = f1)
    print (sum15,file = f1)
   # print ("\n", file = f1)

    sum1=0
    sum2=0
    sum3=0
    sum4=0
    sum5=0
    sum6=0
    sum7=0
    sum8=0
    sum9=0
    sum10=0
    sum11=0
    sum12=0
    sum13=0
    sum14=0
    sum15=0
