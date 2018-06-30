# disfunction = at least 2 bp substitution
import re

chr1 = 'BABA'
chr2 = 'ABABABAB'
chr3 = 'BABABABAABABABABA'
chr4 = 'BABBABA'
chr5 = 'BAbABAbA'
chr6 = 'BBBABABAABAbABABABA'
chr7 = 'AB'
chr8 = 'BBABABABABABABA'
chr9 = 'AB'
chr10 = 'BBABA'
chr11 = 'AAbBB'
chr12 = 'ABAB'
chr13 = 'BAABABABABA'
chr14 = 'B'
chr15 = 'ABABBABABABBAB'
chr16 = 'BABA'
chr17 = 'BBbAABBbABAABBBA'
chr17b = 'AAABbbbbAAB'
chr18 = 'bBABABABAbA'
chr20 = 'BBAbA'
chrX = 'BAABbBAABBBA'
chrY = 'AAAAAAAAAAAAAA'


print("SF1")
for i in [chr1, chr3, chr5, chr6, chr7, chr10, chr12, chr16]:
    print("chr <-"+str(re.findall('A+', i)))
    print(i.count('A'))

for i in [chr1, chr3, chr5, chr6, chr7, chr10, chr12, chr16]:
    print("chr <-"+str(re.findall('B+', i)))
    print(i.count('B'))

print("SF2")
for i in [chr2, chr4, chr8, chr9, chr13, chr15, chr18, chr20]:
    print("chr <-"+str(re.findall('A+', i)))
    print(i.count('A'))

for i in [chr2, chr4, chr8, chr9, chr13, chr15, chr18, chr20]:
    print("chr <-"+str(re.findall('B+', i)))
    print(i.count('B'))

print("SF3")
for i in [chr11, chr17, chr17b, chrX]:
    print("chr <-"+str(re.findall('A+', i)))
    print(i.count('A'))

for i in [chr11, chr17, chr17b, chrX]:
    print("chr <-"+str(re.findall('B+', i)))
    print(i.count('B'))
