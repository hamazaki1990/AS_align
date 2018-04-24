import re

chr1 = 'BABa'
chr2 = 'ABaBABaB'
chr3 = 'BABABaBaaBabABaBA'
chr4 = 'bABBABa'
chr5 = 'BAbABAbA'
chr6 = 'BBbaBabaaBabaBaBABA'
chr7 = 'aB'
chr8 = 'BBAbABaBABABABa'
chr9 = 'AB'
chr10 = 'BBABa'
chr11 = 'Aabbb'
chr12 = 'aBAB'
chr13 = 'BaaBaBABaBa'
chr14 = 'B'
chr15 = 'aBaBBABaBaBBAB'
chr16 = 'BABA'
chr17 = 'BBbAaBBbaBaaBBBa'
chr17b = 'aaaBbbbbaaB'
chr18 = 'bBaBABABAba'
chr20 = 'BbAbA'
chrX = 'BAaBbbAaBBba'
chrY = 'AaaaaAAaaAaaAA'


print("pJa conserved")
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
