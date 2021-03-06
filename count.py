# intact literal A/B
import re

chr1 = 'BaBa'
chr2 = 'ABABABaB'
chr3 = 'BABaBaBaaBabaBaBa'
chr4 = 'bABBaBa'
chr5 = 'BabaBaBa'
chr6 = 'bBbaBabaaBabaBaBaBa'
chr7 = 'aB'
chr8 = 'BBAbABaBABABABA'
chr9 = 'AB'
chr10 = 'BBaBa'
chr11 = 'aabbb'
chr12 = 'aBaB'
chr13 = 'BaABaBABaBA'
chr14 = 'B'
chr15 = 'ABaBBABaBABBAB'
chr16 = 'BaBa'
chr17 = 'bbbAabbbabAabbba'
chr17b = 'aAabbbbbAab'
chr18 = 'bBaBaBABAba'
chr20 = 'bBAbA'
chrX = 'BAaBbbAaBBba'
chrY = 'AAaAaAAaAAaAAa'


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
