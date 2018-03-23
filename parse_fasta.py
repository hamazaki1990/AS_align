from Bio import SeqIO


def AS_makefragment(inputf, start, fraglength):
    rec_iter = SeqIO.parse(inputf, "fasta")
    while True:
        try:
            sat = next(rec_iter)
        except StopIteration:
            break
        else:
            yield sat[start:start+fraglength]


inputf = "ABtype_livingASs_gapremoved.fa"

i = 0
while True:
    if i > 171:
        break
    else:
        outputf = "AS_30bps"+str(i)+".fa"
        with open(outputf, "w") as outfile:
            seq = AS_makefragment(inputf, i, 30)
            SeqIO.write(seq, outfile, "fasta")
        i = i + 10
