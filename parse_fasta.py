from Bio import SeqIO

rec_iter = SeqIO.parse("ABtype_livingASs.fa", "fasta")
outputf = "AS_30bps.fa"


def AS_makefragment(outputf, start, fraglength):
    while True:
        try:
            sat = next(rec_iter)
        except StopIteration:
            break
        else:
            return sat[start:start+fraglength]


with open(outputf, "w") as outfile:
    i = 0
    while True:
        seq = AS_makefragment(outputf, i, 30)
        SeqIO.write(seq, outfile, "fasta")
    i = i + 10
