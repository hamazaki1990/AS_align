from Bio import SeqIO


def AS_makefragment(outputf, start, fraglength):
    with open(outputf, "w") as outfile:
        while True:
            try:
                sat = next(rec_iter)
            except StopIteration:
                break
            else:
                seq = sat[start:start+fraglength]
                SeqIO.write(seq, outfile, "fasta")


i = 0
while True:
    if i > 150:
        break
    else:
        rec_iter = SeqIO.parse("Btype_livingASs.fa", "fasta")
        outputf = "Btype_"+str(i)+".fa"
        AS_makefragment(outputf, i, 55)
    i = i + 45
