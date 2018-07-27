from Bio import SeqIO
import re
import csv


frag_iter = SeqIO.parse("SF3_fragment.fa", "fasta")
inputfile = "chr11_peri1.fa"

outputf = "detect_exactmatch.csv"

with open(outputf, "w") as outfile:
    writer = csv.writer(outfile)
    x = 0
    for fragment in frag_iter:
        for HOR in SeqIO.parse(inputfile, "fasta"):
            # scan for IUPAC; re.I makes search case-insensitive
            match_iter = re.finditer(str(fragment.seq), str(HOR.seq))
            for match in match_iter:
                if match:
                    writer.writerow([fragment.description,
                                     HOR.id, match.start(), fragment.seq])
        x = x + 1
