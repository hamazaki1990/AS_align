import random
import csv


def get_matchlength(l, k):
    mutations = random.sample(range(l), k)
    mutations.sort()
    matchlength = [mutations[i+1] - mutations[i] for i in range(k-1)]
    matchlength.extend([mutations[0], l-mutations[k-1]])
    return matchlength


outputf = "match_distribution.csv"

with open(outputf, "w") as outfile:
    writer = csv.writer(outfile)
    for i in range(10000):
        a = get_matchlength(342, 50)
        writer.writerow(a)
