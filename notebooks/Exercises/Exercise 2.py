from collections import Counter
import csv
import matplotlib.pyplot as plt

human_FASTA = "uniprot-filtered-organism Homo+sapiens+(Human)+[9606] +AND+review--.fasta"
plant_FASTA = "uniprot-filtered-organism__Arabidopsis+thaliana+(Mouse-ear+cress)+[370--.fasta"
yeast_FASTA = "uniprot-filtered-organism Saccharomyces+cerevisiae+(strain+ATCC+204508+%2--.fasta"

def count_AA(file, organism = ""):
    # load the fasta file
    with open(file, "r") as fasta:
        sequence_reader = fasta.readlines()
    fasta.close()

    # remove rows starting with ">"
    to_delete = []
    for pos, row in enumerate(sequence_reader):
        if row[0] == ">":
            to_delete.append(pos)

    for pos in sorted(to_delete, reverse=True):
        sequence_reader.pop(pos)

    sequence = "".join(sequence_reader)
    output = Counter(sequence)

    # save in a csv file
    with open("AA-Count.csv", "w", newline="") as AA_Count:
        csv_writer = csv.writer(AA_Count)
        for aa in output.items():
            csv_writer.writerow(aa)

    aa = []
    freq = []
    for pair in output.most_common():
        aa.append(pair[0])
        freq.append(pair[1])

    # plotting
    plt.bar(aa, freq, width = .5, color = "red")
    plt.ylabel("Frequency")
    plt.xlabel("Amino acid")
    plt.title(organism)
    plt.show()

count_AA(human_FASTA, "Human")
count_AA(plant_FASTA, "Plant")
count_AA(yeast_FASTA, "Yeast")