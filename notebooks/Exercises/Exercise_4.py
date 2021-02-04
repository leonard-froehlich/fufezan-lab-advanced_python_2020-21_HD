import plotly.graph_objects as go
import pandas as pd
from collections import deque
import requests


class Protein():
    def __init__(self, lookup):
        self.lookup = lookup
        self.sequence = None
        self.property = None

    def get_data(self, identifier, protein_name):

        url = 'https://www.uniprot.org/uniprot/' + identifier + '.fasta?fil=reviewed:yes'
        r = requests.get(url)
        sequence = ""
        seq_list = []
        protein_fasta = '../../data/' + protein_name + '.fasta'

        with open(protein_fasta, "wb") as fasta:
            fasta.write(r.content)
            fasta.close()

        with open(protein_fasta) as fasta:
            sequence_reader = fasta.readlines()

        to_delete = []
        for pos, row in enumerate(sequence_reader):
            if row[0] == ">":
                to_delete.append(pos)

        for pos in sorted(to_delete, reverse=True):
            sequence_reader.pop(pos)

        sequence = "".join(sequence_reader)
        sequence = sequence.replace("\n", "")

        self.sequence = sequence

    def map(self, property, sliding_window_size):
        self.property = property
        value_list = []
        for aa in self.sequence:
            value_list.append(self.lookup[property][aa])

        averaged_list = []
        window = deque([], sliding_window_size)
        for value in value_list:
            window.append(value)
            average = sum(window)/len(window)
            averaged_list.append(average)
        return averaged_list




aa_df = pd.read_csv("../../data/amino_acid_properties.csv")
aa_hydropathy = dict()
aa_pI = dict()
aa_Surface = dict()

for pos, aa in enumerate(aa_df["1-letter code"]):
    aa_hydropathy[aa] = aa_df["hydropathy index (Kyte-Doolittle method)"][pos]
    aa_pI[aa] = aa_df["pI"][pos]
    aa_Surface = aa_df["Accessible surface"][pos]
lookup = {"hydropathy":aa_hydropathy, "pI":aa_pI}


GPCR = Protein(lookup)
GPCR.get_data("P32248", "GPCR")

data = [
    go.Scatter(
    x=[i for i in range(len(GPCR.sequence))],
    y=GPCR.map("hydropathy", 1),
    name="No average"),
    go.Scatter(
    x=[i for i in range(len(GPCR.sequence))],
    y=GPCR.map("hydropathy", 5),
    name="Average over 5 amino acids"),
    go.Scatter(
    x=[i for i in range(len(GPCR.sequence))],
    y=GPCR.map("hydropathy", 10),
    name="Average over 10 amino acids"),
    go.Scatter(
    x=[i for i in range(len(GPCR.sequence))],
    y=GPCR.map("hydropathy", 20),
    name="Average over 20 amino acids")
]

fig = go.Figure(data=data)
fig.layout.title.text = "Properties of GPCR 183"
fig.layout.title.font.size = 40
fig.layout.xaxis.title.text = "Amino acid position"
fig.layout.yaxis.title.text = f"{GPCR.property}"
fig.layout.yaxis.title.font.size = 20
fig.add_shape(type="line", x0=0, x1=len(GPCR.sequence), y0=0, y1=0)
fig.show()

