import plotly.graph_objects as go
import pandas as pd
from collections import deque
import numpy as np

aa_df = pd.read_csv("../data/amino_acid_properties.csv")

aa_hydropathy = dict()
for pos, aa in enumerate(aa_df["1-letter code"]):
    aa_hydropathy[aa] = aa_df["hydropathy index (Kyte-Doolittle method)"][pos]

with open("../data/P32249.fasta", "r") as fasta:
    sequence_reader = fasta.readlines()

to_delete = []
for pos, row in enumerate(sequence_reader):
    if row[0] == ">":
        to_delete.append(pos)

for pos in sorted(to_delete, reverse=True):
    sequence_reader.pop(pos)

sequence = "".join(sequence_reader)
sequence = sequence.replace("\n", "")

def get_hydropathy_list(sequence, aa_hydropathy, sliding_window_size = 1):
    hydropathy_list = []
    for aa in sequence:
        hydropathy_list.append(aa_hydropathy[aa])

    averaged_list = []
    window = deque([], maxlen=sliding_window_size)
    for pos, hydropathy in enumerate(hydropathy_list):
        average = 0
        window.append(hydropathy)
        for neighbor in window:
            average += neighbor
        average = average/sliding_window_size
        averaged_list.append(average)
    return averaged_list



data = [
    go.Scatter(
    x=[i for i in range(len(sequence))],
    y=get_hydropathy_list(sequence, aa_hydropathy, 1),
    name="No average"),
    go.Scatter(
    x=[i for i in range(len(sequence))],
    y=get_hydropathy_list(sequence, aa_hydropathy, 5),
    name="Over 5"),
    go.Scatter(
    x=[i for i in range(len(sequence))],
    y=get_hydropathy_list(sequence, aa_hydropathy, 10),
    name="Over 10"),
    go.Scatter(
    x=[i for i in range(len(sequence))],
    y=get_hydropathy_list(sequence, aa_hydropathy, 20),
    name="Over 20")
]
fig = go.Figure(data=data)
fig.layout.title.text = "Amino acid hydrophathy"
fig.layout.title.font.size = 40
fig.layout.xaxis.title.text = "Amino acid in GPCR"
fig.layout.yaxis.title.text = "Hydropathy"
fig.show()
