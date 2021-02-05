import plotly.graph_objects as go
import pandas as pd

aa_df = pd.read_csv("../data/amino_acid_properties.csv")

data = [
    go.Bar(
        x=aa_df["1-letter code"],
        y=aa_df["hydropathy index (Kyte-Doolittle method)"],

        text="Name:" + aa_df['Name'] + "<br />" + \
             "Weight: " + aa_df['Residue Weight'].astype(str) + "<br />" + \
             "Formula:" + aa_df['Residue Formula']

    )
]

fig = go.Figure(data=data)
fig.update_layout(template="plotly_dark", title="AA hydropathy index")

fig.show()