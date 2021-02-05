from collections import deque
import requests
class Protein():
    """"Class that contains different proteins

    Notes:
        Class Protein requires the modules collections and requests

    Attributes:
        lookup (dict): nested dict that contains property and data pairs for all amino acids, e.g. "hydropathy":{"A":1.6, "K":3.2, ...}
        sequence (str): long string that contains the amino acid sequence of the protein
        property (key): tracks which property was chosen from the lookup
    """
    def __init__(self, lookup):
        self.lookup = lookup
        self.sequence = None
        self.property = None

    def get_data(self, identifier, protein_name):
        """"Gets protein sequence from uniprot
        Function takes the protein name and identifier to download the fasta file from uniprot and extracts the amino acid sequence from it

        Args:
            identifier (str): specific identifier for the protein, e.g. P32248
            protein_name (str): name of the protein for to name the downloaded fasta file

        Yields:
            str: sequence of the protein, downloaded from uniprot, sets to the attribute self.sequence
        """

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
        """"Matches the protein's sequence to a given property
        Based on self.lookup, the given property and the size of the sliding window an averaged list of values is returned

        Args:
            property (Key): selects a specific dict of self.lookup to be used for mapping amino acids to values
            sliding_window_size (int): size of the sliding window the values are averaged over

        Returns:
            averaged_list (list): list of the same length as self.sequence with the mapped values from self.lookup[self.property]
        """
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

