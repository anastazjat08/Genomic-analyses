from Bio import SeqIO
from Bio import pairwise2
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import itertools

# This funstion generates dictionary with sequences for alignments
def generate_dict(seqs_file):
    seqs_representatives_dict = {}

    for record in SeqIO.parse(seqs_file, "fasta"):
        seqs_representatives_dict[record.id] = record.seq

    return seqs_representatives_dict

#  This function calculates the percentage identity for all possible pairs of sequences.
def pair_align_all(seqs_dict):

    labels = list(seqs_dict.keys())

    scores = pd.DataFrame(0.0, index=labels, columns=labels)

    for seq1, seq2 in itertools.combinations_with_replacement(labels, 2):
        aln = pairwise2.align.globalxx(seqs_dict[seq1], seqs_dict[seq2], one_alignment_only=True)[0]
        matches = sum(a == b for a, b in zip(aln.seqA, aln.seqB))
        length = max(len(seqs_dict[seq1]), len(seqs_dict[seq2]))
        percent_identity = matches / length * 100
        scores.loc[seq1, seq2] = percent_identity
        scores.loc[seq2, seq1] = percent_identity

    return scores

# This function plots a heatmap of percentage identities
def plot_heatmap(data_scores):
    plt.figure(figsize=(20, 15))
    sns.heatmap(data_scores, fmt=".1f",cmap='viridis',annot=True)

    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.title('Macierz identyczności')
    plt.tight_layout()
    plt.show()
