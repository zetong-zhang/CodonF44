# Author: Zhang Zetong
# Date: 2024-08-17
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import CodonAdaptationIndex

data = pd.read_csv(filepath_or_buffer="F44.CDS.tsv", sep='\t')
record = SeqIO.read(handle='F44.Genome.gb', format='gb')
seqs = []

for i in data.index:
    locus = data['gene_id'][i]

    for feature in record.features:
        if feature.type == 'CDS' and feature.qualifiers['locus_tag'][0] == locus:
            start, end = feature.location.start, feature.location.end
            if (end - start) % 3 == 0:
                seq = record.seq[start:end]
                if feature.strand == -1:
                    seq = seq.reverse_complement()
                seq_record = SeqRecord(seq, id=locus, description='')
                seqs.append(seq_record)
                break

print("Numbers: ", len(seqs))
SeqIO.write(seqs, handle='F44.CDS.fa', format='fasta')

seqs = [r.seq for r in seqs]

codons = {
    "C": ["TGT", "TGC"],
    "D": ["GAT", "GAC"],
    "S": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
    "Q": ["CAA", "CAG"],
    "M": ["ATG"],
    "N": ["AAC", "AAT"],
    "P": ["CCT", "CCG", "CCA", "CCC"],
    "K": ["AAG", "AAA"],
    "*": ["TAG", "TGA", "TAA"],
    "T": ["ACC", "ACA", "ACG", "ACT"],
    "F": ["TTT", "TTC"],
    "A": ["GCA", "GCC", "GCG", "GCT"],
    "G": ["GGT", "GGG", "GGA", "GGC"],
    "I": ["ATC", "ATA", "ATT"],
    "L": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
    "H": ["CAT", "CAC"],
    "R": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
    "W": ["TGG"],
    "V": ["GTA", "GTC", "GTG", "GTT"],
    "E": ["GAG", "GAA"],
    "Y": ["TAT", "TAC"]
}

cai = CodonAdaptationIndex(seqs)

print('Codon CODONS[21][6] = {')
for key, value in codons.items():
    s, i = sum([cai[item] for item in value]), 0
    d = {item: round(cai[item], 2) for item in value}
    d_order = sorted(d.items(), key=lambda x: x[1], reverse=True)

    print("/* %s */\t{" % key, end='')

    for item in d_order:
        i += item[1] / s
        print('{"%s", %.2f, %.2f}' % (item[0], item[1], round(i, 2)), 
              end=', ' if d_order[len(d_order) - 1] != item else '')

    print("}, ")
print("};")
