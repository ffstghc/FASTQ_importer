###############
# Import Data #
###############
project0_data = open("project0.fq", 'r')
raw_data = project0_data.readlines()  # Read all lines to later split them up

col_rows = []  # Sample/sequence names
seq = []  # Sequences
qs_enc = []  # Quality Scores (1-bit encoded)

for i in range(0, len(raw_data), 4):
    col_rows.append(raw_data[i])  # Get sample name lines
    seq.append(raw_data[i + 1])  # Get actual raw sequence lines
    qs_enc.append(raw_data[i + 3])  # Get 1-bit encoded quality score (QS) lines

##########################
# Average Quality Scores #
##########################
qs_dec_combined = []  # All QSs combined to later split them up


def add_qs(pos):
    for ii in qs_enc[pos].strip():  # Loop over all elements in QS sequence
        qss = ord(ii) - 33  # QS from ASCII code
        qs_dec_combined.append(qss)


for j in range(0, len(qs_enc)):
    add_qs(j)

qs_dec_split = [None] * len(qs_enc)  # List of all QSs
qs_means = []  # List of all mean QSs
total_qss = len(qs_enc)  # Amount of QSs
length_qss = int(len(qs_dec_combined) / len(qs_enc))
x = 0
while x < total_qss:
    qs_dec_split[x] = qs_dec_combined[x * length_qss:(x * length_qss + length_qss)]
    qs_means.append(sum(qs_dec_split[x]) / total_qss)
    x += 1

###############
# %GC Content #
###############
print("Input to calculate fraction in sequence:\n(1) Total content of 'G' and 'C'"
      "\n(2) Specify single base or base sequence")
mode = int(input())
amount_c = []
amount_g = []
fractions = []

if mode == 1:
    base1 = "C"
    base2 = "G"
    for i in range(0, len(seq)):
        amount_c.append(seq[i].count("C"))
        amount_g.append(seq[i].count("G"))
        fraction = (amount_c[i] + amount_g[i]) / len(seq[i])  # Total amount of bases divided by sequence length
        fractions.append(fraction)
else:
    print("\nInput nucleoside short code (e.g. 'A', 'T', 'C', 'G' or 'ATC', 'TGC', ...):")
    base1 = input()
    amount_B = []
    for i in range(0, len(seq)):
        amount_B.append(seq[i].count(base1))
        fraction = (amount_B[i] / len(seq[i]))  # Total amount of bases divided by sequence length
        fractions.append(fraction)

###########################
# Combine Data into Table #
###########################
col_names = "Sequence name", "%GC content", "Avg. QS"  # Columns for sequence specific table

for i in range(0, total_qss):
    if mode == 1:
        print("Sample name: {}Average QS: {} | 'GC' content: {}".format(col_rows[i], qs_means[i], fractions[i]))
    else:
        print("Sample name: {}Average QS: {} | Base/sequence content: {}".format(col_rows[i], qs_means[i], fractions[i]))

##########################################
# Translate nucleotides into amino acids #
##########################################
seq_replaced = [None] * len(seq)

for i in range(0, len(seq)):
    seq_replaced[i] = seq[i].replace("T", "U")  # Replace "T" in all sequences by "U"

length_seq = len(seq[0])
total_acids = int(length_seq / 3)
seq_splitted = [None] * total_acids
y = 0
n_seqs = len(seq)

for i in range(0, n_seqs):
    while y < total_acids:
        new = seq_replaced[i][y * 3:(y * 3 + 3)]
        seq_splitted[y] = new
        y += 1

code_sun_dict = {
    'UUU': 'Phe',
    'UUC': 'Phe',
    'UUA': 'Leu',
    'UUG': 'Leu',
    'CUU': 'Leu',
    'CUC': 'Leu',
    'CUA': 'Leu',
    'CUG': 'Leu',
    'AUU': 'Ile',
    'AUC': 'Ile',
    'AUA': 'Ile',
    'AUG': 'Met',
    'GUU': 'Val',
    'GUC': 'Val',
    'GUA': 'Val',
    'GUG': 'Val',
    'UCU': 'Ser',
    'UCC': 'Ser',
    'UCA': 'Ser',
    'UCG': 'Ser',
    'CCU': 'Pro',
    'CCC': 'Pro',
    'CCA': 'Pro',
    'CCG': 'Pro',
    'ACU': 'Thr',
    'ACC': 'Thr',
    'ACA': 'Thr',
    'ACG': 'Thr',
    'GCU': 'Ala',
    'GCC': 'Ala',
    'GCA': 'Ala',
    'GCG': 'Ala',
    'UAU': 'Tyr',
    'UAC': 'Tyr',
    'UAA': 'STOP',
    'UAG': 'STOP',
    'CAU': 'His',
    'CAC': 'His',
    'CAA': 'Gln',
    'CAG': 'Gln',
    'AAU': 'Asn',
    'AAC': 'Asn',
    'AAA': 'Lys',
    'AAG': 'Lys',
    'GAU': 'Asp',
    'GAC': 'Asp',
    'GAA': 'Glu',
    'GAG': 'Glu',
    'UGU': 'Cys',
    'UGC': 'Cys',
    'UGA': 'STOP',
    'UGG': 'Trp',
    'CGU': 'Arg',
    'CGC': 'Arg',
    'CGA': 'Arg',
    'CGG': 'Arg',
    'AGU': 'Ser',
    'AGC': 'Ser',
    'AGA': 'Arg',
    'AGG': 'Arg',
    'GGU': 'Gly',
    'GGC': 'Gly',
    'GGA': 'Gly',
    'GGG': 'Gly'
}

translated_single = []
for i in seq_splitted:
    code = code_sun_dict[i]
    translated_single.append(code)

##########################
# Convert FASTQ to FASTA #
##########################
print("\nInput 'yes' if you want to output a .fasta file of the imported data. Else leave blank:")
x = input()
if x == 'yes':
    textfile = open('Converted.fasta', 'w')

    for i in range(0, len(col_rows)):
        textfile.write(col_rows[i])
        textfile.write(seq[i])
    textfile.close()
    print("\nFile exported")
else:
    print("\nFinished without exporting FASTA data. Bye!")

""""
# Triplet to amino acid over all sequences (NOT WORKING YET)
seqs_translated = []

def translate_seq(pos):
    for iii in seq_replaced[pos].strip():
        seqs_translated.append(code_sun_dict[iii])


for j in range(0, len(seq)):
    translate_seq(j)
"""
