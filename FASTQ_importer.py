###############
# Import Data #
###############
project0_data = open("project0.fq", 'r')
raw_data = project0_data.readlines()  # Read all lines to later split them up

col_rows = []  # Sample/sequence names
seq = []  # Sequences
qs_enc = []  # Quality Scores (1-bit encoded)

for i in range(0, len(raw_data), 4):
    col_rows.append(raw_data[i].split("@")[1].strip())  # Get sample names
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
length_qss = int(len(qs_dec_combined) / len(qs_enc))  # Calculate length of single QS line from combined list

x = 0
while x < total_qss:
    qs_dec_split[x] = qs_dec_combined[x * length_qss:(x * length_qss + length_qss)]  # Split up combined QS list
    qs_means.append((sum(qs_dec_split[x]) / total_qss))  # Calculate mean QS of single sequence
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
        fractions[i] = 100 * round(fractions[i], 4)  # Round and convert to percentage
else:
    print("\nInput nucleoside short code (e.g. 'A', 'T', 'C', 'G' or 'ATC', 'TGC', ...):")
    base1 = input()
    amount_B = []
    for i in range(0, len(seq)):
        amount_B.append(seq[i].count(base1))
        fraction = (amount_B[i] / len(seq[i]))  # Total amount of bases divided by sequence length
        fractions.append(fraction)
        fractions[i] = 100 * round(fractions[i], 4)  # Round and convert to percentage

##########################
# Output Data into Table #
##########################
if mode == 1:
    print("Sample name: | Avg. Q-Score: | GC content:")
    for i in range(0, total_qss):
        print("{}   |   {}  |   {} %".format(col_rows[i], qs_means[i], fractions[i]))
else:
    print("Sample name: | Avg. Q-Score: |", base1, "content:")
    for i in range(0, total_qss):
        print("{}   |   {}  |   {} %".format(col_rows[i], qs_means[i], fractions[i]))

##########################################
# Translate nucleotides into amino acids #
##########################################
print("\n" "Input sample name to translate to amino acid sequence (e.g. 'prov-0098', 'Sample_0001') or "
      "'skip' to skip:")

ind_seq = input()  # Get input for sample to translate or skip translation

if ind_seq == 'skip':
    print("Sample input for triplet to amino acid conversion skipped")
else:
    ind_seq = col_rows.index(ind_seq)  # Get index for sample name

    seq_replaced = seq[ind_seq].replace("T", "U")  # Replace "T" in all sequences by "U"

    length_seq = len(seq[0])
    total_acids = int(length_seq / 3)
    seq_splitted = [None] * total_acids
    y = 0
    while y < total_acids:
        new = seq_replaced[y * 3:(y * 3 + 3)]
        seq_splitted[y] = new
        y += 1

    # Three letter nucleotide to amino acid coding
    code_sun_dict = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met', 'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr', 'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP', 'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys', 'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp', 'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg', 'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }

    translated_single = []
    for i in seq_splitted:
        code = code_sun_dict[i]  # Convert triplet to amino acid
        translated_single.append(code)

    translated_single = ''.join(translated_single).split("STOP")  # Split at stop codons

print("Protein fragments (amino acids):")
for i in range(0, len(translated_single)):
    print(translated_single[i])

##########################
# Convert FASTQ to FASTA #
##########################
print("\nInput 'yes' if you want to output a .FASTA file of all the imported data. Else leave blank:")
x = input()
if x == 'yes':
    textfile = open('Converted.fasta', 'w')  # Create text file

    for i in range(0, len(col_rows)):
        textfile.write(col_rows[i])  # Add sample name
        textfile.write("\n")
        textfile.write(seq[i])  # Add sample sequence
    textfile.close()
    print("\nFASTA data file exported. Bye!")
else:
    print("\nFinished without exporting FASTA data. Bye!")
