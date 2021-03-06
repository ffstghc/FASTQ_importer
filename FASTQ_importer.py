###############
# Import Data #
###############
project0_data = open("project0.fq", 'r')
raw_data = project0_data.readlines()  # Read all lines to later split them up

names = [raw_data[i].split("@")[1].strip() for i in range(0, len(raw_data), 4)]  # Sample/sequence names
seqs = [raw_data[i + 1].strip() for i in range(0, len(raw_data), 4)]  # Sequences
qs_enc = [raw_data[i + 3].strip() for i in range(0, len(raw_data), 4)]  # Quality Scores (1-bit encoded)

##########################
# Average Quality Scores #
##########################
qs_dec_combined = []  # All QSs combined to later split them up


def add_qs(pos):
    for ii in qs_enc[pos]:  # Loop over all elements in QS sequence
        qss = ord(ii) - 33  # QS from ASCII code
        qs_dec_combined.append(qss)  # List of all QS


[add_qs(j) for j in range(0, len(qs_enc))]

qs_dec_split = [None] * len(qs_enc)  # List of all QSs
qs_means = []  # List of all mean QSs
total_qss = len(qs_enc)  # Amount of QSs
length_qss = int(len(qs_dec_combined) / len(qs_enc))  # Calculate length of single QS line from combined list

x = 0
while x < total_qss:
    qs_dec_split[x] = qs_dec_combined[x * length_qss:(x * length_qss + length_qss)]  # Split up combined QS list
    qs_means.append((sum(qs_dec_split[x]) / length_qss))  # Calculate mean QS of single sequence
    qs_means[x] = round(qs_means[x], 4)
    x += 1

###############
# %GC Content #
###############
print("Input to calculate fraction in sequence:\n(1) Total content of 'G' and 'C'."
      "\n(2) Specify single base or base sequence.")

mode = int(input())
amount_c = []
amount_g = []
fractions = []

if mode == 1:

    base1 = "C"
    base2 = "G"

    for i in range(0, len(seqs)):
        amount_c.append(seqs[i].count("C"))
        amount_g.append(seqs[i].count("G"))
        fraction = (amount_c[i] + amount_g[i]) / len(seqs[i])  # Total amount of bases divided by sequence length
        fractions.append(fraction)
        fractions[i] = 100 * round(fractions[i], 4)  # Round and convert to percentage
else:

    print("\nInput nucleoside short code (e.g. 'A', 'T', 'C', 'G' or 'ATC', 'TGC', ...):")

    base1 = input()
    amount_B = []

    for i in range(0, len(seqs)):
        amount_B.append(seqs[i].count(base1))
        fraction = (amount_B[i] / len(seqs[i]))  # Total amount of bases divided by sequence length
        fractions.append(fraction)
        fractions[i] = 100 * round(fractions[i], 4)  # Round and convert to percentage

##########################
# Output Data into Table #
##########################
if mode == 1:

    print("Sample name: | Avg. Q-Score: | GC content:")
    [print("{}\t|\t{}\t|\t{} %".format(names[i], qs_means[i], fractions[i])) for i in range(0, total_qss)]

else:

    print("Sample name: | Avg. Q-Score: |", base1, "content:")
    [print("{}\t|\t{}\t|\t{} %".format(names[i], qs_means[i], fractions[i])) for i in range(0, total_qss)]


##########################################
# Translate nucleotides into amino acids #
##########################################
print("\n" "Input sample name to translate to amino acid sequence (e.g. 'prov-0098', 'Sample_0001') or "
      "'skip' to skip:")

ind_seq = input()  # Get input for sample to translate or skip translation

if ind_seq == 'skip':

    print("Conversion of triplets to amino acid  skipped.")

else:
    ind_seq = names.index(ind_seq)  # Get index for sample name

    seq_replaced = seqs[ind_seq].replace("T", "U")  # Replace "T" in all sequences by "U"

    length_seq = len(seqs[0])
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

    translated_single = [code_sun_dict[i] for i in seq_splitted]  # Convert triplet to amino

    translated_single = ''.join(translated_single).split("STOP")  # Split at stop codons

    print("Protein fragments (translated amino acids):")
    [print(translated_single[i]) for i in range(0, len(translated_single))]


##########################
# Convert FASTQ to FASTA #
##########################
print("\nInput 'yes' if you want to output a .FASTA file of all the imported data. Else leave blank:")
x = input()

if x == 'yes':

    textfile = open('Converted.fasta', 'w')  # Create text file

    for j in names:

        for i in range(0, len(names)):
            textfile.write(">"+names[i])  # Add sample name
            textfile.write("\n")
            textfile.write(seqs[i])  # Add sample sequence
            textfile.write("\n")

    textfile.close()

    print("\nFASTA data file exported. Bye!")

else:

    print("\nFinished without exporting FASTA data. Bye!")
