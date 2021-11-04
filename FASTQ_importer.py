###############
# Import Data #
###############
project0_data = open("project0.fq", 'r')
raw_data = project0_data.readlines()  # Read all lines to later split them up

col_rows = []  # Sample/sequence names
seq = []  # Sequences
qc_enc = []  # QC Scores (1-bit encoded)

for i in range(0, len(raw_data), 4):
    col_rows.append(raw_data[i])  # Get sample name lines
    seq.append(raw_data[i + 1])  # Get actual raw sequence lines
    qc_enc.append(raw_data[i + 3])  # Get 1-bit encoded quality score lines

#####################
# Create Dictionary #
#####################
keys = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',',
        '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8',
        '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D',
        'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
        'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
qc_score_range = list(range(0, 58))
qc_dict = dict(zip(keys, qc_score_range))  # Create dictionary for QC Score decoding

##########################
# Average Quality Scores #
##########################
qc_dec_combined = []  # All qcs combined to later split them up


def add_qc(pos):
    for ii in qc_enc[pos].strip():  # Loop over all elements in first ("[0]") QC sequence
        qc_dec_combined.append(qc_dict[ii])


for j in range(0, len(qc_enc)):
    add_qc(j)

qc_dec_split = [None] * len(qc_enc)  # List of all QCs
qc_means = []
total_qcs = len(qc_enc)
length_qcs = int(len(qc_dec_combined) / len(qc_enc))
x = 0
while x < total_qcs:
    qc_dec_split[x] = qc_dec_combined[x * length_qcs:(x * length_qcs + length_qcs)]
    qc_means.append(sum(qc_dec_split[x]) / total_qcs)
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

for i in range(0, total_qcs):
    if mode == 1:
        print("Sample name: {} | Average QS: {} | 'GC' content: {}".format(col_rows[i], qc_means[i], fractions[i]))
    else:
        print("Sample name: {} | Average QS: {} | Base/sequence content: {}".format(col_rows[i], qc_means[i],
                                                                                    fractions[i]))
