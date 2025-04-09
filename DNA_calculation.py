M_dAMP = 331.2
M_dTMP = 322.2
M_dGMP = 347.2
M_dCMP = 307.2
M_H = 1.0

def complement(seq: str) -> str:
    d = {'A': 'T',
         'C': 'G', 
         'T': 'A', 
         'G': 'C'}
    comp_seq: str = ''
    for nuc in seq:
        comp_seq += d[nuc]
    return comp_seq

def calculate_weight(seq: str) -> float:
    d = {'A': M_dAMP, 'T': M_dTMP, 'G': M_dGMP, 'C': M_dCMP}
    molmass = sum(d[nuc] for nuc in seq) - (len(seq)-1) * M_H
    return molmass

if __name__ == '__main__':
    # For every base one nucleotide base reaction, one H disappears
    sequence: str = ''
    with open('data/standard_fragment_sequence.txt', 'r') as f:
        for line in f:
            sequence += line.strip()

    complementary_sequence: str = complement(sequence)

    M_seq = calculate_weight(sequence)
    M_cseq = calculate_weight(complementary_sequence)
    M_dsDNA = M_seq + M_cseq

    lines = [f'M(seq) = {M_seq} Da\n', f'M(cseq) = {M_cseq} Da\n', f'M(dsDNA) = {M_dsDNA} Da\n']
    for line in lines:
        print(line, end='')

    with open('data/dsDNA_molmass.txt', 'w') as f:
        f.writelines(lines)
