import random
from DNAToolkit import *
from utilities import colored

# Creating a random DNA Sequence for Testing:
# Choice method from random module
# 'Here's a list of nucleotides, give list to function; picks random character from list 50 times'
# Outputs random DNA String from 50 characters; change range to change string length

randDNAStr = ''.join([random.choice(Nucleotides)
                     for nuc in range(50)])

DNAStr = validateSeq(randDNAStr)

print(f'\nSequence: {DNAStr}\n')
print(f'[1] + Sequence Length: {len(DNAStr)}\n')
print(colored(f'[2] + Nucleotide Frequency: {countNucFrequency(DNAStr)}\n'))

print(f'[3] + DNA/RNA Transcription: {transcription(DNAStr)}\n')

print(f"[4] + DNA String + Reverse Complement:\n5' {DNAStr} 3'")
print(f"   {''.join(['|' for c in range(len(DNAStr))])}")
print(f"3' {(reverse_complement(DNAStr))[::-1]} 5' [Complement]")
print(f"5' {reverse_complement(DNAStr)} 3' [Reverse Complement]\n")

print(f'[5] + GC Content: {gc_content(DNAStr)}%\n')
print(f'[6] + GC Content in Subsection k=5: {gc_content_subsec(DNAStr, k=5)}\n')
print(f'[7] + Aminoacids Sequence from DNA: {translate_seq(DNAStr, 0)}\n')
print(f'[8] + Codon Frequency (L): {codon_usage(DNAStr, "L")}\n')
print(f'[9] + Reading_frames:')
for frames in gen_reading_frames(DNAStr):
    print(frames)
print('\n[10] + All proteins in 6 open reading frames:')
for prot in all_proteins_from_orfs(NM_000208_4, 0, 0, True):
    print(f'{prot}')
