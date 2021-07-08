import collections
from structures import *
from collections import Counter
# DNA Toolkit File
# LIST of 4 Nucleotides
Nucleotides = ["A", "C", "G", "T"]
DNA_ReverseComplement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

# Validate Sequence Function; Check to make sure it is a valid DNA String
# Making sure it does not contain any other letters
# Accepts a string, creates a tmpseq which is uppercase, loops through each character to validate it is a nucleotide
# Need an dna_seq.upper() since a != A, etc.
def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

    # Counts Nucleotide function
    # Pass it a string validated beforehand
    # Define a temporary dictionary; a fast hash map
    # Loop through string passed and if it's A, C, G, or T, it will increase the corresponding value

def countNucFrequency(seq):
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict
    #return dict(collections.Counter(seq))
    # for optimization, use a counter ^

def transcription(seq):
    """DNA -> RNA Transcription. Replacing Thymine with Uracil"""
    return seq.replace("T", "U")

def reverse_complement(seq):
    """Swapping Adenine with Thymine and Guanine with Cytosine. Reversing newly generated string"""
    return ''.join([DNA_ReverseComplement[nuc] for nuc in seq][::-1])
    # Pythonic approach. Faster solution
    # mapping = str.maketrans('ATCG', 'TAGC')
    # return seq.translate(mapping)[::-1]
    # [::-1] gives inverse output; TEST == TSET

def gc_content(seq):
    """GC Content in a DNA/RNA Sequence"""
    return round(((seq.count('C') + seq.count('G')) / len(seq) * 100), 6)

def gc_content_subsec(seq, k=20):
    """GC Content in a DNA/RNA sub-sequence length k. k=20 by default"""
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        res.append(gc_content(subseq))
    return res
    # sequence and window of 20
    # defines empty list
    # use a for loop, start at 0, and up to length of sequence minus the k
    # jumps will be the size of the window
    # calculates the subseq by calling gc_content
    # list is populated with values, then returned
    # the higher the GC content, the stronger DNA will be when binding with the primer or the other complementary strands
    # a lower GC content ~40% is desirable
    # GC content determines melting temperature (temperature of separation in single strands)
    # important for determining conditions for efficient PCR reactions

def translate_seq(seq, init_pos=0):
    """Translates a DNA sequence into an aminoacid sequence"""
    return [DNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]
    # accepts a sequence and an initialization position of 0
    # uses IP for generation of reading frames; generate 6 of them
    # uses a dictionary to find corresponding values from the keys provided
    # uses jumps of 3
    # example: ASDFGASKJSDF sifts through ASD, FGA, SKJ, SDF jumps 3, then reads 3, matches with key in DNA Codon table
    # by having init_pos = 0, allows for it to start reading at position 0 unless provided otherwise

def codon_usage(seq, aminoacid):
    """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
    tmpList = []
    for i in range(0, len(seq) - 2,3):
        if DNA_Codons[seq[i:i + 3]] == aminoacid:
            tmpList.append(seq[i:i + 3])
    # instead of returning every single amino acid
    # it checks if it's an aminoacid we want to build a statistic on; if so, appends it to the lsit
    freqDict = dict(Counter(tmpList))
    # use counter to create a dictionary, setting a key as a value, value of key = how many times it appears in the list
    totalWight = sum(freqDict.values())
    # summing up all values in dictionary
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalWight, 2)
        # calculates frequency (percentage)
    return freqDict

def gen_reading_frames(seq):
    """Generate the six reading frames of a DNA sequence, including the reverse complement"""
    frames = []
    # [] creates an empty list
    frames.append(translate_seq(seq, 0))
    frames.append(translate_seq(seq, 1))
    frames.append(translate_seq(seq, 2))
    frames.append(translate_seq(reverse_complement(seq), 0))
    frames.append(translate_seq(reverse_complement(seq), 1))
    frames.append(translate_seq(reverse_complement(seq), 2))
    # append a list of amino acid
    # 0, 1, 2 tells it where to start reading from
    return frames
    # return frames returns a list that contains the lists that we appended to it

def proteins_from_rf(aa_seq):
    """Computes all possible proteins in an aminoacid sequence and returns a list of possible proteins"""
    current_prot = []
    proteins = []
    # accumulates protein in the list, proteins = [] basically is the sifter
    for aa in aa_seq:
        if aa == "_":
            # STOP accumulating amino acids if _ - STOP Codon was found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            # START accumulating amino aacids if M - START Codon was found
            if aa == "M":
                current_prot.append("")
                # makes sure length isn't 0, so it can start accumulating values
            for i in range(len(current_prot)):
                # current_prot measures the amount of strings not characters
                current_prot[i] += aa
    return proteins

def all_proteins_from_orfs(seq, startReadPos = 0, endReadPos = 0, ordered=False):
    """Compure all possible proteins for all open reading frames"""
    """Proein Seaarch Database: https://www.ncbi.nlm.gov/nuccore/NM_001185097.2"""
    """API can be used to pull protein information"""
    # Generate Reading Frames; adding one check; acacepts a sequence, start and end reading position
    # Specified a short parameter
    # ordered=False; be able to order an unordered list
    if endReadPos > startReadPos:
        rfs = gen_reading_frames(seq[startRead: endRead])
    else:
        rfs = gen_reading_frames(seq)
            # determines whether we want to read the whole list or a specific ranged
    res = []
    # creates an empty list for all the proteins we find
    for rf in rfs:
        prots = proteins_from_rf(rf)
        # each reading frame generates a protein from each rf calling the function proteins_from_rf
        for p in prots:
            res.append(p)

    if ordered:
        return sorted(res, key=len, reverse=True)
    return res
    # do we want our output: if yes; use sorted (built in Python function)
    # sort by a length of entries
    # reverse key makes sure list sorted is from longest to shortest string