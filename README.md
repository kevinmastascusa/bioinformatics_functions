# Bioinformatics/Genomics R Functions

This project provides custom R functions designed for tasks commonly performed in bioinformatics and genomics. The functions allow you to calculate GC content, count nucleotides, and translate DNA sequences.

## Functions

### calculate_gc_content

**Description:** Calculates the GC content of a given DNA sequence.

**Input:** 
- `sequence`: A character string representing the DNA sequence.

**Output:** 
- `gc_content`: The GC content of the input sequence as a percentage.

### count_nucleotides

**Description:** Counts the occurrence of each nucleotide (A, T, G, C) in a given DNA sequence.

**Input:** 
- `sequence`: A character string representing the DNA sequence.

**Output:** 
- `nucleotide_counts`: A named vector containing the counts of each nucleotide.

### translate_sequence

**Description:** Translates a given DNA sequence into its corresponding amino acid sequence based on the genetic code.

**Input:** 
- `sequence`: A character string representing the DNA sequence.

**Output:** 
- `translated_sequence`: A character vector representing the translated amino acid sequence.

## Usage

```R
# Include the function definitions
source("bioinformatics_functions.R")

# Calculate the GC content of a DNA sequence
dna_sequence <- "ATGCCGTAATGGCCTAAG"
gc_content <- calculate_gc_content(dna_sequence)
print(gc_content)  # Output: 50

# Count the nucleotides in a DNA sequence
nucleotide_counts <- count_nucleotides(dna_sequence)
print(nucleotide_counts)
# Output: A  C  G  T 
#         5  4  5  4 

# Translate a DNA sequence to amino acids
translated_sequence <- translate_sequence(dna_sequence)
print(translated_sequence)
# Output: [1] "M" "P"
