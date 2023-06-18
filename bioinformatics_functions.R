######################################################
#  bioinformatics_functions.R
#
#  Description: This file contains functions for calculating GC content, counting nucleotides, and translating DNA sequences.
#  Author: Kevin Zhou Mastascusa
#  Date: 6.18.2023
#
#  This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with this program.
#  If not, see <http://www.gnu.org/licenses/>.
#
######################################################



# Function: calculate_gc_content
calculate_gc_content <- function(sequence) {
  sequence <- toupper(sequence)
  gc_count <- sum(strsplit(sequence, "")[[1]] %in% c("G", "C"))
  total_count <- nchar(sequence)
  gc_content <- gc_count / total_count * 100
  return(gc_content)
}

# Function: count_nucleotides
count_nucleotides <- function(sequence) {
  sequence <- toupper(sequence)
  nucleotide_counts <- table(strsplit(sequence, "")[[1]])
  return(nucleotide_counts)
}

# Function: translate_sequence
# Function: translate_sequence
translate_sequence <- function(sequence) {
  sequence <- toupper(sequence)
  genetic_code <- c("TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
                    "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
                    "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG",
                    "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
                    "TAT", "TAC", "CAT", "CAC", "CAA", "CAG", "AAT", "AAC",
                    "AAA", "AAG", "GAT", "GAC", "GAA", "GAG", "TGT", "TGC",
                    "TGG", "CGT", "CGC", "CGA", "CGG", "AGT", "AGC", "AGA",
                    "AGG", "GGT", "GGC", "GGA", "GGG",
                    "TAA", "TAG", "TGA")  # Stop codons added
  amino_acids <- c("F", "F", "L", "L", "L", "L", "L", "L",
                   "I", "I", "I", "M", "V", "V", "V", "V",
                   "S", "S", "S", "S", "P", "P", "P", "P",
                   "T", "T", "T", "T", "A", "A", "A", "A",
                   "Y", "Y", "H", "H", "Q", "Q", "N", "N",
                   "K", "K", "D", "D", "E", "E", "C", "C",
                   "W", "R", "R", "R", "R", "S", "S", "R",
                   "R", "G", "G", "G", "G",
                   "*", "*", "*")  # Stop codons represented as "*"
  codons <- strsplit(sequence, "(?<=.{3})", perl = TRUE)[[1]]
  translated_sequence <- unname(sapply(codons, function(codon) {  # Remove names from the vector
    index <- match(codon, genetic_code)
    if (!is.na(index)) {
      amino_acids[index]
    } else {
      "X"
    }
  }))

  # Find the position of the first stop codon
  stop_codon_index <- match("*", translated_sequence)

  # If a stop codon was found, cut the sequence at that point
  if (!is.na(stop_codon_index)) {
    translated_sequence <- translated_sequence[1:(stop_codon_index-1)]
  }

  return(translated_sequence)
}



# Test: calculate_gc_content
test_calculate_gc_content <- function(sequence) {
  gc_content <- calculate_gc_content(sequence)
  expected_gc_content <- sum(strsplit(sequence, "")[[1]] %in% c("G", "C")) / nchar(sequence) * 100
  tolerance <- 0.01  # Adjust the tolerance level as needed

  diff <- abs(gc_content - expected_gc_content)
  print(paste("GC Content:", gc_content))
  print(paste("Expected GC Content:", expected_gc_content))
  print(paste("Difference:", diff))
  print(paste("Tolerance:", tolerance))

  stopifnot(diff <= tolerance)
}

# Test: count_nucleotides
test_count_nucleotides <- function() {
  sequence <- "ATGCCGTAATGGCCTAAG"
  nucleotide_counts <- count_nucleotides(sequence)
  expected_counts <- c(A = 4, T = 5, G = 4, C = 5)
  stopifnot(all(sort(nucleotide_counts) == sort(expected_counts)))
}

# Test: translate_sequence
# Test: translate_sequence
# Test: translate_sequence
# Test: translate_sequence
# Test: translate_sequence
# Test: translate_sequence
test_translate_sequence <- function() {
  sequence <- "ATGCCGTAATGGCCTAAG"
  translated_sequence <- translate_sequence(sequence)
  print(translated_sequence)
  expected_sequence <- c("M", "P")
  stopifnot(identical(translated_sequence, expected_sequence))
}



# Run the test again
test_translate_sequence()

# Run all tests
test_calculate_gc_content("ATGCCGTAATGGCCTAAG")
test_count_nucleotides()
test_translate_sequence()
# TODO: Checklist
# - [ ] Update the `translate_sequence` function to handle stop codons correctly
# - [ ] Update the `test_translate_sequence` function to use the expected sequence with stop codons
# - [ ] Test the `translate_sequence` function to ensure it returns the expected translated sequence
# - [ ] Review and update the comments and documentation for the functions
# - [ ] Handle edge cases and error cases in the functions (e.g., empty string, invalid characters, sequence length not multiple of 3)
# - [ ] Perform additional tests and validation for all functions to ensure correct behavior in various scenarios
# - [ ] Refactor the code or functions for better readability or performance if necessary
# - [ ] Update or add any additional functions or features as needed

# End of TODO checklist
