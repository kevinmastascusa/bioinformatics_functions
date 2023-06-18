library(shiny)

# UI
ui <- fluidPage(
  titlePanel("Bioinformatics/Genomics R Functions"),
  sidebarLayout(
    sidebarPanel(
      textAreaInput("dna_sequence", "DNA Sequence:", value = "ATGCCGTAATGGCCTAAG"),
      actionButton("calculate_gc", "Calculate GC Content"),
      actionButton("count_nucleotides", "Count Nucleotides"),
      actionButton("translate_sequence", "Translate Sequence")
    ),
    mainPanel(
      verbatimTextOutput("gc_content_output"),
      verbatimTextOutput("nucleotide_counts_output"),
      verbatimTextOutput("translated_sequence_output")
    )
  )
)
