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

# Server
server <- function(input, output) {

  # Calculate GC content
  observeEvent(input$calculate_gc, {
    dna_sequence <- input$dna_sequence
    gc_content <- calculate_gc_content(dna_sequence)
    output$gc_content_output <- renderPrint({
      paste("GC Content:", gc_content)
    })
  })

  # Count nucleotides
  observeEvent(input$count_nucleotides, {
    dna_sequence <- input$dna_sequence
    nucleotide_counts <- count_nucleotides(dna_sequence)
    output$nucleotide_counts_output <- renderPrint({
      paste("Nucleotide Counts:")
      print(nucleotide_counts)
    })
  })

  # Translate sequence
  observeEvent(input$translate_sequence, {
    dna_sequence <- input$dna_sequence
    translated_sequence <- translate_sequence(dna_sequence)
    output$translated_sequence_output <- renderPrint({
      paste("Translated Sequence:")
      print(translated_sequence)
    })
  })
}

# Run the Shiny app
shinyApp(ui, server)
