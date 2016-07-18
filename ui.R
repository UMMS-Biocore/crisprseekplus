
library(shiny)
library(CRISPRseek)
library("BSgenome.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(shinyjs)
library(DT)
library(GUIDEseq)

shinyUI(fluidPage(
  shinyjs::useShinyjs(),
  
  conditionalPanel(
    condition <- "input.goButton > 0",
    uiOutput("loading")
  ),
  uiOutput("logo"),
  titlePanel(h1("CRISPRSeeker")),
  br(),
  br(),
  br(),
  sidebarLayout(
    #File Upload
  sidebarPanel(
    radioButtons("chooseAction", "Choose which function you would like to do",
                 choices = list("Off Target Analysis" = 1,
                            "Compare 2 Sequences" = 2, "GUIDEseq" = 3), selected = "1"),
    
      conditionalPanel(
      condition <- "input.chooseAction == 1",
      paste("Off Target Analysis File Upload"),
      br(),
      fileInput("file1", label = "Choose Input File",
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
      fileInput("file2", label = "Choose Pattern File",
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv'))
    ),
    conditionalPanel(
      condition <- "input.chooseAction == 2",
      paste("Compare 2 Sequences File Upload"),
      br(),
      fileInput("file3", label = "Choose Input File 1",
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
      fileInput("file4", label = "Choose Input File 2",
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv'))
    ),
    conditionalPanel(
      condition <- "input.chooseAction == 3",
      paste("GUIDE-Seq File Upload"),
      br(),
      fileInput("file5", label = "Choose UMI File",
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
      fileInput("file6", label = "Choose Input Alignment File",
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
      fileInput("file7", label = "Choose Input gRNA File",
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv'))
    ),
    
    helpText(a(strong("Help Page"), 
               href="http://crisprseeker.readthedocs.io/en/develop/index.html",
               TARGET="_blank")),
    helpText(a("What is Off Target Analysis?", 
                  href="http://crisprseeker.readthedocs.io/en/latest/quickstart.html#what-is-off-target-analysis",
               TARGET="_blank")),
    helpText(a("What is Compare 2 Sequences?", 
               href="http://crisprseeker.readthedocs.io/en/latest/quickstart.html#what-is-compare-2-sequences",
               TARGET="_blank")),
    helpText(a("What is GUIDEseq?", 
               href="http://crisprseeker.readthedocs.io/en/latest/quickstart.html#what-is-guide-seq-analysis",
               TARGET="_blank")),
    helpText(a("How To Use the Interface", 
               href=" http://crisprseeker.readthedocs.io/en/develop/quickstart.html#using-the-interface",
               TARGET="_blank")),
    helpText(a("About Files", 
               href=" http://crisprseeker.readthedocs.io/en/develop/quickstart.html#about-files",
               TARGET="_blank")),
    br(),
    
    downloadButton("downloadData", "Download Output")

    ),

    #Data Entry
    mainPanel( 
    tabsetPanel(
      tabPanel("Submissions Panel",
    wellPanel(
    fluidRow(
    conditionalPanel(
      condtion <- ("input.chooseAction == 1"),
        column(4,
             radioButtons("radio1", "Find gRNAS with RE cut only?",
                                choices = list("Yes" = 1, "No" = 2),
                                selected = 1))),
    conditionalPanel(
      condtion <- ("input.chooseAction == 1"),
      column(4,
             radioButtons("radio2", "Find  paired gRNAS only?",
                          choices = list("Yes" = 1, "No" = 2),
                          selected = 2))),
    conditionalPanel(
      condtion <- ("input.chooseAction == 2"),
      column(4,
             radioButtons("radio3", "Find gRNAS with RE cut only?",
                          choices = list("Yes" = 1, "No" = 2),
                          selected = 2))),
    conditionalPanel(
      condtion <- ("input.chooseAction == 2"),
      column(4,
             radioButtons("radio4", "Find  paired gRNAS only?",
                          choices = list("Yes" = 1, "No" = 2),
                          selected = 2))),
      column(4,
             selectInput("organism", "Organism:",
                         choices = list("hg19" = 1, "mm10" = 2,
                                        "ce6" = 3, "rn5" = 4, 
                                        "dm3" = 5), selected = 1)),
      conditionalPanel(
        condition <- "input.chooseAction == 3",
        column(4,
        numericInput("minR1", "Max mapped R1 bp length to be considered for downstream analysis",
                     value = 30)),
        column(4,
               numericInput("minR2", "Max mapped R2 bp length to be considered for downstream analysis",
                            value = 30)))),
    br(),
      fluidRow(
        column(4,
               numericInput("mismatch", "Max mismatch:", min = 0, value = 3)),
        conditionalPanel(
          condtion <- ("input.chooseAction == 1"),
          column(4,
               textInput("chromSearch", "Chromosome to search:", value = "chrX"))
      ),
      conditionalPanel(
        condtion <- ("input.chooseAction == 3"),
        column(4,
               radioButtons("sameChrom", "Paired reads aligned to same chromosome?",
                            choices = list("Yes" = 1, "No" = 2),
                            selected = 1))
      )
      )),

    
  
      
    #Advanced Options
    wellPanel(
    fluidRow(
    column(6, actionLink("advanced", "Advanced Options"))),
    conditionalPanel(
      condtion <- "input.advanced > 0",
      br(),
      fluidRow(
        column(12,
               sliderInput("overlapgRNA", "Set overlap positions of gRNA and
                           restriction enzyme cut site", min = 0, max = 50,                    c(17, 18)))),
      fluidRow(
        column(4,
               numericInput("gRNASize", "Enter gRNA size", value = 20)),
        conditionalPanel(
          condtion <- "input.chooseAction == 1 || input.chooseAction == 2",
          column(4,
                 numericInput("baseBefore", "Number bases before gRNA for efficiency", value = 4)),
          column(4,
                 numericInput("baseAfter", "Number bases after PAM for efficiency", value = 3))),
          conditionalPanel(
            condition <- "input.chooseAction == 3",
            column(4,
                 numericInput("pamSize", "Set PAM length", value = 3)))
        ),
      conditionalPanel(
        condtion <- "input.chooseAction == 1 || input.chooseAction == 2",
        fluidRow(
          column(4,
                 numericInput("minGap", "Minimum distance between two oppositely oriented 
                            gRNAs to be validly paired", value = 0)),
          column(4,
                 numericInput("maxGap", "Maximum distance between two oppositely oriented 
                            gRNAs to be validly paired", value = 20)),
          conditionalPanel(
            condtion <- "input.chooseAction == 1",
            column(4,
                   radioButtons("annExon", "Indicate whether off target is inside an exon?",
                                choices = list("Yes" = 1, "No" = 2),
                                selected = 1))),
          conditionalPanel(
            condtion <- "input.chooseAction == 2",
            column(4,
                   checkboxGroupInput("searchDir", "Search direction",
                                      choices = list("Both" = "both", "1to2" = "1to2", "2to1" = "2to1"),
                                      selected = c("both", "1to2", "2to1"))))
        )),
      conditionalPanel(
        condtion <- "input.chooseAction == 1 || input.chooseAction == 2",
        fluidRow(
          column(4,
                 numericInput("pamSize", "Set PAM length", value = 3)),
          column(4,
                 numericInput("temp", "Set temperature (celsius)", value = 37)),
          conditionalPanel(
            condtion <- "input.chooseAction == 1",
            column(4,
                   radioButtons("multicore", "Enable parallel processing?",
                                choices = list("Yes" = 1, "No" = 2),
                                selected = 2))))),
      fluidRow(
        conditionalPanel(
          condition <- "input.chooseAction == 1",
          column(4,
                 numericInput("REPatSize1", "Minimum RE Pattern Size", value = 4)),
          column(4,
                 radioButtons("findgRNA1", "Find gRNAs from input sequences?",
                              choices = c("Yes" = 1, "No" = 2), selected = 1)),
          column(4,
                 radioButtons("annPaired1", "Annotate paired information?",
                              choices = c("Yes" = 1, "No" = 2), selected = 1))),
        conditionalPanel(
          condition <- "input.chooseAction == 2",
          column(4,
                 numericInput("REPatSize2", "Minimum RE Pattern Size", value = 6)),
          column(4,
                 radioButtons("findgRNA2", "Find gRNAs from input sequences?",
                              choices = c("Yes (both inputs)" = 1, "No (both inputs)" = 2,
                                          "Yes(input1)/No(input 2)" = 3, 
                                          "No(input1)/Yes(input2)" = 4), selected = 1)),
          column(4,
                 radioButtons("annPaired2", "Annotate paired information?",
                              choices = c("Yes" = 1, "No" = 2), selected = 2))
          
        )
      ),
      
      fluidRow(
        column(4,
              textInput("PAMSeq", "Enter PAM Sequence", "NGG")),
        column(4,
               radioButtons("overwriteFile", "Overwrite the existing files in output directory?",
                            choices = c("Yes" = 1, "No" = 2), selected = 1))
      )
      )#Well Panel for advanced settings
    ),
    fluidRow(
      column(4, actionButton("goButton", "Submit"))
    ),
    br(),
    fluidRow(
      column(4, actionButton("resetFields", "Reset all fields to defaults"))
    ),
    uiOutput('output1')
    ),#Tab for Submission Panel
    
    tabPanel("Data Table",
             conditionalPanel(
               condtion <- "input.chooseAction == 1",
               titlePanel("RE Cut Details")
             ),
             conditionalPanel(
               condtion <- "input.chooseAction == 2",
               titlePanel("Scores for 2 Input Sequences")
             ),
             conditionalPanel(
               condtion <- "input.chooseAction == 3",
               titlePanel("gRNA Peaks")
             ),
             DT::dataTableOutput("tables")
             )#Tab for Data Table Panel
    )#Tabset Panel
  )#mainPanel
  )#sidebarLayout
))


