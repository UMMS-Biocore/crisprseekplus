
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
  titlePanel(h1("CRISPRseekPlus")),
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
                           '.csv')),
        conditionalPanel(
          condition <- "input.scoringmethod == 2",
          fileInput("mismatchActivityFile", label = "Choose Mismatch Activity File",
                    accept=c('text/csv', 
                             'text/comma-separated-values,text/plain', 
                             '.csv')))),
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
                           '.csv')),
        conditionalPanel(
          condition <- "input.scoringmethod == 2",
          fileInput("mismatchActivityFile", label = "Choose Mismatch Activity File",
                    accept=c('text/csv', 
                             'text/comma-separated-values,text/plain', 
                             '.csv')))
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
                           '.csv'))),
      
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
                 fluidRow(
                   column(4, textInput("runNum", "Enter Run Name for Output File Title:", ""))
                 ),
                 wellPanel(
                   fluidRow(
                     column(6, actionLink("mainPanel", "Main Options"))),
                   conditionalPanel(
                     condition <- "input.mainPanel%2 == 1",
                     br(), br(),
                     fluidRow(
                       conditionalPanel(
                         condition <- "input.chooseAction == 1 || input.chooseAction == 2",
                         column(3,
                                selectInput("fileFormat", "Format of input file",
                                            choices = list("fasta" =1, "fastq" =2, "bed" = 3), selected =1)),
                         conditionalPanel(
                           condtion <- "input.chooseAction == 2",
                           column(3,
                                  selectInput("fileFormat2", "Format of input file",
                                              choices = list("fasta" =1, "fastq" =2, "bed" = 3), selected =1))),
                         conditionalPanel(
                           condtion <- "input.fileFormat == 3 || input.fileFormat2 == 3",
                           column(3,
                                  radioButtons("fileHeader", "Does input file contain a header?",
                                               choices = list("Yes" = 1, "No" =2), selected = 2)
                           )),
                         column(3,
                                selectInput("gRNAexport", "Export potential gRNAs in which format?",
                                            choices = list("fasta" =1, "genbank" =2, "both" = 3, "none" = 4), selected =3))
                       )),
                     fluidRow(
                       column(4,
                              selectInput("organism", "Organism:",
                                          choices = list("hg19" = 1, "mm10" = 2,
                                                         "ce6" = 3, "rn5" = 4, 
                                                         "dm3" = 5), selected = 1)),
                       column(4,
                              radioButtons("overwriteFile", "Overwrite the existing files in output directory?",
                                           choices = c("Yes" = 1, "No" = 2), selected = 1))),
                     br(),
                     fluidRow(
                       conditionalPanel(
                         condtion <- ("input.chooseAction == 1"),
                         column(4,
                                textInput("chromSearch", "Chromosome to search:", value = "chrX")),
                         column(4,
                                radioButtons("annPaired1", "Annotate paired information?",
                                             choices = c("Yes" = 1, "No" = 2), selected = 1))
                       ),
                       conditionalPanel(
                         condtion <- ("input.chooseAction == 2"),
                         column(4,
                                radioButtons("annPaired2", "Annotate paired information?",
                                             choices = c("Yes" = 1, "No" = 2), selected = 2))),
                       conditionalPanel(
                         condtion <- ("input.chooseAction == 3"),
                         column(4,
                                numericInput("umicol", "Index of column containing the umi/first 
                                             few bases of sequence R1 reads",
                                             value = 2)),
                         column(4,
                                radioButtons("umiheader", "Does input umi file contain a header line?",
                                             choices = list("Yes" = 1, "No" = 2),
                                             selected = 2)),
                         column(4,
                                numericInput("readID", "Index of umi file column containing the read identifier",
                                             value = 1))),
                       conditionalPanel(
                         condtion <- "input.chooseAction == 1",
                         column(4,
                                radioButtons("multicore", "Enable parallel processing?",
                                             choices = list("Yes" = 1, "No" = 2),
                                             selected = 2))))
                   )), #FIRST WELL PANEL
                 
                 
                 
                 
                 #################### GRNA PANEL #######################
                 wellPanel(
                   fluidRow(
                     column(6, actionLink("gRNApanel", "gRNA Advanced Options"))),
                   conditionalPanel(
                     condition <- "input.gRNApanel%2 == 1",
                     br(), br(),
                     fluidRow(
                       conditionalPanel(
                         condtion <- ("input.chooseAction == 1"),
                         column(4,
                                radioButtons("radio1", "Find gRNAS with RE cut only?",
                                             choices = list("Yes" = 1, "No" = 2),
                                             selected = 1)),
                         column(4,
                                radioButtons("radio2", "Find  paired gRNAS only?",
                                             choices = list("Yes" = 1, "No" = 2),
                                             selected = 2))),
                       conditionalPanel(
                         condtion <- ("input.chooseAction == 2"),
                         column(4,
                                radioButtons("radio3", "Find gRNAS with RE cut only?",
                                             choices = list("Yes" = 1, "No" = 2),
                                             selected = 2)),
                         column(4,
                                radioButtons("radio4", "Find  paired gRNAS only?",
                                             choices = list("Yes" = 1, "No" = 2),
                                             selected = 2))),
                       column(4,
                              numericInput("gRNASize", "Enter gRNA size", value = 20))),
                     fluidRow(
                       column(12,
                              sliderInput("overlapgRNA", "Set overlap positions of gRNA and
                                          restriction enzyme cut site", min = 0, max = 50,                   
                                          c(17, 18)))),
                     fluidRow(
                       conditionalPanel(
                         condtion <- "input.chooseAction == 1 || input.chooseAction == 2",
                         column(4,
                                numericInput("baseBefore", "Number bases before gRNA for efficiency", value = 4)),
                         column(4,
                                numericInput("minGap", "Minimum distance between two oppositely oriented 
                                             gRNAs to be validly paired", value = 0)),
                         column(4,
                                numericInput("maxGap", "Maximum distance between two oppositely oriented 
                                             gRNAs to be validly paired", value = 20)))),
                     fluidRow(
                       conditionalPanel(
                         condtion <- "input.chooseAction == 2",
                         column(4,
                                radioButtons("removeDetails", "Remove detailed gRNA information?",
                                             choices = list("Yes(both files)" = 1, "No(both files)" = 2,
                                                            "Yes(file1)/No(file2)" = 3, "No(file1)/Yes(file2)" = 4),
                                             selected = 2)),
                         column(4,
                                radioButtons("findgRNA2", "Find gRNAs from input sequences?",
                                             choices = c("Yes (both inputs)" = 1, "No (both inputs)" = 2,
                                                         "Yes(input1)/No(input 2)" = 3, 
                                                         "No(input1)/Yes(input2)" = 4), selected = 1)),
                         column(4,
                                checkboxGroupInput("searchDir", "Search direction",
                                                   choices = list("Both" = "both", "1to2" = "1to2", "2to1" = "2to1"),
                                                   selected = c("both", "1to2", "2to1")))),
                       conditionalPanel(
                         condition <- "input.chooseAction == 1",
                         column(4,
                                radioButtons("findgRNA1", "Find gRNAs from input sequences?",
                                             choices = c("Yes" = 1, "No" = 2), selected = 1)),
                         column(4,
                                radioButtons("usescore", "Use score to indicate gRNA efficacy?",
                                             choices = list("Yes" = 1, "No" = 2), selected = 1)),
                         conditionalPanel(
                           condition <- "input.usescore == 1",
                           column(3, 
                                  radioButtons("efficacyFromIS", "Use input sequences to calculate gRNA efficacy?",
                                               choices = list("Yes" = 1, "No" = 2), selected = 2))))),
                     fluidRow(
                       conditionalPanel(
                         condtion <- "input.chooseAction == 1 || input.chooseAction == 2",
                         column(4,
                                textInput("gRNAname", "Enter a prefix used when assigning a name to found gRNAs", "gRNA"))),
                       conditionalPanel(
                         condition <- "input.chooseAction == 1",
                         column(3,
                                numericInput("upstreamSearch","Upstream offset from the bed input starts for gRNA search",
                                             value = 0)),
                         column(3,
                                numericInput("downstreamSearch","Downstream offset from the bed input ends for gRNA search",
                                             value = 0))))
                                )), ###gRNA Settings
                 
                 
                 
                 
                 
                 #################### PAM PANEL #######################
                 wellPanel(
                   fluidRow(
                     column(6, actionLink("PAMpanel", "PAM Advanced Options"))),
                   conditionalPanel(
                     condition <- "input.PAMpanel%2 == 1",
                     br(), br(),
                     fluidRow(
                       column(4,
                              numericInput("pamSize", "Set PAM length", value = 3)),
                       column(4, 
                              numericInput("PAMmismatch", "Number degenerative bases in the PAM",
                                           value = 2)),
                       conditionalPanel(
                         condition <- "input.chooseAction == 1 || input.chooseAction == 2",
                         column(4,
                                numericInput("baseAfter", "Number bases after PAM for efficiency", value = 3)))
                     ),
                     fluidRow(
                       column(4,
                              textInput("PAMSeq", "Enter PAM Sequence after gRNA", "NGG")),
                       conditionalPanel(
                         condition <- "input.chooseAction == 1 || input.chooseAction == 2",
                         column(4,
                                textInput("PAMpattern", "Enter PAM pattern for off target", "N[A|G]G$"))
                       ),
                       conditionalPanel(
                         condition <- "input.chooseAction == 3",
                         column(4,
                                textInput("PAMpattern", "Enter PAM pattern for off target", "(NAG|NGG|NGA)$"))))
                     
                     
                   )), ###PAM Settings
                 
                 
                 #################### MORE ADVANCED PANEL #######################
                 wellPanel(
                   fluidRow(
                     column(6, actionLink("advanced", "Other Advanced Options"))),
                   conditionalPanel(
                     condition <- "input.advanced%2 == 1",
                     br(),
                     
                     conditionalPanel(
                       condtion <- "input.chooseAction == 1 || input.chooseAction == 2",
                       fluidRow(
                         column(4,
                                numericInput("temp", "Set temperature (celsius)", value = 37)),
                         conditionalPanel(
                           condtion <- "input.chooseAction == 1",
                           column(4,
                                  radioButtons("annExon", "Indicate whether off target is inside an exon?",
                                               choices = list("Yes" = 1, "No" = 2),
                                               selected = 1)),
                           column(4,
                                  numericInput("REPatSize1", "Minimum RE Pattern Size", value = 4))),
                         conditionalPanel(
                           condition <- "input.chooseAction == 2",
                           column(4,
                                  numericInput("REPatSize2", "Minimum RE Pattern Size", value = 6))))),
                     
                     conditionalPanel(
                       condition <- "input.chooseAction == 3",
                       fluidRow(
                         column(4,
                                numericInput("minR1", "Max mapped R1 bp length to be considered for downstream analysis",
                                             value = 30)),
                         column(4,
                                numericInput("minR2", "Max mapped R2 bp length to be considered for downstream analysis",
                                             value = 30)),
                         column(4,
                                radioButtons("applyMinMapped", "Apply min. mapped length requirements to both R1 and R2",
                                             choices = list("Yes" = 1, "No" = 2),
                                             selected = 2))),
                       fluidRow(
                         column(4,
                                radioButtons("concordantStrand", 
                                             "R1 and R2 aligned to same strand or opposite strand?",
                                             choices = list("Same Strand" = 1, "Opposite Strand" = 2),
                                             selected = 2)),
                         column(4,
                                numericInput("maxPairedDistance", "Maximum bp distance between paired R1 and R2 reads",
                                             value = 1000)),
                         column(4,
                                numericInput("minMapQuality", "Minimum mapping quality of acceptable alignments",
                                             value = 30))),
                       fluidRow(
                         column(4,
                                numericInput("distInterChrom", 
                                             "Distance value to assign to paired reads that are alighned to different
                                             chromosome",
                                             value = -1)),
                         column(4,
                                radioButtons("sameChrom", "Paired reads aligned to same chromosome?",
                                             choices = list("Yes" = 1, "No" = 2),
                                             selected = 1)),
                         column(4,
                                numericInput("minReads", "Minimum number of reads to be considered a peak",
                                             value = 20L))),
                       fluidRow(
                         column(4,
                                numericInput("max.P", 
                                             "Maximum p-value to be considered as significant",
                                             value = 0.05)),
                         column(4,
                                selectInput("stat", "Statistical test choice",
                                            choices = list("poisson" = 1, "nbinom" = 2),
                                            selected = 1)),
                         column(4,
                                numericInput("distThreshold", "Maximum gap allowed between plus and negative strand peak",
                                             value = 40))),
                       fluidRow(
                         radioButtons("adjustMethods", "Adjustment method for multiple comparisons",
                                      choices = list("None" = 1, "BH" = 2, "holm" = 3, "hochberg" =4,
                                                     "hommel"= 5, "bonferroni" =6, "BY" = 7, "fdr" = 8),
                                      selected = 1, inline = TRUE))),
                     br(),
                     br(),
                     
                     conditionalPanel(
                       condition <- "input.chooseAction == 1",
                       fluidRow(
                         column(4, 
                                numericInput("minScore", "Minimum score of an off target to be included in the final output",
                                             value = 0)),
                         column(4, 
                                numericInput("top.N", "Top N off targets to be included in the final output",
                                             value = 1000)),
                         column(4, 
                                numericInput("topNscore", "Top N off target used to calculate the total off target score",
                                             value = 10))),
                       fluidRow(
                         column(3, 
                                radioButtons("fetchSeq", "Fetch flank sequence of off target?",
                                             choices = list("Yes" = 1, "No" = 2), selected = 1)),
                         conditionalPanel(
                           condition <- "input.fetchSeq == 1",
                           column(3,
                                  numericInput("up.stream", "Upstream offset from the off target start", 200)),
                           column(3,
                                  numericInput("down.stream", "Downstream offset from the off target start", 200))))),
                     
                     fluidRow(
                       conditionalPanel(
                         condition <- "input.chooseAction == 2",
                         column(3,
                                numericInput("up.stream", "Upstream offset from the off target start", 0)),
                         column(3,
                                numericInput("down.stream", "Downstream offset from the off target start", 0))),
                       
                       conditionalPanel(
                         condition <- "input.chooseAction == 3",
                         column(3,
                                numericInput("up.stream", "Upstream offset from the off target start", 50)),
                         column(3,
                                numericInput("down.stream", "Downstream offset from the off target start", 50)))),
                     
                     conditionalPanel(
                       condition <- "input.chooseAction == 1 || input.chooseAction == 2",
                       fluidRow(
                         selectInput("scoringmethod", "Select method to use for offtarget cleavage rate
                                     estimation", choices = list("Hsu-Zhang" = 1, "CFDscore"= 2),
                                     selected = 1))),
                     fluidRow(
                       conditionalPanel(
                         condition <- "input.scoringmethod == 1",
                       column(8, 
                              print(strong("Weights")), br(),
                              tags$textarea(id="weight", rows = 2, cols = 100, 
                              "0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583"))
                     )),
                     
                     conditionalPanel(
                       condition <- "input.scoringmethod == 2",
                       br(),
                       paste("Please see side panel for Mismatch Activity File Upload"),
                       br(), br(),
                       fluidRow(
                         column(6,
                                numericInput("subPamPos1", "Start Position of Sub PAM", value = 22)),
                         column(6,
                                numericInput("subPamPos2", "End Position of Sub PAM", value = 23))),
                       print("Sub PAM Activity Input:"),
                       br(),
                       column(3, numericInput("AA", "AA", value = 0)),
                       column(3, numericInput("AC", "AC", value = 0)),
                       column(3, numericInput("AG", "AG", value = 0.259259259)),
                       column(3, numericInput("AT", "AT", value = 0)),
                       column(3, numericInput("CA", "CA", value = 0)),
                       column(3, numericInput("CC", "CC", value = 0)),
                       column(3, numericInput("CG", "CG", value = 0.107142857)),
                       column(3, numericInput("CT", "CT", value = 0)),
                       column(3, numericInput("GA", "GA", value = 0.0694444440)),
                       column(3, numericInput("GC", "GC", value = 0.022222222)),
                       column(3, numericInput("GG", "GG", value = 1)),
                       column(3, numericInput("GT", "GT", value = 0.016129032)),
                       column(3, numericInput("TA", "TA", value = 0)),
                       column(3, numericInput("TC", "TC", value = 0)),
                       column(3, numericInput("TG", "TG", value = 0.038961039)),
                       column(3, numericInput("TT", "TT", value = 0))),
                     
                     fluidRow(
                       conditionalPanel(
                         condition <- "input.chooseAction == 3",
                         column(4,
                                numericInput("window", "Window size to calculate coverage", value = 20)),
                         column(4,
                                numericInput("stepSize", "Step size to calculate coverage", value = 20)),
                         column(4,
                                numericInput("BGWindow", "Background Window Size", value = 5000)))),
                     fluidRow(
                       column(4,
                              numericInput("mismatch", "Maximum mismatch:", min = 0, value = 3)))
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
        ),#Tab for Data Table Panel
        
        tabPanel("View Input Files",
                 conditionalPanel(
                   condition <- "input.chooseAction == 1",
                   radioButtons("fileViewChoices", "Choose file to view",
                                choices = list("Input File" = 1, "Pattern File" =2,
                                               "Mismatch Activity File" = 3), selected = 1),
                   conditionalPanel(
                     condition <- "input.fileViewChoices == 1",
                     textOutput("fileInput1")
                   ),
                   conditionalPanel(
                     condition <- "input.fileViewChoices == 2",
                     DT::dataTableOutput("fileInput2")
                   )
                 ),
                 
                 conditionalPanel(
                   condition <- "input.chooseAction == 2",
                   radioButtons("fileViewChoices2", "Choose file to view",
                                choices = list("Input File 1" = 1, "Input File 2" =2,
                                               "Mismatch Activity File" = 3), selected = 1),
                   
                   conditionalPanel(
                     condition <- "input.fileViewChoices2 == 1",
                     textOutput("fileInput3")
                   ),
                   
                   conditionalPanel(
                     condition <- "input.fileViewChoices2 == 2",
                     textOutput("fileInput4")
                   )
                 ),
                 
                 conditionalPanel(
                   condition <- "input.chooseAction == 3",
                   radioButtons("fileViewChoices3", "Choose file to view",
                                choices = list("UMI File" = 1, "Alignment File" = 2, "gRNA File" = 3), selected =1),
                   conditionalPanel(
                     condition <- "input.fileViewChoices3 == 1",
                     DT::dataTableOutput("fileInput5")
                   ),
                   conditionalPanel(
                     condition <- "input.fileViewChoices3 == 2",
                     print("Alignment file too large to view")
                   ),
                   conditionalPanel(
                     condition <- "input.fileViewChoices3 == 3",
                     textOutput("fileInput7")
                   )
                 )))#Tabset Panel
)#mainPanel
)#sidebarLayout
)
)


