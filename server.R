library(shiny)
library(CRISPRseek)
library("BSgenome.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(shinyjs)
library(DT)
library(GUIDEseq)

shinyServer(function(input, output) {
  #' getLoadingMsg
  #'
  #' @note \code{getLoadingMsg}
  #' @return loading msg
  #' @examples
  #'    x <- getLoadingMsg()
  #' @export
  #'
  getLoadingMsg <- function() {
    imgsrc <- "www/images/loading.gif"
    a <- list(
      tags$head(tags$style(type = "text/css", "
            #loadmessage {
            position: fixed;
            top: 0px;
            left: 200px;
            width: 70%;
            height: 100;
            padding: 5px 0px 5px 0px;
            text-align: center;
            font-weight: bold;
            font-size: 100%;
            color: #000000;
            opacity: 0.8;
            background-color: #FFFFFFF;
            z-index: 100;
            }")),
      conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
                       tags$div("Please wait! Loading...", id = "loadmessage",
                                tags$img(src = imgsrc
                                ))))
  }
  
  
  #' getLogo
  #'
  #' Generates and displays the logo to be shown within DEBrowser.
  #'
  #' @note \code{getLogo}
  #' @return return logo
  #' @examples
  #'     x <- getLogo()
  #' @export
  #'
  getLogo <- function(){
    imgsrc <- "www/images/logo.png"
    a<-list(img(src=imgsrc, align = "left"))
  }
  
  output$loading <- renderUI({
    getLoadingMsg()
  })
  output$logo <- renderUI({
    getLogo()
  })
  
  #Enable or disable to download button depending on if 
  #analysis is complete
  disableDownload <- function() {
    toggleState(
      id = "downloadData",
      condition = input$goButton > 0
    )
  }
  
output$output1 <- renderUI({
  
  #InputFilePath
  if(is.null(input$file1)) {
    inputFilePath <- system.file("extdata", "inputseq.fa", package = "CRISPRseek")
  }
  else {
    inputFilePath <- input$file1$datapath
  }
  #REpatternFile
  if(is.null(input$file2)) {
    REpatternFile <- system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek")
  }
  else {
    REpatternFile <- input$file2$datapath
  }
  #InputFile1Path for C2S
  if(is.null(input$file3)) {
    inputFile1Path <- system.file("extdata", "rs362331T.fa", package = "CRISPRseek")
  }
  else {
    inputFile1Path <- input$file3.datapath
  }
  #InputFile2Path for C2S
  if(is.null(input$file4)) {
    inputFile2Path <- system.file("extdata", "rs362331C.fa", package = "CRISPRseek")
  }
  else {
    inputFile2Path <- input$file4.datapath
  }
  
  if(is.null(input$file5)) {
    umifile <- system.file("extdata", "UMI-HEK293_site4_chr13.txt", package = "GUIDEseq")
  }
  else {
    umifile <- input$file5$datapath
  }
  
  if(is.null(input$file6)) {
    bamfile <- system.file("extdata","bowtie2.HEK293_site4_chr13.sort.bam", package = "GUIDEseq")
  }
  else {
    bamfile <- input$file6$datapath
  }
  
  if(is.null(input$file7)) {
    gRNA.file <- system.file("extdata","gRNA.fa", package = "GUIDEseq")
  }
  else {
    gRNA.file <- input$file7$datapath
  }

  
  #Output Directory
  isolate(  
      outputDir <- tempdir()
    )
    #Find gRNAs With RE Cut Only?
  isolate(
    if(input$chooseAction == 1) {
    if(input$radio1 == 2)
      findgRNAsWithREcutOnly<- FALSE
    else (
      findgRNAsWithREcutOnly<- TRUE
    )
  }
  else {
    if(input$radio3 == 2)
      findgRNAsWithREcutOnly<- FALSE
    else (
      findgRNAsWithREcutOnly<- TRUE
      )
    })

    #Find Paired gRNA only?
  isolate(
  if(input$chooseAction == 1) {
    if(input$radio2 == 1)
      findPairedgRNAOnly <- TRUE
    else (
      findPairedgRNAOnly <- FALSE
    )}
  else {
    if(input$radio4 == 1)
      findPairedgRNAOnly <- TRUE
    else (
      findPairedgRNAOnly <- FALSE
    )})
    
    #Input Organism: updates txdb, BSGenomeName, and orgAnn
    orgAnn <- org.Hs.egSYMBOL
    BSgenomeName <- Hsapiens
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    isolate(
    if(input$organism == "mm10") {
      library("BSgenome.Mmusculus.UCSC.mm10")
      library(TxDb.Mmusculus.UCSC.mm10.knownGene)
      library(org.Mm.eg.db)
      BSgenomeName= Mmusculus
      txdb=TxDb.Mmusculus.UCSC.mm10.knownGene
      orgAnn = org.Mm.egSYMBOL 
    }
    
    else if(input$organism == "ce6") {
      library("BSgenome.Celegans.UCSC.ce6")
      library(TxDb.Celegans.UCSC.ce6.ensGene)
      library(org.Ce.eg.db)
      BSgenomeName= "Celegans"
      txdb=TxDb.Celegans.UCSC.ce6.ensGene
      orgAnn = org.Ce.egSYMBOL 
    }
    else if(input$organism == "rn5") {
      library("BSgenome.Rnorvegicus.UCSC.rn5")
      library(TxDb.Rnorvegicus.UCSC.rn5.refGene)
      library(org.Rn.eg.db)
      BSgenomeName= Rnorvegicus
      txdb=TxDb.Rnorvegicus.UCSC.rn5.refGene
      orgAnn = org.Rn.egSYMBOL
    }
    
    else if(input$organism == "dm3") {
      library("BSgenome.Dmelanogaster.UCSC.dm3")
      library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
      library(org.Dm.eg.db)
      BSgenomeName= Dmelanogaster
      txdb=TTxDb.Dmelanogaster.UCSC.dm3.ensGene
      orgAnn = org.Mm.egSYMBOL 
    })
    #Chromosome to Search?
    isolate(
    chromToSearch <- input$chromSearch
    )
    
    #Advanced Settings 
    
    #Change default gRNA and RE cut size overlap positions
    isolate(
    overlap.gRNA.positions <- c(input$overlapgRNA[1], input$overlapgRNA[2]))
    #Change the input size of the gRNA
    isolate(
    gRNA.size <- input$gRNASize)
    #Change number of bases before the gRNA for calculating gRNA efficiency
    isolate(  
    baseBeforegRNA <- input$baseBefore)
    #Change number of bases before PAM for calculating gRNA efficiency
    isolate(   
    baseAfterPAM <- input$baseAfter)
    #Change min. distance between oppositely orientied gRNA to be valid paired gRNA
    isolate(
    min.gap <- input$minGap)
    #Change max. distance between oppositely orientied gRNA to be valid paired gRNA
    isolate(
    max.gap <- input$maxGap)
    #Choose whether or not to indicate whether the off target is inside and exon or not
    isolate(  
    if(input$annExon == 2)
      annotateExon <- FALSE
      else (
        annotateExon <- TRUE
      ))
    #Change PAM length
    isolate(
    PAM.size <- input$pamSize)
    #Change temperature
    isolate(
    temperature <- input$temp)
    #Enable multicore processing?
    isolate(
      if(input$multicore == 1)
      enable.multicore <- TRUE
    else(
      enable.multicore <- FALSE
    ))
    #Search direction 
    searchDirection <- c()
    isolate(
    append(searchDirection, input$searchDir)
    )

    #Find gRNA from input?
    isolate(
    if(input$chooseAction == 1) {
      if(input$findgRNA1 == 1){
        findgRNAs <- TRUE
      }
      else {
        findgRNAs <- FALSE
      }
    }
    else {
      if(input$findgRNA2 == 1) {
        findgRNAs <- c(TRUE, TRUE)
      }
      else if(input$findgRNA2 == 2) {
        findgRNAs <- c(FALSE, FALSE)
      }
      else if(input$findgRNA2 == 3) {
        findgRNAs <- c(TRUE, FALSE)
      }
      else if(input$findgRNA2 == 4) {
        findgRNAs <- c(FALSE, TRUE)
      }
    })
    #Change default PAM Sequence
    isolate(
      PAM <- input$PAMSeq
    )
    #Overwrite existing filesin the output directory each time analysis is run?
    if(input$overwriteFile == 1) {
      overwrite <- TRUE
    }
    else {
      overwrite <- FALSE
    }
    #Max mapped R1 bp length to be considered for downstream analysis
    min.R1.mapped <- input$minR1
    min.R2.mapped <- input$minR2
    #Specify whether the paired reads are required to align to the same chromosme
    if(input$sameChrom == 1) {
      same.chromosome <- TRUE
    }
    else {
      same.chromosome <- FALSE
    }
    
    #Run off target analysis
    resultsOTA <- eventReactive(input$goButton, {
      offTargetAnalysis(inputFilePath = inputFilePath, findgRNAsWithREcutOnly=findgRNAsWithREcutOnly,
                        REpatternFile =REpatternFile,findPairedgRNAOnly=findPairedgRNAOnly,
                        BSgenomeName = BSgenomeName, txdb=txdb,
                        orgAnn = orgAnn,max.mismatch = input$mismatch, chromToSearch = chromToSearch,
                        outputDir = outputDir, overwrite = overwrite, 
                        overlap.gRNA.positions=overlap.gRNA.positions, gRNA.size = gRNA.size,
                        baseBeforegRNA =baseBeforegRNA, baseAfterPAM=baseAfterPAM,
                        max.gap = max.gap, annotateExon = annotateExon, enable.multicore = enable.multicore,
                        PAM.size = PAM.size, temperature = temperature, minREpatternSize = input$REPatSize1,
                        findgRNAs = findgRNAs, PAM = PAM)
      
    })
    
    #run compare 2 sequences
    resultsC2S <- eventReactive(input$goButton, {
      compare2Sequences(inputFile1Path = inputFile1Path, inputFile2Path =inputFile2Path, REpatternFile = REpatternFile, outputDir = outputDir, 
                        overwrite = overwrite, BSgenomeName = BSgenomeName, findgRNAsWithREcutOnly=findgRNAsWithREcutOnly, 
                        findPairedgRNAOnly=findPairedgRNAOnly,  overlap.gRNA.positions=overlap.gRNA.positions, 
                        gRNA.size = gRNA.size, baseBeforegRNA =baseBeforegRNA, baseAfterPAM=baseAfterPAM, 
                        min.gap = min.gap, max.gap = max.gap, PAM.size = PAM.size, temperature = temperature,
                        minREpatternSize = input$REPatSize2, findgRNAs = findgRNAs, max.mismatch = input$mismatch,
                        searchDirection = searchDirection, PAM = PAM)
    })
    
    #run GUIDEseqAnalysis
    resultsGSA <- eventReactive(input$goButton, {
      GUIDEseqAnalysis(alignment.inputfile = bamfile, umi.inputfile = umifile,
                       BSgenomeName = BSgenomeName, gRNA.file = gRNA.file,
                       outputDir = outputDir, overlap.gRNA.positions = overlap.gRNA.positions,
                       PAM.size = PAM.size, gRNA.size = gRNA.size,PAM = PAM, max.mismatch = input$mismatch,
                       overwrite = overwrite, min.R1.mapped = min.R1.mapped, min.R2.mapped = min.R2.mapped,
                       same.chromosome = same.chromosome)})
    
    
    setwd(outputDir)
    disableDownload()
    if(isolate(input$chooseAction == 1)) {
       resultsOTA()
     }
     else if(isolate(input$chooseAction == 2)) {
       resultsC2S()
     }
    else if(isolate(input$chooseAction == 3)) {
      resultsGSA()
    }

    #Data Tables
    output$tables <- DT::renderDataTable(DT::datatable({
      if(input$goButton < 1) {
        return()
      }
      else {
        if(input$chooseAction == 1) {
          data <- read.table(paste0(outputDir, "/RECutDetails.xls"),
                             header = TRUE)
        }
        else if(input$chooseAction == 2) {
          data <- read.table(paste0(outputDir, "/scoresFor2InputSequences.xls"),
                             colClasses = c(NA, NA, NA, NA, NA, "NULL", "NULL", "NULL", "NULL", "NULL", 
                                            "NULL", "NULL", "NULL", "NULL"), header = TRUE)
          
        }
        else if(input$chooseAction == 3) {
          data <- read.table(paste0(outputDir, "/gRNA-peaks.xls"), header = TRUE)
        }
      }
    }))#Data Table
    #Download output as zip file
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("crisprSeekOutput","zip", sep=".")
      },
      content = function(fname) {
        crisprOutput <- outputDir
        
        zip(zipfile = fname, files = crisprOutput)
      },
      contentType = "application/zip"
    )
      return()
  })

#Reset all fields
observeEvent(input$resetFields, {
  reset("givenOutputDir") 
  reset("radio1") 
  reset("radio2") 
  reset("radio3") 
  reset("radio4")
  reset("organism")
  reset("mismatch") 
  reset("chromSearch") 
  reset("overlapgRNA") 
  reset("gRNASize")
  reset("baseBefore") 
  reset("baseAfter") 
  reset("minGap") 
  reset("maxGap") 
  reset("annExon") 
  reset("searchDir") 
  reset("pamSize") 
  reset("REPatSize1") 
  reset("REPatSize2") 
  reset("findgRNA1") 
  reset("annPaired1") 
  reset("temp") 
  reset("findgRNA2")
  reset("annPaired2") 
  reset("multicore")
  reset("PAMSeq")
  reset("overwriteFile")
  reset("minR1")
  reset("minR2")
  reset("sameChrom")})
})




