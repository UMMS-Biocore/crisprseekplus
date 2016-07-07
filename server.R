library(shiny)
library(CRISPRseek)
library("BSgenome.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(shinyjs)

shinyServer(function(input, output) {
  
  inputFilePath <- system.file("extdata", "inputseq.fa", package = "CRISPRseek")
  inputFile1Path <- system.file("extdata", "rs362331T.fa", package = "CRISPRseek")
  inputFile2Path <- system.file("extdata", "rs362331C.fa", package = "CRISPRseek")
  REpatternFile <- system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek")

  output$contents1 <- renderTable({
    inFile <- input$file1
    if (is.null(inFile))
      return(read.csv(inputFilePath))
    read.csv(inFile$datapath)
  })
  
  output$contents2 <- renderTable({
    inFile2 <- input$file2
    if (is.null(inFile2))
      return(read.csv(REpatternFile))
    read.csv(inFile2$datapath)
  })
  
  output$contents3 <- renderTable({
    inFile3 <- input$file3
    if (is.null(inFile3))
      return(read.csv(inputFile1Path))
    
    read.csv(inFile3$datapath)
  })
  
  output$contents3 <- renderTable({
    inFile4 <- input$file4
    if (is.null(inFile4))
      return(read.csv(inputFile2Path))
    
    read.csv(inFile4$datapath)
  })
  
  
  
output$output1 <- renderUI({
    #Output Directory
  isolate(  
  if(input$givenOutputDir != "")
      outputDir <- input$givenOutputDir
    else (
      outputDir <- "~/CRISPRout"
    ))
    
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

    #Run off target analysis
    resultsOTA <- eventReactive(input$goButton, {
        offTargetAnalysis(inputFilePath = inputFilePath, findgRNAsWithREcutOnly=findgRNAsWithREcutOnly,
                        REpatternFile =REpatternFile,findPairedgRNAOnly=findPairedgRNAOnly,
                        BSgenomeName = BSgenomeName, txdb=txdb,
                        orgAnn = orgAnn,max.mismatch = input$mismatch, chromToSearch = chromToSearch,
                        outputDir = outputDir,overwrite = TRUE, 
                        overlap.gRNA.positions=overlap.gRNA.positions, gRNA.size = gRNA.size,
                        baseBeforegRNA =baseBeforegRNA, baseAfterPAM=baseAfterPAM,
                        max.gap = max.gap, annotateExon = annotateExon, enable.multicore = enable.multicore,
                        PAM.size = PAM.size, temperature = temperature, minREpatternSize = input$REPatSize1,
                        findgRNAs = findgRNAs)
      
    })
    
    #run compare 2 sequences
    resultsC2S <- eventReactive(input$goButton, {
      compare2Sequences(inputFile1Path, inputFile2Path, REpatternFile = REpatternFile, outputDir = outputDir, 
                        overwrite = TRUE, BSgenomeName = BSgenomeName, findgRNAsWithREcutOnly=findgRNAsWithREcutOnly, 
                        findPairedgRNAOnly=findPairedgRNAOnly,  overlap.gRNA.positions=overlap.gRNA.positions, 
                       gRNA.size = gRNA.size, baseBeforegRNA =baseBeforegRNA, baseAfterPAM=baseAfterPAM, 
                       min.gap = min.gap, max.gap = max.gap, PAM.size = PAM.size, temperature = temperature,
                        minREpatternSize = input$REPatSize2, findgRNAs = findgRNAs, max.mismatch = input$mismatch,
                       searchDirection = searchDirection)
    })
    
    if(isolate(input$chooseAction == 1)) {
       resultsOTA()
     }
     else if(isolate(input$chooseAction == 2)) {
       resultsC2S()
     }
  
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
  reset("multicore")})





})
