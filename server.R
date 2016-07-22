library(shiny)
library(CRISPRseek)
library("BSgenome.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(shinyjs)
library(DT)
library(GUIDEseq)

shinyServer(function(input, output) {

options( shiny.maxRequestSize = 1000 * 1024 ^ 2,
           shiny.fullstacktrace = FALSE, shiny.trace=FALSE, 
           shiny.autoreload=TRUE)
  
  output$loading <- renderUI({
    getLoadingMsg()
  })
  output$logo <- renderUI({
    getLogo()
  })
  
  #' getLoadingMsg
  #'
  #' @note \code{getLoadingMsg}
  #' @return loading msg
  #' @examples
  #'     x <- getLoadingMsg()
  #' @export
  #'
  getLoadingMsg <- function() {
    imgsrc <- "images/loading.gif"
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
                                  ))),
        
      conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",
                         tags$div("Analysis complete! Click 'Download Output'", id = "loadmessage" ))
      )
}
  
  #' getLogo
  #'
  #'
  #' @note \code{getLogo}
  #' @return return logo
  #' @examples
  #'     x <- getLogo()
  #' @export
  #'
  getLogo <- function(){
    imgsrc <- "images/logo.png"
    a<-list(img(src=imgsrc, align = "right"))
  }
  
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
  #Mismatch Activity FIle
  if(is.null(input$mismatchActivityFile)) {
    mismatch.activity.file <- system.file("extdata", 
                                          "NatureBiot2016SuppTable19DoenchRoot.csv", 
                                          package = "CRISPRseek")
  }
  else {
    mismatch.activity.file <- input$mismatchActivityFile$datapath
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
  
  isolate(
    if(input$fileFormat == 1) {
      format <- "fasta"
    }
    else if(input$fileFormat == 2) {
      format <- "fastq"
    }
    else if(input$fileFormat == 3) {
      format <- "bed"
    }
  )
  isolate(
    if(input$chooseAction == 2) {
    if(input$fileFormat2 == 1) {
      format <- append(format, "fasta")
    }
    else if(input$fileFormat2 == 2) {
      format <- append(format, "fastq")
    }
    else {
      format <- append(format, "bed")
    }}
  )
  
  isolate(
    if(input$gRNAexport == 1) {
      exportAllgRNAs <- "fasta"
    }
    else if(input$gRNAexport == 2) {
      exportAllgRNAs <- "genbank"
    }
    else if(input$gRNAexport == 3) {
      exportAllgRNAs <- "all"
    }
    else {
      exportAllgRNAs <- "no"
    }
  )
  isolate(
    if(input$fileHeader == 1) {
      header = TRUE
    }
    else {
      header = FALSE
    }
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
  
  #Annotate paired?
  isolate(
    if(input$annPaired1 == 1 || input$annPaired2 == 1) {
      annotatePaired <- TRUE
    }
    else {
      annotatePaired <- FALSE
    }
  )

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
    searchDirection <- append(searchDirection, input$searchDir)
    ) 
    #grna name prefix
    isolate(
      gRNA.name.prefix <- input$gRNAname
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
    isolate(
    if(input$overwriteFile == 1) {
      overwrite <- TRUE
    }
    else {
      overwrite <- FALSE
    })
    #Max mapped R1 bp length to be considered for downstream analysis
    isolate(
    min.R1.mapped <- input$minR1)
    isolate(
    min.R2.mapped <- input$minR2)
    #Specify whether the paired reads are required to align to the same chromosme
    isolate(
    if(input$sameChrom == 1) {
      same.chromosome <- TRUE
    }
    else {
      same.chromosome <- FALSE
    })
    #Does UMI File contain a header?
    isolate(
      if(input$umiheader == 1) {
        umi.header <- TRUE
      }
      else {
        umi.header <- FALSE
      }
    )
    #Index of the column containing the read identified in UMI file
    isolate(
      read.ID.col <- input$readID
    )
    #Index of column containing umi or umi plus the first few bases of sequence from R1 reads
    isolate(
      umi.col <- input$umicol
    )
    #specify whether R1/R2 should be aligned to the same or opposite strands
    isolate(
      if(input$concordantStrand == 1) {
        concordant.strand <- FALSE
      }
      else {
        concordant.strand <- TRUE
      }
    )
    #Maximum distance allowed between R1/R2 reads
    isolate(
      max.paired.distance <- input$maxPairedDistance
    )
    #Minimum mapping quality of acceptable alignments
    isolate(
      min.mapping.quality <- input$minMapQuality
    )
    #Distance value to be assigned to paired reads aligned to different chromosomes
    isolate(
      distance.inter.chrom <- input$distInterChrom
    )
    #apply minimum mapped length requirement to both r1/r2 reads
    isolate(
      if(input$applyMinMapped == 1) {
        apply.both.min.mapped <- TRUE
      }
      else {
        apply.both.min.mapped <- FALSE
      }
    )
    #min number of reads to be considered as a peak
    isolate(
      min.reads <- input$minReads
    )
    #maximum pcalue to be considered as significant
    isolate(
      maxP <- input$max.P
    )
    #Statistical test
    isolate(
      if(input$stat == 1) {
        stats <- "poisson"
      }
      else {
        stats <- "nbinom"
      }
    )
    #Max gap allowed between the plus strand and the negative strand peak
    isolate(
      distance.threshold <- input$distThreshold
    )
    #PAM sequenece after the gRNA
    isolate(
      PAM.pattern <- input$PAMpattern
    )
    #number of degenerative bases in the PAM sequence
    isolate(
      allowed.mismatch.PAM <- input$PAMmismatch
    )
    #weights to be used in SPcas9 system
    v <- as.numeric(unlist(strsplit(input$weight,",")))
    isolate(
      if(length(v) < gRNA.size) {
        x <- (gRNA.size - length(v))
        padZeros <- vector("numeric", length = x)
        weights <- append(padZeros, v)
      }
      else {
        weights <- v
      })
    #min score of an off target to be included in final output
    isolate(
      min.score <- input$minScore)
    #top N off targets to be included in the final output
    isolate(
      topN <- input$top.N)
    #top N off target used to calculate the total off target score
    isolate(
    topN.OfftargetTotalScore <- input$topNscore)
    #Fetch flank sequence of off target or not
    isolate(
      if(input$fetchSeq == 1) {
        fetchSequence <- TRUE
      }
      else {
        fetchSequence <- FALSE
      })
    #upstream offset from the off target start
    isolate(
      upstream <- input$up.stream
    )
    #downstream offset from the off target end
    isolate(
      downstream <- input$down.stream
    )
    # display in gray scale with the darkness indicating the gRNA efficacy
    isolate(
      if(input$usescore == 1) {
        useScore <- TRUE
      }
      else {
        useScore <- FALSE
      })
    
    isolate(
      if(input$efficacyFromIS == 1) {
        useEfficacyFromInputSeq <- TRUE
      }
      else {
        useEfficacyFromInputSeq <- FALSE
      })
    #Indicates which method to use for offtarget cleavage rate estimation
    isolate(
      if(input$scoringmethod == 1) {
        scoring.method <- "Hsu-Zhang"
      }
      else {
        scoring.method <- "CFDscore"
      })
    #Sub PAM activity
    AA = input$AA
    AC = input$AC
    AG = input$AG
    AT = input$AT
    CA = input$CA
    CC = input$CC
    CG = input$CG
    CT = input$CT
    GA = input$GA
    GC = input$GC
    GG = input$GG
    GT = input$GT
    TA = input$TA
    TC = input$TC
    TG = input$TG
    TT = input$TT
    isolate(
      subPAM.activity <- hash::hash("AA" = AA, "AC" = AC, "AG" = AG, "AT" = AT,
                                    "CA" = CA, "CC" = CC, "CG" = CG, "CT" = CT,
                                    "GA" = GA, "GC" = GC, "GG" = GG, "GT" = GT,
                                    "TA" = TA, "TC" = TC, "TG" = TG, "TT" = TT))
    #upstream offset from the bed input starts to search for gRNAs
    isolate(
      upstream.search <- input$upstreamSearch
    )
    #downstream offset from the bed input ends to search for gRNAs
    isolate(
      downstream.search <- input$downstreamSearch
    )
    #The start and end positions of the sub PAM.
    isolate(
      subPAM.position <- c(input$subPamPos1, input$subPamPos2)
    )
    #whether to remove the detailed gRNA information such as efficacy file/restriction enzyme cut sites
    isolate(
      if(input$removeDetails == 1) {
        removegRNADetails <- c(TRUE, TRUE)
      }
      else if(input$removeDetails == 2) {
        removegRNADetails <- c(FALSE, FALSE)
      }
      else if(input$removeDetails == 3) {
        removegRNADetails <- c(TRUE, FALSE)
      }
      else if(input$removeDetails == 4) {
        removegRNADetails <- c(FALSE, TRUE)
      }
    )
    #window size to calculate coverage
    isolate(
      window.size <- input$window
    )
    #step size to calculate coverage
    isolate(
      step <- input$stepSize
    )
    #window size to calculate local background
    isolate(
      bg.window.size <- input$BGWindow
    )
    #Adjustment method for multiple comparisons
    isolate(
      if(input$adjustMethods == 1) {
        p.adjust.methods <- "none"
      }
      else if(input$adjustMethods == 2) {
        p.adjust.methods <- "BH"
      }
      else if(input$adjustMethods == 3) {
        p.adjust.methods <- "holm"
      }
      else if(input$adjustMethods == 4) {
        p.adjust.methods <- "hochberg"
      }
      else if(input$adjustMethods == 5) {
        p.adjust.methods <- "hommel"
      }
      else if(input$adjustMethods == 6) {
        p.adjust.methods <- "bonferroni"
      }
      else if(input$adjustMethods == 7) {
        p.adjust.methods <- "BY"
      }      else if(input$adjustMethods == 8) {
        p.adjust.methods <- "fdr"
      }
    )
   #doesnt work
    # isolate(
     # if(input$PAMlocation == 1) {
      #  PAM.location <- '3prime'
      #}
      #else {
      #  PAM.location <- '5prime'
     # }
   # )
    #Chromosome to exlude searchin in off target analysis
  isolate(
    if(input$chromExclude == 1) {
      chromToExclude <- ""
    }
    else if(input$chromExclude == 2) {
      chromToExclude <- "chr17_ctg5_hap1"
    }
    else if(input$chromExclude == 3) {
      chromToExclude <- "chr4_ctg9_hap1"
    }
    else if(input$chromExclude == 4) {
      chromToExclude <- "chr6_apd_hap1"
    }
    else if(input$chromExclude == 5) {
      chromToExclude <- "chr6_cox_hap2"
    }
    else if(input$chromExclude == 6) {
      chromToExclude <- "chr6_dbb_hap3"
    }
    else if(input$chromExclude == 7) {
      chromToExclude <- "chr6_mann_hap4"
    }
    else if(input$chromExclude == 8) {
      chromToExclude <- "chr6_mcf_hap5"
    }
    else if(input$chromExclude == 9) {
      chromToExclude <- "chr6_qbl_hap6"
    }
    else if(input$chromExclude == 10) {
      chromToExclude <- "chr6_ssto_hap7"
    }
  )
    
    #Run off target analysis
    resultsOTA <- eventReactive(input$goButton, {
      offTargetAnalysis(inputFilePath = inputFilePath, findgRNAsWithREcutOnly=findgRNAsWithREcutOnly,
                        REpatternFile = REpatternFile, findPairedgRNAOnly = findPairedgRNAOnly, 
                        annotatePaired = annotatePaired, BSgenomeName = BSgenomeName, txdb=txdb,
                        orgAnn = orgAnn,max.mismatch = input$mismatch, chromToSearch = chromToSearch,
                        outputDir = outputDir, overwrite = overwrite, allowed.mismatch.PAM = allowed.mismatch.PAM,
                        overlap.gRNA.positions=overlap.gRNA.positions, gRNA.size = gRNA.size,
                        baseBeforegRNA =baseBeforegRNA, baseAfterPAM=baseAfterPAM,
                        max.gap = max.gap, annotateExon = annotateExon, enable.multicore = enable.multicore,
                        PAM.size = PAM.size, temperature = temperature, minREpatternSize = input$REPatSize1,
                        findgRNAs = findgRNAs, PAM = PAM, PAM.pattern = PAM.pattern, format = format,
                        header = header, exportAllgRNAs = exportAllgRNAs, gRNA.name.prefix = gRNA.name.prefix,
                        min.score = min.score, topN = topN, topN.OfftargetTotalScore = topN.OfftargetTotalScore,
                        fetchSequence = fetchSequence, upstream = upstream, downstream = downstream, weights = weights,
                        useScore = useScore, useEfficacyFromInputSeq = useEfficacyFromInputSeq, 
                        scoring.method = scoring.method, subPAM.activity = subPAM.activity,
                        mismatch.activity.file = mismatch.activity.file, upstream.search = upstream.search,
                        downstream.search = downstream.search, subPAM.position = subPAM.position, chromToExclude = chromToExclude)
      
    })
    
    #run compare 2 sequences
    resultsC2S <- eventReactive(input$goButton, {
      compare2Sequences(inputFile1Path = inputFile1Path, inputFile2Path =inputFile2Path, REpatternFile = REpatternFile, outputDir = outputDir, 
                        overwrite = overwrite, BSgenomeName = BSgenomeName, findgRNAsWithREcutOnly=findgRNAsWithREcutOnly, 
                        findPairedgRNAOnly=findPairedgRNAOnly, annotatePaired = annotatePaired,  overlap.gRNA.positions=overlap.gRNA.positions, 
                        gRNA.size = gRNA.size, baseBeforegRNA =baseBeforegRNA, baseAfterPAM=baseAfterPAM, 
                        min.gap = min.gap, max.gap = max.gap, PAM.size = PAM.size, temperature = temperature,
                        minREpatternSize = input$REPatSize2, findgRNAs = findgRNAs, max.mismatch = input$mismatch,
                        searchDirection = searchDirection, PAM = PAM, PAM.pattern = PAM.pattern, allowed.mismatch.PAM = allowed.mismatch.PAM,
                        scoring.method = scoring.method, subPAM.activity = subPAM.activity, weights = weights, 
                        mismatch.activity.file = mismatch.activity.file, subPAM.position = subPAM.position, format = format,
                        header = header, removegRNADetails = removegRNADetails, exportAllgRNAs = exportAllgRNAs, gRNA.name.prefix = gRNA.name.prefix,
                        upstream = upstream, downstream = downstream
                        )
    })
    
    #run GUIDEseqAnalysis
    resultsGSA <- eventReactive(input$goButton, {
      GUIDEseqAnalysis(alignment.inputfile = bamfile, umi.inputfile = umifile,
                       BSgenomeName = BSgenomeName, gRNA.file = gRNA.file,
                       outputDir = outputDir, overlap.gRNA.positions = overlap.gRNA.positions,
                       PAM.size = PAM.size, PAM.pattern = PAM.pattern, gRNA.size = gRNA.size,PAM = PAM, max.mismatch = input$mismatch,
                       overwrite = overwrite, min.R1.mapped = min.R1.mapped, min.R2.mapped = min.R2.mapped,
                       same.chromosome = same.chromosome, umi.header = umi.header, read.ID.col = read.ID.col,
                       umi.col = umi.col, concordant.strand = concordant.strand, max.paired.distance = max.paired.distance,
                       min.mapping.quality = min.mapping.quality, distance.inter.chrom = distance.inter.chrom,
                       apply.both.min.mapped = apply.both.min.mapped, min.reads = min.reads, maxP = maxP, stats = stats,
                       distance.threshold = distance.threshold, weights = weights, window.size = window.size,
                       step = step, bg.window.size = bg.window.size, p.adjust.methods = p.adjust.methods,
                       allowed.mismatch.PAM = allowed.mismatch.PAM, upstream = upstream, downstream = downstream, )})
    
    
    setwd(outputDir)
    
    #Toggle for download button
    disableDownload()
    
    #Run Program
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
        #OTA data table example
        if(input$chooseAction == 1) {
          data <- read.table(paste0(outputDir, "/RECutDetails.xls"),
                             header = TRUE)
        }
        #C2S data table example
        else if(input$chooseAction == 2) {
          data <- read.table(paste0(outputDir, "/scoresFor2InputSequences.xls"),
                             colClasses = c(NA, NA, NA, NA, NA, "NULL", "NULL", "NULL", "NULL", "NULL", 
                                            "NULL", "NULL", "NULL", "NULL"), header = TRUE)
        }
        #GSA data table example
        else if(input$chooseAction == 3) {
          data <- read.table(paste0(outputDir, "/gRNA-peaks.xls"), header = TRUE)
        }
      }
    }))#Data Table
    #Download output as zip file
    output$downloadData <- downloadHandler(
      filename = function() {
        paste0("crisprSeekPlus_", input$runNum,".zip")
      },
      content = function(fname) {
        crisprOutput <- list.files(path = outputDir)
        zip(zipfile = fname, files = crisprOutput)
      },
      contentType = "application/zip"
    )
      return()
  })

#Reset all fields
observeEvent(input$resetFields, {
  reset("fileFormat")
  reset("fileHeader")
  reset("gRNAexport")
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
  reset("sameChrom")
  reset("umiheader")
  reset("readID")
  reset("umicol")
  reset("concordantStrand")
  reset("maxPairedDistance")
  reset("minMapQuality")
  reset("distInterChrom")
  reset("applyMinMapped")
  reset("minReads")
  reset("max.P")
  reset("stat")
  reset("distThreshold")
  reset("PAMpattern")
  reset("weight")
  reset("gRNAname")
  reset("PAMmismatch")
  reset("minScore")
  reset("top.N")
  reset("topNscore")
  reset("fetchSeq")
  reset("up.stream")
  reset("down.stream")
  reset("usescore")
  reset("efficacyFromIS")
  reset("scoringmethod")
  reset("AA")
  reset("AC")
  reset("AG")
  reset("AT")
  reset("CA")
  reset("CC")
  reset("CG")
  reset("CT")
  reset("GA")
  reset("GC")
  reset("GG")
  reset("GT")
  reset("TA")
  reset("TC")
  reset("TG")
  reset("TT")
  reset("upstreamSearch")
  reset("downstreamSearch")
  reset("subPamPos1")
  reset("subPamPos2")
  reset("removeDetails")
  reset("window")
  reset("BGWindow")
  reset("stepSize")
  reset("adjustMethods")
  reset("PAMlocation")
  reset("weight")
  reset("chromExlude")})


##################### VIEW INPUT FILES #####################
  output$fileInput1 <- renderPrint(
    if(is.null(input$file1)) {
      print(read.csv(system.file("extdata", "inputseq.fa", package = "CRISPRseek")))
    }
    else {
      print(read.csv(input$file1$datapath))
    }
  )

  output$fileInput2 <- DT::renderDataTable(DT::datatable({
    if(is.null(input$file2)) {
    data <- read.csv(system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek"), header = TRUE)
    }
    else {
      data <- read.csv(input$file2$datapath)
    }
  }))
  output$mismatchActivityInput <- DT::renderDataTable(DT::datatable({
    if(is.null(input$mismatchActivityFile)) {
      data <- read.csv(system.file("extdata", 
                                   "NatureBiot2016SuppTable19DoenchRoot.csv", 
                                   package = "CRISPRseek"))
    }
    else {
      data <- read.csv(input$mismatchActivityFile$datapath)
    }
  }))
  
  output$fileInput3 <- renderPrint(
    if(is.null(input$file3)) {
      print(read.csv(system.file("extdata", "rs362331T.fa", package = "CRISPRseek")))
    }
    else {
      print(read.csv(input$file3$datapath))
    }
  )
  output$fileInput4 <- renderPrint(
    if(is.null(input$file4)) {
      print(read.csv(system.file("extdata", "rs362331C.fa", package = "CRISPRseek")))
    }
    else {
      print(read.csv(input$file4$datapath))
    }
  )
  
  output$fileInput5 <- DT::renderDataTable(DT::datatable({
    if(is.null(input$file5)) {
          data <- read.csv(system.file("extdata", "UMI-HEK293_site4_chr13.txt", package = "GUIDEseq"))
        }
        else {
          data <- read.csv(input$file5$datapath)
        }
      }))
  output$fileInput7 <- renderPrint(
    if(is.null(input$file7)) {
      print(read.csv(system.file("extdata","gRNA.fa", package = "GUIDEseq")))
    }
    else {
      print(read.csv(input$file7$datapath))
    }
  )

  
})




