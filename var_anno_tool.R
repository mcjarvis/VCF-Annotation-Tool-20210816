###
### 20210809 Tempus Code Challenge: VCF annotation tool
###

setwd("~/Desktop/Projects/tempus_vcf_anno/")

# install required packages and reference genomes from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MutationalPatterns")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

# load required packages
library(BSgenome)
library(MutationalPatterns)
library(ggplot2)
library(gridExtra)
library(VariantAnnotation)
library(httr)
library(shiny)

# load ref genome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

##
## UI Script
##
ui <- fluidPage(
  
  # Inputs
  titlePanel("VCF Annotation Tool"),
  tags$p(style = "font-family:Arial", "Use this Shiny app to annotate a VCF of interest."),
  tags$hr(),
  
  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = "upload", 
                label = "Upload VCF File",
                multiple = F,
                accept = c('.vcf')),
      actionButton(inputId = "go", 
                   label = "Annotate VCF", 
                   icon("arrow-alt-circle-right"), 
                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
      downloadButton("download", "Download annotated VCF file"),
    ),
  
    # Outputs
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Anno. Preview", 
                           tags$hr(),
                           column(12, tableOutput(outputId = "preview"))),
                  tabPanel("Input VCF Stats", 
                           tags$hr(),
                           column(12, tableOutput(outputId = "stats"))),
                  tabPanel("SBSs", 
                           tags$hr(),
                           column(6, plotOutput(outputId = "sbs")),
                           column(6, plotOutput(outputId = "trinuc")))
                
      )
    )
  )
)


##
## Server Script
## 
server <- function(input, output, session) {
  
  ## get expression groups dataframe
  final_table <- eventReactive(input$go, {
    input$upload
    
    vcf_grange <- readVcf(input$upload$datapath, verbose = FALSE, "hg19")
    vcf_df <- as.data.frame(read.table(input$upload$datapath, header = F))
    
    # get INFO column information
    info <- as.data.frame(info(vcf)[ ,c("TYPE", "DP", "AC", "AF")])
    row.names(info) <- NULL
    
    # combine chr/pos columns with INFO columns
    vcf_df1 <- cbind(vcf_df, info)
    
    # get most deletarious alt variant per ref pos
    varRanking <-function(variantType)
    {
      variantWeight <- list("snp"=0, "mnp"=1, "complex"=2, "ins"=3, "del"=3)
      variantClass <- list("snp" = "Single_Nucleotide_Polymorphism", "mnp" = "Multi_Nucleotide_Polymorphism", "complex" = "Complex", "ins" = "Insertion", "del"="Deletion")
      variantWT = vector()
      for (tp in variantType)
      {
        variantWT <- append(variantWT, as.numeric(variantWeight[tp]))
        if (tp == "ins" || tp == "del") return (c(variantType <- ifelse(tp == "ins", variantClass$ins, variantClass$del), which.max(variantWT)))
      }
      if (variantWT[which.max(variantWT)] == 2) return (c(variantType <- variantClass$complex, which.max(variantWT)))
      else if (variantWT[which.max(variantWT)] == 1) return (c(variantType <- variantClass$mnp, which.max(variantWT)))
      else if (variantWT[which.max(variantWT)] == 0) return (c(variantType <- variantClass$snp, which.max(variantWT)))
    }
    
    # get subfinal table of values for vcf
    subfinal_table <- data.frame()
    for (n in 1:nrow(vcf_df1))
    {
      chr <- as.character(data[n,]$V1)
      pos <- as.character(data[n,]$V2)
      ref <- as.character(data[n,]$V4)
      alt <- as.character(data[n,]$V5)
      info <- as.character(data[n,]$V8)
      variantType <- strsplit(as.character(unlist(vcf_df1[n, "TYPE"])), ",")[[1]]
      variannos <- strsplit(as.character(unlist(vcf_df1[n, "AC"])), ",")[[1]]
      variannct <- strsplit(as.character(unlist(vcf_df1[n, "AF"])), ",")[[1]]
      type <- varRanking(variantType)[1]
      depth <- as.numeric(unlist(vcf_df1[n, "DP"]))
      af_count <- variannos[as.numeric(varRanking(variantType)[2])]
      af_percent <- as.numeric(variannct[as.numeric(varRanking(variantType)[2])])*100
      ExAC_id <- paste0(chr,"-",pos,"-",ref,"-",alt)
      subfinal_table <- rbind(subfinal_table, cbind(chr,pos,ref,alt,type,depth,af_count,af_percent,ExAC_id))
    }
    
    # get info from ExAC API
    ExACretrieve <- httr::POST(url="http://exac.hms.harvard.edu/rest/bulk/variant", body=jsonlite::toJSON(as.character(subfinal_table$ExAC_id)), encode = "json")
    ExAC_API_data <- content(ExACretrieve)
    
    ExAC_Results <- vector()
    for (n in 1:nrow(subfinal_table))
    {
      ExACinfo <- jsonContent[[as.character(subfinal_table$ExAC_id[n])]]
      exACall <- as.character(subfinal_table$ExAC_id[n])
      alleleFreq <- if (!(is.null(ExACinfo$variant$allele_freq))) paste(unlist(ExACinfo$variant$allele_freq),collapse = ",") else "NA"
      varConseq <- if (!(is.null(names(ExACinfo$consequence)))) paste(unlist(names(ExACinfo$consequence)), collapse = ",") else "NA"
      rsid <- if (!(is.null(ExACinfo$variant$rsid))) paste(unlist(ExACinfo$variant$rsid), collapse = ",") else "NA"
      ensgid <- if (!(is.null(ExACinfo$variant$genes))) paste(unlist(ExACinfo$variant$genes), collapse = ",") else "NA"
      enstid <- if (!(is.null(ExACinfo$variant$transcripts))) paste(unlist(ExACinfo$variant$transcripts), collapse = ",") else "NA"
      ExAC_Results <- rbind(ExAC_Results, cbind(exACall, alleleFreq, varConseq, rsid, ensgid, enstid))
    }
    colnames(ExAC_Results) <- c("ExAC_id", "Allele_Frequency","Variant_Consequences", "Reference_SNP_ID", "Ensembl_Gene_ID", "Ensembl_Transcript_ID")
    
    # get final table
    final_table <- merge(subfinal_table, ExAC_Results, by="ExAC_id")
    final_table <- final_table[ ,c(2:14,1)]
    final_table
  })
  
  stats_table <- eventReactive(input$go, {
    input$upload
    table1 <- as.data.frame(table(final_table()[,"type"]))
    colnames(table1) <- c("Variant Type", "Number of Variants")
    table2 <- as.data.frame(aggregate(as.numeric(depth) ~ as.factor(type), final_table(), mean))
    colnames(table2) <- c("Variant Type", "Mean Read Depth")
    table3 <- merge(table1, table2, by = "Variant Type")
    table3
  })
  
  sbs_profile <- eventReactive(input$go, {
    input$upload
    vcfs <- read_vcfs_as_granges(input$upload$datapath, "sample", ref_genome)
    type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
    sbs_plot <- plot_spectrum(type_occurrences, CT = F, legend = F)
    sbs_plot
  })
  
  sbs_trinuc <- eventReactive(input$go, {
    input$upload
    vcfs <- read_vcfs_as_granges(input$upload$datapath, "sample", ref_genome)
    mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
    trinuc_plot <- plot_96_profile(mut_mat, ymax = 0.1)
    trinuc_plot
  })
  
  # render objects
  output$preview <- renderTable({
    head(final_table(), 5)},
    align = "c"
  )
  
  output$stats <- renderTable({
    head(stats_table(), 5)},
    align = "c"
  )
  
  output$sbs <- renderPlot({
    sbs_profile()
  })
  
  output$trinuc <- renderPlot({
    sbs_trinuc()
  })
  
  output$download <- downloadHandler(
    filename = function() {paste0("output.csv")},
    content = function(file) {
      write.csv(final_table(), file)
    }
  )
}


##
## Run the shiny app
## 
shinyApp(ui = ui, server = server)
