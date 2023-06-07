library(SummarizedExperiment)
#devtools::install_github (repo = "BioinformaticsFMRP/TCGAbiolinks")
#library(TCGAbiolinks)
library(survival)
library(survminer)

#library("RTCGA.clinical")
#library("RTCGA.mRNA")
library("tidyverse")
library("dplyr")
observe({
  observeEvent(input$runOS, {
    
    ## ensemblID to gene symbol
    geneprobe = read.table("gencode.v22.annotation.gene.probeMap",sep = "\t",header = T)
    
    # Create checkData object
    OSgene <- input$OSsearchGene
    cancerType = input$Tabs_OScancer
    # load cancer Rdata
    load(paste0("Rdata/TCGA-",cancerType,".htseq_fpkm-uq.Rdata"))
    data_fpclin <- data_fpclin
    
    if (OSgene %in% colnames(data_fpclin)) {
      checkData <- data_fpclin[, c("sample", "OS", "OS.time", OSgene)]
      checkData$group <- ifelse(checkData[, OSgene] > apply(checkData[OSgene], 2, median), "high", "low")
      
      # Render data table
      output$testB <- renderDataTable({
        datatable(checkData)
      })
      
      # Create survival plot
      dd_sfit <- survfit(Surv(OS.time, OS) ~ group, data = checkData)
      OSout <- ggsurvplot(dd_sfit, data = checkData, conf.int =F, pval = TRUE, risk.table = TRUE, title = "Overall survival", palette =c('red','black'), xlab = "Time (Days)")
      output$OSplot <- renderPlot({
        print(OSout)
        #ggsurvplot(dd_sfit, data = checkData, conf.int =F, pval = TRUE, risk.table = TRUE, title = "Overall survival", palette =c('red','black'))
      })
      output$DownloadOSplotTab <- renderUI({
        downloadButton('downloadOSplot', 'Download overall survival')
      })
      output$downloadOSplot <- downloadHandler(
        filename = 'OverallSurvival.pdf',
        content = function(file){
          pdf(file)
          print(OSout, newpage=FALSE)
          dev.off()
        })
      
      # Create subset of data for 5-year survival plot
      if (max(checkData$OS.time) > 1501) {
        checkData5Year <- checkData[which(checkData$OS.time <= 1500), ]
        checkData5Year$group <- ifelse(checkData5Year[, OSgene] > apply(checkData5Year[OSgene], 2, median), "high", "low")
        dd5_sfit <- survfit(Surv(OS.time, OS) ~ group, data = checkData5Year)
        OS5out <- ggsurvplot(dd5_sfit, data = checkData5Year, conf.int =F, pval = TRUE, risk.table = TRUE, title = "5-year survival", palette =c('red','black'), xlab = "Time (Days)")
        output$OSplot5year <- renderPlot({
          print(OS5out)
        })
        output$DownloadOS5plotTab <- renderUI({
          downloadButton('downloadOS5plot', 'Download 5-year survival')
        })
        output$downloadOS5plot <- downloadHandler(
          filename = '5-year-Survivalplot.pdf',
          content = function(file){
            pdf(file)
            print(OS5out, newpage=FALSE)
            dev.off()
          })
        
      }
      
      # Create subset of data for 10-year survival plot
      if(max(checkData$OS.time) > 3001){
        checkData10Year = checkData[which(checkData$OS.time<=3000),]
        checkData10Year$group = ifelse(checkData10Year[,OSgene] > apply(checkData10Year[OSgene],2,median),"high","low")
        dd10_sfit <- survfit(Surv(OS.time, OS) ~ group, data = checkData10Year)
        OS10out <- ggsurvplot(dd10_sfit, data = checkData10Year, conf.int =F, pval = TRUE, risk.table = TRUE, title = "10-year survival", palette =c('red','black'), xlab = "Time (Days)")
        
        output$OSplot10year <- renderPlot({
          print(OS10out)
        })
        output$DownloadOS10plotTab <- renderUI({
          downloadButton('downloadOS10plot', 'Download 10-year survival')
        })
        output$downloadOS10plot <- downloadHandler(
          filename = '10-year-Survivalplot.pdf',
          content = function(file){
            pdf(file)
            print(OS10out, newpage=FALSE)
            dev.off()
          })
      }

    } else {
      showModal(modalDialog(
        title = "Error",
        "The gene symbol you queried does not exist in the matrix.\n
    Please enter an alias or other gene symbol!",
    easyClose = TRUE,
    footer = NULL
      ))
    }


  }) ##end of observeEvent
})
