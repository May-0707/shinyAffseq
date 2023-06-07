library(AnnotationHub)
library(org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)
library(Rgraphviz)
library(enrichplot)
library(cowplot)

library(GSEABase)

observe({
  
  observeEvent(input$runGSEA, {
  if(input$Tabs_gseaset=="(human) hallmark gene sets"){
    gmtfile ='geneSets/h.all.v2022.1.Hs.symbols.gmt'
  }else if(input$Tabs_gseaset=="(human) C1 positional gene sets"){
    gmtfile ='geneSets/c1.all.v2022.1.Hs.symbols.gmt'
  }else if(input$Tabs_gseaset=="(human) C2 curated gene sets"){
    gmtfile ='geneSets/c2.all.v2022.1.Hs.symbols.gmt'
  }else if(input$Tabs_gseaset=="(human) C3 regulatory target gene sets"){
    gmtfile ='geneSets/c3.all.v2022.1.Hs.symbols.gmt'
  }else if(input$Tabs_gseaset=="(human) C4 computational gene sets"){
    gmtfile ='geneSets/c4.all.v2022.1.Hs.symbols.gmt'
  }else if(input$Tabs_gseaset=="(human) C5 ontology gene sets"){
    gmtfile ='geneSets/c5.all.v2022.1.Hs.symbols.gmt'
  }else if(input$Tabs_gseaset=="(human) C6 oncogenic signature gene sets"){
    gmtfile ='geneSets/c6.all.v2022.1.Hs.symbols.gmt'
  }else if(input$Tabs_gseaset=="(human) C7 immunologic signature gene sets"){
    gmtfile ='geneSets/c7.all.v2022.1.Hs.symbols.gmt'
  }else if(input$Tabs_gseaset=="(human) C8 cell type signature gene sets"){
    gmtfile ='geneSets/c8.all.v2022.1.Hs.symbols.gmt'
  }else if(input$Tabs_gseaset=="(mouse) hallmark gene sets"){
    gmtfile ='geneSets/mh.all.v2022.1.Mm.symbols.gmt'
  }else if(input$Tabs_gseaset=="(mouse) M1: positional gene sets"){
    gmtfile ='geneSets/m1.all.v2022.1.Mm.symbols.gmt'
  }else if(input$Tabs_gseaset=="(mouse) M2: curated gene sets"){
    gmtfile ='geneSets/m2.all.v2022.1.Mm.symbols.gmt'
  }else if(input$Tabs_gseaset=="(mouse) M3: regulatory target gene sets"){
    gmtfile ='geneSets/m3.all.v2022.1.Mm.symbols.gmt'
  }else if(input$Tabs_gseaset=="(mouse) M5: ontology gene sets"){
    gmtfile ='geneSets/m5.all.v2022.1.Mm.symbols.gmt'
  }else if(input$Tabs_gseaset=="(mouse) M8: cell type signature gene sets"){
    gmtfile ='geneSets/m8.all.v2022.1.Mm.symbols.gmt'
  }
  
  
  geneset <- read.gmt(gmtfile)

  geneFCdataMat <- input$geneFCdata
  dataMat <- read.table(geneFCdataMat$datapath, header = T,sep = "\t")
  
  colnames(dataMat) = c("genes","avg_logFC")
  dataMat = dataMat[!duplicated(dataMat$genes),]
  rownames(dataMat) = dataMat$genes
  resOrdered = dataMat[order(-dataMat$avg_logFC),]
  geneList= resOrdered$avg_logFC 
  names(geneList)= toupper(rownames(resOrdered))
  geneList=sort(geneList,decreasing = T)
  
  #output$testA <- renderDataTable({datatable({resOrdered})})
  
  egmt <- GSEA(geneList, TERM2GENE=geneset, 
               minGSSize = 1,
               pvalueCutoff = 0.99,
               verbose=FALSE)
  output$gseaplot <- renderPlot({
    print(gseaplot2(egmt,geneSetID = 1, pvalue_table=T))
    
  })
  
  output$DownloadgseaTab <- renderUI({
    downloadButton('downloadDataGSEA', 'Download Table')
    
  })
  output$DownloadgseaPicTab <- renderUI({
    downloadButton('downloadGSEAPic', 'Download GSEA Picture')
  })
  output$downloadDataGSEA <- downloadHandler(
    filename = 'GSEA.xls',
    content = function(file){
      write.table(egmt@result, file,sep = "\t")
    })
  gseapic1 = gseaplot2(egmt,geneSetID = 1, pvalue_table=T)
  output$downloadGSEAPic <- downloadHandler(
    filename = 'GSEA.pdf',
    content = function(file){
      ggsave(file,plot = ggpubr::ggarrange(gseapic1[[1]],gseapic1[[2]],gseapic1[[3]],nrow = 3, ncol = 1))
    })
  
  ### gseKEGG
  organ = str_extract(input$Tabs_gseaset,"(?<=\\().+?(?=\\))")
  if(organ=="human"){
    gsek_orga = "hsa"
    selOrgDb = "org.Hs.eg.db"
  }else if(organ=="mouse"){
    gsek_orga = "mmu"
    selOrgDb = "org.Mm.eg.db"
  }
  # change SYMBOL to ENTREZID
  data.df <-data.frame(SYMBOL=rownames(resOrdered),logFC=resOrdered$avg_logFC)
  gene.df <-bitr(data.df$SYMBOL,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb=selOrgDb)
  data.df <- data.df[which(gene.df$SYMBOL %in% data.df$SYMBOL),]
  data.df_new <-merge(data.df,gene.df,by.y="SYMBOL")
  data.df_newGenelist = data.df_new$logFC
  names(data.df_newGenelist)=data.df_new$ENTREZID
  data.df_newGenelist=sort(data.df_newGenelist,decreasing=T)
  
  KEGG_gseresult <- gseKEGG(data.df_newGenelist,organism = gsek_orga, minGSSize = 3, maxGSSize = 500, pvalueCutoff=0.9,eps=0)
  KEGG_gseresultOrder = KEGG_gseresult[order(KEGG_gseresult$enrichmentScore,decreasing=T)]
  #write.csv(KEGG_gseresultOrder, file ="gsedownDEGs_KEGGresult.csv", row.names =TRUE)
  
  ###分面点图激活和抑制
  gsekeggOut = dotplot(KEGG_gseresult,showCategory = 15)#,split=".sign"+facet_wrap(~.sign,scales="free")
  output$gseKEGGout <- renderPlot({
    print(gsekeggOut)
  })
  output$DownloadgsekeggTab <- renderUI({
    downloadButton('downloadDataGSEAkegg', 'Download gseKEGG Table')
  })
  output$downloadDataGSEAkegg <- downloadHandler(
    filename = 'gseKEGG.xls',
    content = function(file){
      write.table(KEGG_gseresultOrder, file,sep = "\t")
    })
  output$DownloadgsekeggPicTab <- renderUI({
    downloadButton('downloadgseakeggPic', 'Download gseKEGG Picture')
  })
  output$downloadgseakeggPic <- downloadHandler(
    filename = 'gseKEGG.pdf',
    content = function(file){
      ggsave(file,plot = gsekeggOut)
    })
  
  
 })## end observeEvent
})## end observe



