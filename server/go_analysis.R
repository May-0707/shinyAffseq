library(AnnotationHub)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(clusterProfiler)
library(Rgraphviz)
library(enrichplot)
library(cowplot)
library(shinyjs)



### select species
observe({
output$species <- renderUI({
  selectInput('selspecies', 'select species',
              choices = c("human","mouse"), 
              selected = NULL, multiple = FALSE, selectize = FALSE
              )})
observeEvent(input$runGO, {
  withProgress(message = 'Calculation in progress',
               #style = shinyOptions(progress.style="notification",'position: fixed; top: 50%; left: 50%; transform: translate(-50%, -50%);'),
               detail = 'This may take a while...', value = 0, {
                 tags$div(
                   style = 'position: fixed; top: 50%; left: 50%; transform: translate(-50%, -50%);',
                   tags$div(
                     style = 'display: inline-block; margin-top: -50px; margin-left: -50px; position: absolute; top: 50%; left: 50%;',
                     tags$div(
                       style = 'width: 100px; height: 100px; border: 10px solid #f3f3f3; border-top: 10px solid #3498db; border-radius: 50%; animation: spin 2s linear infinite;',
                     ),
                     tags$h3('This may take a while...')
                   )
                 )
  selOrgDb = ifelse(input$selspecies=="human", "org.Hs.eg.db","org.Mm.eg.db")
  #selOrgDb = "org.Hs.eg.db"
  if(isTruthy(input$Gene4Golist)){
    genelist = input$Gene4Golist #as.list(str_split(input$Gene4Golist, "\r\n"))

    name_ID = bitr(as.list(str_split(input$Gene4Golist,"\n") %>% unlist %>% str_trim()), fromType="SYMBOL", toType = c("ENTREZID"), OrgDb=selOrgDb)
  }else{
    inFile <- input$upload
    genelist <- read.table(inFile$datapath, header = F,sep = "\t")
    name_ID = bitr(as.list(genelist$V1), fromType="SYMBOL", toType = c("ENTREZID"), OrgDb=selOrgDb)
  }
  
  #output$test <- renderText({str_split(input$Gene4Golist, "\r")})
  #output$test2 <- renderText({as.list(str_split(input$test,"\n") %>% unlist %>% str_trim()) })
  #name_ID = bitr(as.list(genelist$V1), fromType="SYMBOL", toType = c("ENTREZID"), OrgDb=selOrgDb)
  ##output$test <- renderDataTable({datatable({name_ID},options = list(dom = "f",ordering = F),rownames = T)})
  ego_BP <- enrichGO(gene = name_ID$ENTREZID,
                     OrgDb=selOrgDb,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 3,
                     pvalueCutoff = 0.1,
                     qvalueCutoff = 0.1,
                     readable = TRUE)
  
  
  bp.filter<- clusterProfiler::simplify(ego_BP,cutoff=0.7,by="p.adjust",select_fun=min)
  bp.dot <-dotplot(bp.filter,showCategory=30,color="pvalue",title="EnrichmentGO_BP")
  bp <- bp.dot  + scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=35))
  }) ##wait
  output$BP <- renderPlot({
    res = 300
    print(bp)
  })
  
  ego_CC <- enrichGO(gene = name_ID$ENTREZID,
                     OrgDb=selOrgDb,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 3,
                     pvalueCutoff = 0.1,
                     qvalueCutoff = 0.1,
                     readable = TRUE)
  CC.filter<- clusterProfiler::simplify(ego_CC,cutoff=0.7,by="p.adjust",select_fun=min)
  CC.dot <-dotplot(CC.filter,showCategory=30,color="pvalue",title="EnrichmentGO_CC")
  CC <- CC.dot  + scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=35))
  #CC<-plot_grid(CC.barNew,CC.dotNew,ncol=2)
  
  output$CC <- renderPlot({
    res = 300
    print(CC)
  })
  
  ego_MF <- enrichGO(gene = name_ID$ENTREZID,
                     OrgDb=selOrgDb,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     minGSSize = 3,
                     pvalueCutoff = 0.1,
                     qvalueCutoff = 0.1,
                     readable = TRUE)
  MF.filter<- clusterProfiler::simplify(ego_MF,cutoff=0.7,by="p.adjust",select_fun=min)
  #MF.bar <-barplot(MF.filter,showCategory=30,color="pvalue")
  #MF.barNew = MF.bar  + scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=35))
  MF.dot <-dotplot(MF.filter,showCategory=30,color="pvalue",title="EnrichmentGO_MF")
  MF <- MF.dot  + scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=35))
  #MF<-plot_grid(MF.barNew,MF.dotNew,ncol=2)
  #ggsave("ego_MF_bardot.pdf",MF,width=15,height=13)
  #write.csv(MF.filter,file = "ego_MF.csv")
  
  output$MF <- renderPlot({
    res = 300
    print(MF)
  })
  
  output$DownloadGOTab <- renderUI({
    downloadButton('downloadDataBPCCMF', 'Download GO Table')
  })
  output$DownloadGOPic <- renderUI({
    downloadButton("downloadDataPic", "Download GO Picture")
  })
  
  
  
  GO_table = as.data.frame(rbind(bp.filter@result, CC.filter@result, MF.filter@result))
  
  ### download bp cc mf table
  output$downloadDataBPCCMF <- downloadHandler(
    filename = 'BPCCMF.csv',
    content = function(file) {
      write.csv(GO_table, file)
    })
  
  #GO_Pic = plot_grid(bp,CC,MF,ncol=3)
  GO_Pic = ggpubr::ggarrange(bp,CC,MF,nrow = 1, ncol = 3, labels = c('BP', 'CC', 'MF'))#, widths=c(3,1)
  ### download bp cc mf picture
  output$downloadDataPic <- downloadHandler(
    filename = 'BPCCMF.pdf',
    content = function(file){
      ggsave(file,plot = GO_Pic,width =840, height = 297, units="mm")
    })
  
  ### KEGG
  kk <- enrichKEGG(gene = name_ID$ENTREZID,
                   organism ="human",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   minGSSize = 1,
                   use_internal_data =FALSE)
  sig.kegg<-filter(kk,pvalue<0.05)
  kegg.dot<-dotplot(sig.kegg,showCategory=30,color ="pvalue",title="Enrichment_KEGG")
  
  
  output$KEGGout <- renderPlot({
    res = 300
    print(kegg.dot)
    #ggplot(kk,aes(x=-log(p.adjust),y=reorder(Description,-log(p.adjust)),fill=GeneRatio))+geom_bar(stat="identity")
    
  })
  
  ### download KEGG table
  output$downloadDataKEGGTab <- renderUI({
    downloadButton("downloadDataKEGG", "Download KEGG Table")
  })
  
  kegg_table = as.data.frame(kk@result)
  output$downloadDataKEGG <- downloadHandler(
    filename = 'KEGG.csv',
    content = function(file) {
      write.csv(kegg_table, file)
    })
  output$downloadKEGGPicTab <- renderUI({
    downloadButton("downloadKEGGPic", "Download KEGG Picture")
  })
  output$downloadKEGGPic <- downloadHandler(
    filename = 'KEGG.pdf',
    content = function(file){
      ggsave(file,plot = kegg.dot)
    })
  
})

})




