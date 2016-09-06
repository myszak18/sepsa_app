library("DESeq2")
library("ggrepel")
library("shiny")
library("biomaRt")
library("GOstats")
library("org.Mm.eg.db")


directory<-"data/"
sampleFiles <- grep(".counts",list.files(directory),value=TRUE)
sampleCondition<-sampleFiles
sampleCondition[grep("K",sampleCondition)]<-"control"
sampleCondition[grep("CLP",sampleCondition)]<-"sepsis"
name_files<-strsplit(sampleFiles,".",fixed=TRUE)
names<-unlist(lapply(name_files,'[[',1))
organ<-strsplit(sampleFiles,"[0-9]{4}",perl=TRUE)
organ<-unlist(lapply(organ,'[[',1))
sampleTable <- data.frame(sampleName = names,
                          fileName = sampleFiles,
                          condition = sampleCondition,
                          organ=organ)


# Define server logic for random distribution application
shinyServer(function(input, output) {
  
  # Reactive expression to generate the requested distribution. This is 
  # called whenever the inputs change. The output expressions defined 
  # below then all used the value computed from this expression
  
  data <- eventReactive(input$goButton,{
    sampleTable<-subset(sampleTable,organ %in% input$organ)
    ddsHTSeq<-NULL
    if (length(input$organ)>1){
    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                           directory = directory,
                                           design= ~ organ+condition)
    }
    else{
    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                             directory = directory,
                                             design= ~ condition)
    }
    dds <- DESeq(ddsHTSeq)
    rld <- rlog(dds)
    data_plot <- plotPCA(rld, intgroup=c("condition","organ"), returnData=TRUE)
    percentVar <- round(100 * attr(data_plot, "percentVar"))
    res <- results(dds)
    resOrdered <- res[order(res$padj),]
    #summary(res)
    resOrdered$gene<-rownames(resOrdered)
    count<-as.data.frame(counts(dds))
    count$gene<-rownames(count)
    resOrdered$FC<-2^resOrdered[["log2FoldChange"]]
    all_results<-merge(count,as.data.frame(resOrdered),by="gene")
    all_results<-all_results[order(all_results$padj,-abs(all_results$log2FoldChange)),]
    ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
    BM_res=getBM(attributes=c("entrezgene","ensembl_gene_id","external_gene_name","gene_biotype"),mart=ensembl)
    all_results<-merge(as.data.frame(all_results),BM_res[,c(2:4)],by.x="gene",by.y="ensembl_gene_id")
    ret_data<-list(data_plot,all_results,percentVar)
    ret_data
  })
  goData<-eventReactive(input$goButton,{
    if (input$GO_type!="none"){
      result<-data.frame(data()[[2]])
      ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
      BM_res=getBM(attributes=c("entrezgene","ensembl_gene_id","external_gene_name","gene_biotype"),mart=ensembl)
      BM_res<-BM_res[BM_res$ensembl_gene_id %in% result[["gene"]],]
      universe<-unique(BM_res[["entrezgene"]])
      up_reg<-result[result$FC < 0.5 & result$padj<0.05,]
      BM_up<-BM_res[BM_res$external_gene_name %in% up_reg[["external_gene_name"]],]
      entrez_up<-unique(BM_up[["entrezgene"]])
      hgCutoff <- 0.99
      params <- new("GOHyperGParams",geneIds=entrez_up,
                    universeGeneIds=universe,
                    annotation="org.Mm.eg.db",
                    ontology=input$GO_type,
                    pvalueCutoff=hgCutoff,
                    conditional=TRUE,
                    testDirection="over")
      hgOver <- hyperGTest(params)
      df_up <- as.data.frame(summary(hgOver))
      df_up$qValue<-p.adjust(df_up[["Pvalue"]],method='fdr')
      df_up$log_odds_ratio<-log(df_up$OddsRatio)
      down_reg<-result[result$FC > 2 & result$padj<0.05,]
      BM_down<-BM_res[BM_res$external_gene_name %in% down_reg[["external_gene_name"]],]
      entrez_down<-unique(BM_down[["entrezgene"]])
      params <- new("GOHyperGParams",
                    geneIds=entrez_down,
                    universeGeneIds=universe,
                    annotation="org.Mm.eg.db",
                    ontology=input$GO_type,
                    pvalueCutoff=hgCutoff,
                    conditional=TRUE,
                    testDirection="over")
      
      hgOver <- hyperGTest(params)
      df_down <- as.data.frame(summary(hgOver))
      df_down$qValue<-p.adjust(df_down[["Pvalue"]],method='fdr')
      df_down$log_odds_ratio<-log(df_down$OddsRatio)
      go_list<-list(df_up,df_down)
      go_list
    }
  })
  clickData<-eventReactive(input$plot_click,{
    dfGO<-nearPoints(data.frame(goData()[[2]])[1:50,], input$plot_click, threshold = 10, maxpoints = 1)
    termGO<-paste("GO",input$GO_type,"ID",sep="")
    choiceGo<-dfGO[[termGO]][1]
    urlGO<-paste("http://www.ebi.ac.uk/QuickGO/GAnnotation?tax=10090&goid=%20",choiceGo,"%20&format=tsv",
                 sep="")
    infoGO<-read.delim(url(urlGO))
    ifnoGO<-infoGO[infoGO$GO.ID==choiceGo,]
    genes<-unique(infoGO[["Symbol"]])
    df_genes_go<-data.frame(data()[[2]])
    df_genes_go<-df_genes_go[df_genes_go$padj<0.05,]
    rel_genes_go<-df_genes_go[df_genes_go$external_gene_name %in% genes,]
    rel_genes_go
  })
  # Generate a plot of the data. Also uses the inputs to build the 
  # plot label. Note that the dependencies on both the inputs and
  # the data reactive expression are both tracked, and all expressions 
  # are called in the sequence implied by the dependency graph
  # Generate an HTML table view of the data
  output$table <- renderDataTable({
    data.frame(data()[[2]])
  },options=list(lengthMenu = c(20, 30, 50), pageLength = 10, orderClasses = TRUE))
  output$downloadData <- downloadHandler(
    filename = function() { paste(input$organ, '.txt', sep='') },
    content = function(file) {
      write.table(data()[[2]], file,sep="\t",row.names=FALSE)
    }
  )
  output$plot<-renderPlot(
    {
      ggplot(data.frame(data()[[1]]), aes(PC1, PC2, color=condition,shape=organ)) +
        geom_point(size=4) +
        xlab(paste0("PC1: ",data()[[3]][1],"% variance")) +
        ylab(paste0("PC2: ",data()[[3]][2],"% variance")) + scale_colour_brewer(palette="Set1") +
        geom_text_repel(aes(label=name), size=3)+scale_shape_manual(values = c(15:18)[1:length(input$organ)])
    }
  )
  output$plot_go_up<-renderPlot(
    {
      ggplot(data.frame(goData()[[1]])[1:50,], aes(x=log_odds_ratio, y=Term, size=Count,color=qValue)) + 
        geom_point()+
        scale_color_gradient2(midpoint=mean(data.frame(goData()[[1]])[["qValue"]][1:50]), low="red", mid="white",
                              high="blue", space ="Lab" )
    }
  )
  output$plot_go_down<-renderPlot(
    {
      ggplot(data.frame(goData()[[2]])[1:50,], aes(x=log_odds_ratio, y=Term, size=Count,color=qValue)) + 
        geom_point()+
        scale_color_gradient2(midpoint=mean(data.frame(goData()[[2]])[["qValue"]][1:50]), low="red", mid="white",
                              high="blue", space ="Lab" )
    }
  )
  output$info<-renderDataTable({
    data.frame(clickData())
  })
})

