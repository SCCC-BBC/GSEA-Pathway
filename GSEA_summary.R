#!/usr/bin/env Rscript
## Collect arguments
args = commandArgs(trailingOnly=T)

## Default setting when no arguments passed
if(length(args) < 1) {
  args = c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      GSEA result summary 
 
      Arguments:
      --indir=somevalue   - character, input directory containing all GSEA output folders 
      --outdir=indir  - character, output directory
      --FDR_cut=0.25  - number between 0 and 1
      --pval_cut=0.05 - number between 0 and 1, either FDR_cut or pval_cut should be specified
      --help              - print this text
 
      Example:
      Rscript GSEA_summary.R --indir=/home/jamesban/Desktop/Projects/Yuqi/GSEA_platform_hs/GSEA

      Output:
	    single GSEA results summary file
      \n")
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs = function(x) strsplit(sub("^--", "", x), "=")
argsDF = as.data.frame(do.call("rbind", parseArgs(args)))
argsL = as.list(as.character(argsDF$V2))
names(argsL) = argsDF$V1

# setup arguments and folders
if(is.null(argsL$indir)) {
  stop("Use --help to view usage example", call.=F)
}

if(is.null(argsL$outdir)) {
  cat("Default output directory: indir\n")
  out_dir = argsL$indir
}else{
  out_dir = argsL$outdir
}

if(is.null(argsL$pval_cut)&is.null(argsL$FDR_cut)){
  cat("Default cutoff use FDR 0.25\n", call.=F)
  order_by = "FDR"
  Cutoff = 0.25
}else if(!is.null(argsL$FDR_cut)) {
  order_by = "FDR"
  Cutoff = as.numeric(argsL$FDR_cut)
}else{
  order_by = "pvalue"
  Cutoff = as.numeric(argsL$pval_cut)
}

require(gplots)
require(RColorBrewer)
require(gtable)
require(grid)
require(cowplot)
require(igraph)

GSEA_summarize = function(GSEA_dir, output_dir, Cutoff = 0.25){
  File_gsea = dir(GSEA_dir)
  if(length(grep("Preranked", File_gsea))>0){
    Rank = c("gsea_report_for_na_pos_", "gsea_report_for_na_neg_")
  }else{
    Rank = c("gsea_report_for_high_", "gsea_report_for_low_")
  }
  cat("NAME\tCOMPARISON\tDIRECTION\tNES\tPVAL\tFDR\tSIZE\tMEMBER_overlap\n",file=paste0(output_dir, "/summary.gsea.txt"),append=F)
  for(j in 1:length(File_gsea)){
    input_dir = paste0(GSEA_dir, "/", File_gsea[j])
    if(length(dir(input_dir))==0){
      next
    }
    data.files = c(grep(".xls", grep(Rank[1], dir(input_dir), value=T), value = T), 
                   grep(".xls", grep(Rank[2], dir(input_dir), value=T), value = T))
    # get member genes of pathways
    gmt.dir = paste0(input_dir, "/edb/gene_sets.gmt")
    data.gmt = readLines(gmt.dir)
    pathwayName = sapply(strsplit(data.gmt, "\\t"), `[`, 1)
    pathway2gene = sapply(strsplit(data.gmt, "\\t"), `[`, -c(1,2))
    pathway2gene = unlist(lapply(pathway2gene, function(h)paste(h, collapse = ",")))
    for(i in 1:length(Rank)){
      if(i==1){
        reg_ty = "UP"
      }else{
        reg_ty = "DN"
      }
      file.dir = paste0(input_dir, "/", data.files[i])
      Temp = read.table(file.dir, header=T, sep="\t")
      if(order_by == "pvalue"){
        Temp = Temp[order(Temp$NOM.p.val, Temp$FDR.q.val), ]
        Temp = Temp[Temp$NOM.p.val<Cutoff, ]
      }else{
        Temp = Temp[order(Temp$FDR.q.val, Temp$NOM.p.val), ]
        Temp = Temp[Temp$FDR.q.val<Cutoff, ]
      }
      if(nrow(Temp)==0){
        next
      }
      Temp1 = data.frame(NAME = Temp$NAME, COMPARISON = unlist(strsplit(File_gsea[j], split=".Gsea."))[1], DIRECTION = reg_ty, NES = Temp$NES, PVAL = Temp$NOM.p.val, 
                         FDR = Temp$FDR.q.val, SIZE = Temp$SIZE, MEMBER_overlap = pathway2gene[match(as.character(Temp$NAME), pathwayName)])
      write.table(Temp1, file=paste0(output_dir, "/summary.gsea.txt"), append=T, row.name=F, col.name=F, sep="\t", quote=F)
    }
  }
}

GSEA_summarize(argsL$indir, out_dir, Cutoff = Cutoff)


# pathway heatmap
# input from GSEA summary file
Data = read.table(paste0(out_dir, "/summary.gsea.txt"), sep="\t",header=T)

if(any(grepl("_BP", Data$COMPARISON))){
  Data$COMPARISON = gsub("_BP", "", Data$COMPARISON)
}
if(any(grepl("_CP", Data$COMPARISON))){
  Data$COMPARISON = gsub("_CP", "", Data$COMPARISON)
}

pathway.unique = unique(Data$NAME)
comp.unique = unique(Data$COMPARISON)
cat("Choose column order for heatmap: \n")
for(i in 1:length(comp.unique)){
  cat(i,"---",comp.unique[i],"\n")
}
input.con = file("stdin")
x <- readLines(input.con,1)
close(input.con)
order.comp = rep(0, length(comp.unique))
for(i in 1:length(comp.unique)){
  order.comp[i] = substr(x, i, i)
}
comp.unique = comp.unique[as.numeric(order.comp)]
print(comp.unique)
#comp.unique = c("ctrl.1vsmupa.1", "ctrl.2vsmupa.2", "ctrl.5vsmupa.5", "ctrl.2vshupa.2", 
#                "ctrl.5vshupa.5",  "mupa.2vshupa.2", "mupa.5vshupa.5")

if(length(comp.unique)<=1){
  cat("There is no time course sample, skip heatmap for pathway...\n")
}else{
  Data.heatmap = data.frame(matrix(0, nrow=length(pathway.unique), ncol=length(comp.unique)))
  for(i in 1:nrow(Data.heatmap)){
    Temp.match.id = which(Data$NAME==pathway.unique[i])
  
    for(j in 1:length(Temp.match.id)){
      Temp.comp.id = match(Data$COMPARISON[Temp.match.id[j]], comp.unique)
      if(Data$DIRECTION[Temp.match.id[j]] == "DN"){
        if(Data$PVAL[Temp.match.id[j]]<0.01){
          Score = -2
        }else{
          Score = -1
        }
      }else{
        if(Data$PVAL[Temp.match.id[j]]<0.01){
          Score = 2
        }else{
          Score = 1
        }
      }
      Data.heatmap[i, Temp.comp.id] = Score
    }
  }
  colnames(Data.heatmap) = comp.unique
  rownames(Data.heatmap) = pathway.unique
  
  hmcol<-rev(colorRampPalette(brewer.pal(10, "RdBu"))(5))
  lmat <- rbind( c(5,3,4), c(2,1,4) )
  lhei <- c(1, 6)
  lwid <- c(1, 6, 1)
  
  pdf(file=paste0(out_dir,"/pathway_heatmap.pdf"), width=8, height=16)
  par(mar=c(10,0,4,4))
  heatmap.2(as.matrix(Data.heatmap), col=hmcol, Colv = FALSE, dendrogram ="none",
            scale="none", trace="none", key = F, 
            lmat = lmat, lhei = lhei, lwid = lwid,
            # draw grid or not
            #sepwidth=c(0,0),sepcolor="black",colsep=1:ncol(Data.heatmap),rowsep=1:nrow(Data.heatmap),
            cexRow = 0.2 + 1/log10(nrow(Data.heatmap)), margins=c(8,16), srtCol=45,
            main="Gene set enrichment map")
  par(lend=2)           # square line ends for the color legend
  legend("topright",      # location of the legend on the heatmap plot
        legend = c("Upregulated", "Weak upregulated", "No differential", "Weak downregulated", "Downregulated"), # category  labels
        col = rev(hmcol),  # color key
        lty= 1,             # line style
        lwd = 15            # line width
  )
  dev.off()
}
  
# plot barchar separately
barchart_gsea = function(input_dir, out_file, nTop = 15, cutoff = 0.25, method = "FDR"){
  if(length(grep("Preranked", input_dir))>0){
    Rank = c("gsea_report_for_na_pos_", "gsea_report_for_na_neg_")
  }else{
    Rank = c("gsea_report_for_high_", "gsea_report_for_low_")
  }
  data.files = c(grep(".xls", grep(Rank[1], dir(input_dir), value=T), value = T), 
                 grep(".xls", grep(Rank[2], dir(input_dir), value=T), value = T))
  
  for(i in 1:2){
    if(i==1){
      title_highlow = "Pathways enriched by up-regulated genes       "
      outfile_ext = "_upreg.png"
    }else{
      title_highlow = "Pathways enriched by down-regulated genes       "
      outfile_ext = "_dnreg.png"
    }
    file.dir = paste0(input_dir, "/", data.files[i])
    Temp = read.table(file.dir, header=T, sep="\t")
    if(method == "FDR"){
      Temp = Temp[order(Temp$FDR.q.val), ]
      Temp$Score = ifelse(-log(Temp$FDR.q.val)>=10, exp(-10), Temp$FDR.q.val)
      xLABEL = "-log(FDR.q.val)"
    }else{
      # use nominal p value
      Temp = Temp[order(Temp$NOM.p.val), ]
      Temp$Score = ifelse(-log(Temp$NOM.p.val)>=10, exp(-10), Temp$NOM.p.val)      
      xLABEL = "-log(p.value)"
    }
    Temp$Score = -log10(Temp$Score)
    Temp$Order = nrow(Temp):1
    Temp.NAME = rep(0, nrow(Temp))
    Temp$NAME = as.character(Temp$NAME)
    for(i in 1:nrow(Temp)){
      if(nchar(Temp$NAME[i])>70){
        Temp.NAME[i] = paste0(substr(Temp$NAME[i], start=1, stop=70), "...")
      }else{
        Temp.NAME[i] = as.character(Temp$NAME[i])
      }
      if(nchar(Temp.NAME[i])>42){
        Temp.NAME[i] = paste0(substr(Temp.NAME[i], start=1, stop=42), 
                              "-\n-", substr(Temp.NAME[i], start=43, stop=nchar(Temp.NAME[i])))
      }
    }
    
    Temp$NAME = Temp.NAME
    Temp$NAME = factor(Temp$NAME, levels = as.character(Temp$NAME[nrow(Temp):1]))
    
    #theme_set(theme_cowplot(font_size=12))
    CalcFudgeAxis = function( y1, y2=y1) {
      Cast2To1 = function(x) (ylim1[1]+(ylim1[2]-ylim1[1])/(ylim2[2]-ylim2[1])*(x-ylim2[1])) # x gets mapped to range of ylim1
      ylim1 <- c(0,max(y1))
      y2 = c(y2, 230)
      ylim2 <- c(0,max(y2))    
      yf <- Cast2To1(y2)
      labelsyf <- pretty(y2)  
      return(list(
        yf=yf[-length(yf)],
        labels=labelsyf[-length(labelsyf)],
        breaks=Cast2To1(labelsyf[-length(labelsyf)])
      ))
    }
    if(nrow(Temp)<nTop){
      Temp.selected = Temp
      nTop = nrow(Temp)
    }else{
      Temp.selected = Temp[1:nTop, ]
    }
    
    y1_max = 3
    Temp.FudgeAxis <- CalcFudgeAxis(seq(length=10, from=0, to=y1_max), Temp.selected$SIZE)
    
    p1 = ggplot(Temp.selected) + 
      geom_bar(aes(y=Score, x=NAME), stat = "identity", width=0.7, color="grey30", fill="grey30") +
      #geom_text(aes(y=0, x=Order, label = NAME), hjust=1, vjust=0, color="black", size=3) + 
      coord_flip(ylim=c(0,y1_max)) + xlab("") + ylab(xLABEL) + theme_bw() +
      theme(axis.ticks.y = element_blank(), axis.text.y = element_text(color="black", size=10), plot.margin=unit(c(0,0,0,0),"mm")) + 
      geom_point(aes(y= Temp.FudgeAxis$yf, x=NAME), color="orange") + 
      #geom_line(aes(y=Temp.FudgeAxis$yf, x=NAME), color="orange", alpha=0.5) + 
      geom_hline(yintercept = -log10(cutoff), color="red", linetype = "dashed") + 
      geom_text(aes( y=-log10(cutoff), x=NAME[nTop], label = "Threshold", vjust = 1, hjust=0), size = 4, color="red") 
    
    p2 = p1 + with(Temp.FudgeAxis, scale_y_continuous( breaks=breaks, labels=labels)) + 
      theme(axis.text.x = element_text(color="orange"),
            axis.title.x= element_text(color="orange")) +
      ylab("Number of genes") 
    
    p3 = ggdraw(switch_axis_position(p1 + theme_bw() + ggtitle(title_highlow) +
                                       theme(axis.ticks.y = element_blank(), 
                                             axis.text.y = element_text(color="black", size=10),
                                             plot.margin=unit(c(0,0,0,0),"mm")), axis = 'x', keep='x'))
    
    # extract gtable
    g1 <- ggplot_gtable(ggplot_build(p3))
    g2 <- ggplot_gtable(ggplot_build(p2))
    
    # save
    png(paste0(out_file, outfile_ext), width=800, height=600, res=103)
    
    ## what I dd is to plot a new viewports.
    grid.newpage()
    grid.draw(g1)
    vp=viewport(x=0.5,y=0.446,height=0.89,width=1)
    pushViewport(vp)
    grid.draw(g2)
    upViewport()
    #
    dev.off()
  }
}

GSEA_folder = grep(".Gsea.", dir(out_dir), value = T)

for(i in 1:length(GSEA_folder)){
  gsea.name = unlist(strsplit(GSEA_folder[i], split = ".Gsea."))[1]
  barchart_gsea(paste0(out_dir, "/", GSEA_folder[i]), paste0(out_dir, "/", gsea.name), nTop = 15, cutoff=Cutoff, method=order_by)
}

# data from GSEA outcome
# construct adjacency matrix
#   by calculating the jaccord index
# reduce to top10 or less than 10 pathways
Data.network = Data[order(Data$PVAL), ]

for(k in 1:length(comp.unique)){
  Data.match.id = grep(comp.unique[k], Data.network$COMPARISON)
  if(length(Data.match.id)>10){
    Data.match.id = Data.match.id[1:10]
  }

  network.matrix = matrix(0, nrow = length(Data.match.id), ncol = length(Data.match.id))
  for(i in 1:(length(Data.match.id)-1)){
    gMember.1 = unlist(strsplit(as.character(Data.network$MEMBER_overlap[Data.match.id[i]]), split = ","))
    for(j in (i+1):length(Data.match.id)){
      gMember.2 = unlist(strsplit(as.character(Data.network$MEMBER_overlap[Data.match.id[j]]), split = ","))
      gMember.int = intersect(gMember.1, gMember.2)
      network.matrix[i,j] = length(gMember.int)/(length(gMember.1)+length(gMember.2)-length(gMember.int))
    }
  }
  network.matrix = network.matrix + t(network.matrix)
  colnames(network.matrix) = substr(as.character(Data.network$NAME[Data.match.id]), 1, 40)
  rownames(network.matrix) = substr(as.character(Data.network$NAME[Data.match.id]), 1, 40)
  g.network.matrix = graph.adjacency(network.matrix, mode="undirected", weighted=T)
  
  CastY2X = function(xmin, xmax, y){
    ymax = max(y) 
    ymin = min(y)
    (xmin+(xmax-xmin)/(ymax-ymin)*(y-ymin[1])) # x gets mapped to range of ylim1
  }
  lay=layout.kamada.kawai(g.network.matrix)
  #lay=layout.kamada.kawai(g.network.matrix)
  # first plot
  E(g.network.matrix)$color = "grey"
  E(g.network.matrix)$width = CastY2X(1, 5, network.matrix[get.edgelist(g.network.matrix)])
  #V(g.network.matrix)$label = as.character(Data$NAME[Data.match.id])
  V(g.network.matrix)$label.cex = 1
  V(g.network.matrix)$size = CastY2X(10, 30, Data.network$SIZE[Data.match.id])
  V(g.network.matrix)$color = ifelse(Data.network$DIRECTION[Data.match.id]=="UP", "red", "green")
  V(g.network.matrix)$shape = "circle"
  
  png(paste0(out_dir, "/", comp.unique[k], ".png"), heigh=800, width=800)
  plot(g.network.matrix, layout=lay, main = "Pathways connected by common genes", asp=0.9)
  legend("topright",      # location of the legend on the heatmap plot
        legend = c("Upregulated", "Downregulated"), # category  labels
        col = c("red", "green"),  # color key
        pch= 19,            # line style
        cex= 1.5
  )
  dev.off()
}
