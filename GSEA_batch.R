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
      GSEA batch submit
 
      Arguments:
      --indir=somevalue   - character, input directory containing all .cls, .gct files 
      --outdir=indir  - character, output directory
      --pathway_type=CP - character, 'BP' or 'CP'
      --prerank=FALSE - logic, indicate whether the type of analysis is preranked
      --permute_method=gene_set  - character, either 'gene_set' , 'phenotype'
      --help              - print this text
 
      Example:
      Rscript GSEA_batch.R --indir=/home/jamesban/Desktop/Projects/Yuqi --outdir=/home/jamesban/Desktop/Projects/Yuqi/GSEA --pathway_type=CP --prerank=FALSE --permute_method=gene_set

      Output:
	GSEA results
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

if(is.null(argsL$pathway_type)) {
	cat("Default canonical pathways for GSEA\n")
	gmt_file = "/home/jamesban/Desktop/Projects/GSEA/c2.cp.v5.1.symbols.gmt"
	p_type = "CP"
}else if(argsL$pathway_type=="CP"){
	gmt_file = "/home/jamesban/Desktop/Projects/GSEA/c2.cp.v5.1.symbols.gmt"
	p_type = "CP"
}else if(argsL$pathway_type=="BP"){
	gmt_file = "/home/jamesban/Desktop/Projects/GSEA/c5.bp.v5.1.symbols.gmt"
	p_type = "BP"
}else{
	stop("Do not regonize the parameter pathway_type, use --help to view usage example", call.=F)
}

if(is.null(argsL$prerank)){
	cat("Default analysis input is .cls file.\n")
	argsL$prerank = FALSE
}else{
	Prerank = (argsL$prerank=="TRUE")
}

if(is.null(argsL$permute_method)) {
	cat("Default permute method: gene_set for GSEA\n")
	p_method = "gene_set"
}else if(argsL$permute_method=="gene_set"){
	p_method = "gene_set"
}else if(argsL$permute_method=="phenotype"){
	p_method = "phenotype"
}else{
	stop("Do not regonize the parameter permute_method, use --help to view usage example", call.=F)
}

if(is.null(argsL$outdir)) {
	cat("Default output directory: indir\n")
	out_dir = argsL$indir
}else{
	out_dir = argsL$outdir
}

#dir for cls and gct files
gsea_indir = argsL$indir # "/Users/aiminyan/Desktop/Project/McCarthy/GSEA/"

if(!Prerank){
	filename_prefix = gsub(".gct", "", grep(".gct", dir(gsea_indir),value = T))

# 
	gct.files = paste0(gsea_indir,"/", filename_prefix, ".gct")
	cls.files = paste0(gsea_indir,"/", filename_prefix, ".cls")
}else{
        filename_prefix = gsub(".rnk", "", grep(".rnk", dir(gsea_indir), value = T))
	rnk.files = paste0(gsea_indir,"/", filename_prefix, ".rnk")
}
dir.create(out_dir, showWarnings = FALSE)

#code lines
if(!Prerank){
	command_line = paste0("java -Xmx4096m -cp /home/jamesban/Softwares/gsea2-2.2.2.jar xtools.gsea.Gsea -res ", gct.files, " -cls ", cls.files, "#high_versus_low -gmx ", gmt_file, " -collapse false -mode Max_probe -norm meandiv -nperm 10000 -permute ", p_method, " -rnd_type no_balance -scoring_scheme weighted -rpt_label ", filename_prefix, "_",p_type, " -metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 40 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 10 -zip_report false -out ", out_dir, " -gui false","\n")
}else{
	command_line = paste0("java -Xmx4096m -cp /home/jamesban/Softwares/gsea2-2.2.2.jar xtools.gsea.GseaPreranked -rnk ", rnk.files, " -gmx ", gmt_file, " -collapse false -mode Max_probe -norm meandiv -nperm 10000 -scoring_scheme weighted -rpt_label ", filename_prefix, "_",p_type, " -include_only_symbols true -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 10 -zip_report false -out ", out_dir, " -gui false","\n")
}

write.table("#!/bin/bash\n", paste0(gsea_indir, "/command_gsea.txt"), append=F, row.names=FALSE,col.name=FALSE, quote=FALSE)
for(i in 1:length(command_line)){
  write.table(command_line[i], paste0(gsea_indir, "/command_gsea.txt"), append=T, row.names=FALSE,col.name=FALSE, quote=FALSE)
}
