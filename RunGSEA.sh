while read line; do

f=`echo "$line"`

gct_name=$(echo "$f" | cut -f1)

cls_name=$(echo "$f" | cut -f2)


gct_dir_name=$(dirname "$gct_name")
gct_file_name=$(basename "$gct_name")

cls_dir_name=$(dirname "$cls_name")
cls_file_name=$(basename "$cls_name")


echo "$gct_dir_name" 
echo "$cls_dir_name"

echo "$gct_file_name"
echo "$cls_file_name"

Cutoff_name=`echo "$cls_file_name" | awk -F"." '{print $1}' | awk -F"_" '{print $1}'`

echo "$Cutoff_name"

subtype1_name=`echo "$cls_file_name" | awk -F"." '{print $1}' | awk -F"_" '{print $2}' | awk -F"vs" '{print $1}'`
subtype2_name=`echo "$cls_file_name" | awk -F"." '{print $1}' | awk -F"_" '{print $2}' | awk -F"vs" '{print $2}'`


echo "$subtype1_name"
echo "$subtype2_name"

java -Xmx1024m -cp ~/Downloads/gsea2-2.2.2.jar xtools.gsea.Gsea -res ~/GSEA_Input/"$gct_file_name" \
-cls ~/GSEA_Input/"$cls_file_name" \
-gmx ~/Downloads/c2.cp.v5.1.symbols.gmt.txt \
-collapse true -mode Max_probe -norm meandiv -nperm 5000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label "$Cutoff_name"_"$subtype1_name"_vs_"$subtype2_name"_gene_set_permutation_5000_analysis \
-metric Signal2Noise -sort real -order descending -chip ~/Downloads/HG_U133A_2.chip \
-include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false \
-set_max 5000 -set_min 15 -zip_report false -out /Users/axy148/gsea_home/output/mar22 -gui false

#java -Xmx512m xtools.gsea.Gsea -res /Volumes/Bioinformatics$/2015/Sophia/Results/5_St34.gct -cls /Volumes/Bioinformatics$/2015/Sophia/Results/st3vsst4_phenotype.cls#st3_versus_st4 -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.v5.1.symbols.gmt -collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label st3_vs_st4_gene_set_permutation_5000_analysis -metric Signal2Noise -sort real -order descending -chip gseaftp.broadinstitute.org://pub/gsea/annotations/HG_U133A_2.chip -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 5000 -set_min 15 -zip_report false -out /Users/axy148/gsea_home/output/mar22 -gui false

#cat<<EO
#EOF> ~/Script_bash/Run_"$fbname"_process_vcf.sh
#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J ANNOVAR
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

#wme $f)/$(basename "$f" | cut -d. -f1)_str.avinputc -l /scratch/projects/bbc/Project/zLi/Mutation/S2vsS6/results/mutect_output_passed.vcf.recode.vcf
#wc -l "$f" 

#tail -n +52 "$f" | cut -f1,2,4,5,10-11 > /scratch/projects/bbc/Project/zLi/Mutation/"$fbname"/results/Str_processed.vcf

#wc -l /scratch/projects/bbc/Project/zLi/Mutation/"$fbname"/results/Str_processed.vcf

#R --slave --args /scratch/projects/bbc/Project/zLi/Mutation/"$fbname"/results/Str_processed.vcf "$fbname" /scratch/projects/bbc/Project/zLi/Mutation/"$fbname"/results/Str_processed_vcf_freq_3.txt < ~/Script_bash/Run_process_one_vcf_files.R

#java -Xmx512m xtools.gsea.Gsea -res /Volumes/Bioinformatics$/2015/Sophia/Results/5_St34.gct -cls /Volumes/Bioinformatics$/2015/Sophia/Results/st3vsst4_phenotype.cls#st3_versus_st4 -gmx gseaftp.broadinsti\
#tute.org://pub/gsea/gene_sets/c2.cp.v5.1.symbols.gmt -collapse true -mode Max_probe -norm meandiv -nperm 5000 -permute phenotype -rnd_type no_balance -scoring_scheme weighted -rpt_label st3_vs_st4_phenot\
#ype_permutation_5000_analysis -metric Signal2Noise -sort real -order descending -chip gseaftp.broadinstitute.org://pub/gsea/annotations/HG_U133A_2.chip -include_only_symbols true -make_sets true -median \
#false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 5000 -set_min 15 -zip_report false -out /Users/axy148/gsea_home/output/mar22 -gui false

#cut -f1,2,10-11 "$f" | head -100 | cut -f2 | cut -d"," -f2,3,4,5 | cut -d"," -f1 | cut -d":" -f1

#mutect_output_passed.vcf.recode.vcf

#cut -f1-2,4- "$f" > /scratch/projects/bbc/Project/zLi/Mutation/S7vsS9/results/mutect_output_2.txt

#cut -f1-2,4- "$f" > $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3.txt

#perl /scratch/projects/bbc/NGS-tools/Convert_mut_2_vcf.pl $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3.txt $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3.vcf

#$annovar/convert2annovar.pl -format vcf4old $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3.vcf > $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3_mutect.avinput

#$annovar/annotate_variation.pl -out $(dirname $f)/myanno_mutect -build mm10 $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3_mutect.avinput $annovar/mm10db/

#-out $(dirname $f)/myanno -remove -protocol refGene -operation g -nastring . -vcfinput

#EOF
#bsub -P ANNOVAR < ~/Script_bash/Run_"$fbname"_process_vcf.sh
done < $1
