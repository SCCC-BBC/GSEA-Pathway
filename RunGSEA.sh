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


done < $1
