#!/bin/bash

# Set compare group name
compare_groups="female_AD_vs_CT"

# Path definitions
coloc_dir="/work/aliu10/AD_sc_Project/results/Post_GWAS/${compare_groups}/coloc"
smr_dir="/work/aliu10/AD_sc_Project/SMR/smr_files"
eqtl_file="/work/aliu10/AD_sc_Project/GTEx_Analysis_v10_eQTL_updated/Brain_Cortex.v10.eGenes.txt.gz"
bfile_prefix="/work/aliu10/AD_sc_Project/susieR/1000Genome/Plink_files/allchr_EUR_renamed"
gwas_summary="${smr_dir}/mygwas.ma"
eqtl_prefix="${smr_dir}/myeqtl"
probe_list="${smr_dir}/myprobe_${compare_groups}.list"
snp_list="${smr_dir}/mysnp_${compare_groups}.list"
smr_out="${smr_dir}/mysmr_${compare_groups}"

mkdir -p "$smr_dir"

# Extract gene list
ls "${coloc_dir}"/*.csv | sed 's#.*/##' | cut -d'_' -f1 | sort -u > "$probe_list"

# Extract SNP list
tmp_snp_all="all_coloc_snps.tsv"
> "$tmp_snp_all"

for file in "${coloc_dir}"/*.csv; do
  awk -F',' 'NR > 1 {gsub(/"/, "", $2); gsub(/"/, "", $3); print $2 "\n" $3}' "$file" >> "$tmp_snp_all"
done

sort "$tmp_snp_all" | uniq > all_coloc_snps_unique.txt

# Build SNP chr:pos mapping
zcat "$eqtl_file" | \
awk 'NR > 1 && $20 != "NA" && $15 != "NA" && $16 != "NA" {
  gsub(/^chr/, "", $15);
  print $20 "\t" $15 ":" $16
}' | sort -u > snp_map.tsv

# Join to get chr:pos formatted SNPs
join -t $'\t' -1 1 -2 1 \
  <(LC_ALL=C sort all_coloc_snps_unique.txt) \
  <(LC_ALL=C sort snp_map.tsv) \
  > ${smr_dir}/coloc_snps_chrpos_${compare_groups}.txt

cut -f2 ${smr_dir}/coloc_snps_chrpos_${compare_groups}.txt > "$snp_list"

# Clean up temp
rm "$tmp_snp_all" all_coloc_snps_unique.txt snp_map.tsv

# Run SMR
./smr --bfile "$bfile_prefix" --gwas-summary "$gwas_summary" --beqtl-summary "$eqtl_prefix" \
  --extract-probe "$probe_list" \
  --maf 0.05 \
  --diff-freq 0.5 \
  --diff-freq-prop 1 \
  --peqtl-smr 0.05 \
  --peqtl-heidi 0.05 \
  --out "$smr_out" \
  --thread-num 64
  
# Change rsID
 
# Define input/output
smr_result="${smr_out}.smr"
snp_map_file="${smr_dir}/coloc_snps_chrpos_${compare_groups}.txt"
out_file="${smr_out}_rsID.smr"

# Step 1: 
snp_map_rev="${snp_map_file%.txt}_rev.txt"
awk -F'\t' '{print $2 "\t" $1}' "$snp_map_file" > "$snp_map_rev"

# Step 2: 
awk -F'\t' -v OFS='\t' -v map="$snp_map_rev" '
  BEGIN {
    while ((getline < map) > 0) {
      chrpos2rs[$1] = $2
    }
  }
  NR == 1 {
    print; next
  }
  {
    if ($5 in chrpos2rs) $5 = chrpos2rs[$5]
    print
  }
' "$smr_result" > "$out_file"

echo "Replaced topSNP column in $smr_result and saved to $out_file"

rm "$snp_map_rev" 
mv "$out_file" "$smr_result"
