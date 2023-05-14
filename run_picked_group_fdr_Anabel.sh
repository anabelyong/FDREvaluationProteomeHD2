cat combined_targets.tsv <(tail -n+2 combined_decoys.tsv) > combined_targets+decoys.tsv
sed -i 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys.tsv

python -um picked_group_fdr --perc_evidence combined_targets+decoys.tsv --protein_groups_out proteinGroups.txt --method picked_protein_group_no_remap | tee proteinGroups.log

#modified shell script picked group FDR package. Modify them according to range of subset. 
#Classic_no_grouping_no_remap, 100, 
cat subset_100_target_rawfiles.tsv <(tail -n+2 subset_100_decoy_rawfiles.tsv) > combined_targets+decoys_100.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_100.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_100.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_100.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_100.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_100.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_200_target_rawfiles.tsv <(tail -n+2 subset_200_decoy_rawfiles.tsv) > combined_targets+decoys_200.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_200.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_200.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_200.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_200.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_200.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_300_target_rawfiles.tsv <(tail -n+2 subset_300_decoy_rawfiles.tsv) > combined_targets+decoys_300.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_300.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_300.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_300.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_300.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_300.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_400_target_rawfiles.tsv <(tail -n+2 subset_400_decoy_rawfiles.tsv) > combined_targets+decoys_400.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_400.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_400.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_400.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_400.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_400.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log


cat subset_500_target_rawfiles.tsv <(tail -n+2 subset_500_decoy_rawfiles.tsv) > combined_targets+decoys_500.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_500.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_500.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_500.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_500.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_500.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_600_target_rawfiles.tsv <(tail -n+2 subset_600_decoy_rawfiles.tsv) > combined_targets+decoys_600.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_600.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_600.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_600.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_600.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_600.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_700_target_rawfiles.tsv <(tail -n+2 subset_700_decoy_rawfiles.tsv) > combined_targets+decoys_700.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_700.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_700.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_700.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_700.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_700.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_800_target_rawfiles.tsv <(tail -n+2 subset_800_decoy_rawfiles.tsv) > combined_targets+decoys_800.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_800.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_800.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_800.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_800.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_800.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log


cat subset_900_target_rawfiles.tsv <(tail -n+2 subset_900_decoy_rawfiles.tsv) > combined_targets+decoys_900.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_900.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_900.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_900.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_900.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_900.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log


cat subset_1000_target_rawfiles.tsv <(tail -n+2 subset_1000_decoy_rawfiles.tsv) > combined_targets+decoys_1000.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_1000.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1000.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1000.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1000.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1000.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log


cat subset_1100_target_rawfiles.tsv <(tail -n+2 subset_1100_decoy_rawfiles.tsv) > combined_targets+decoys_1100.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_1100.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1100.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1100.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1100.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1100.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log


cat subset_1200_target_rawfiles.tsv <(tail -n+2 subset_1200_decoy_rawfiles.tsv) > combined_targets+decoys_1200.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_1200.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1200.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1200.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1200.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1200.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log


cat subset_1300_target_rawfiles.tsv <(tail -n+2 subset_1300_decoy_rawfiles.tsv) > combined_targets+decoys_1300.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_1300.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1300.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1300.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1300.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1300.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log


cat subset_1400_target_rawfiles.tsv <(tail -n+2 subset_1400_decoy_rawfiles.tsv) > combined_targets+decoys_1400.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_1400.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1400.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1400.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1400.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1400.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log


cat subset_1500_target_rawfiles.tsv <(tail -n+2 subset_1500_decoy_rawfiles.tsv) > combined_targets+decoys_1500.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_1500.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1500.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1500.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1500.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1500.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_1600_target_rawfiles.tsv <(tail -n+2 subset_1600_decoy_rawfiles.tsv) > combined_targets+decoys_1600.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_1600.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1600.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1600.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1600.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1600.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_1700_target_rawfiles.tsv <(tail -n+2 subset_1700_decoy_rawfiles.tsv) > combined_targets+decoys_1700.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_1700.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1700.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1700.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1700.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1700.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_1800_target_rawfiles.tsv <(tail -n+2 subset_1800_decoy_rawfiles.tsv) > combined_targets+decoys_1800.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_1800.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1800.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1800.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1800.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1800.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_1900_target_rawfiles.tsv <(tail -n+2 subset_1900_decoy_rawfiles.tsv) > combined_targets+decoys_1900.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_1900.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1900.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1900.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_1900.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_1900.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_2000_target_rawfiles.tsv <(tail -n+2 subset_2000_decoy_rawfiles.tsv) > combined_targets+decoys_2000.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_2000.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_2000.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_2000.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_2000.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_2000.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_3000_target_rawfiles.tsv <(tail -n+2 subset_3000_decoy_rawfiles.tsv) > combined_targets+decoys_3000.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_3000.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_3000.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_3000.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_3000.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_3000.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_4000_target_rawfiles.tsv <(tail -n+2 subset_4000_decoy_rawfiles.tsv) > combined_targets+decoys_4000.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_4000.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_4000.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_4000.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_4000.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_4000.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_5000_target_rawfiles.tsv <(tail -n+2 subset_5000_decoy_rawfiles.tsv) > combined_targets+decoys_5000.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_5000.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_5000.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_5000.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_5000.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_5000.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log

cat subset_6000_target_rawfiles.tsv <(tail -n+2 subset_6000_decoy_rawfiles.tsv) > combined_targets+decoys_6000.tsv
sed -i ''-e 's/,/\t/g;s/rev_/REV__/g;s/_;/_HUMAN\t/g' combined_targets+decoys_6000.tsv
python -um picked_group_fdr --perc_evidence combined_targets+decoys_6000.tsv --protein_groups_out proteinGroups.txt --method classic_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_6000.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
python -um picked_group_fdr --perc_evidence combined_targets+decoys_6000.tsv --protein_groups_out proteinGroups.txt --method classic_rescued_subset_no_remap   --peptide_protein_map pep_to_prot_mapping_6000.txt --special-aas '' --enzyme trypsinp | tee proteinGroups.log
