# Evaluation of Different FDR strategies with newly curated ProteomeHD2 dataset

THis repository includes all the code generated in order to analyse the performance of different FDR methods in terms of target protein retrieval, isoform retrieval, microprotein retrieval, etc.
Code is as explained below.

## Code and their respective outputs
These are mostly executed on Juoyter Notebook where the figures and outputs can be visualised immediately after running it.

## Getting Started

### Raw Files and Random Raw Files Generation OPEN TO VISUALISE THE GENERATION OF THESE FILES!
These include all the code required to subset all the target and decoy files which are produced before downstream processing in pgFDR, ppFDR and cFDR methods. Thee are preliminary sources required. 

*Random_5000_Files_Generator.ipynb: Jupyter notebook containing code for generating random subsets of files.<br>
*RANDOM1000Files3MethodsAnalysis.ipynb: Jupyter notebook containing code for analyzing proteomics data using three methods on a random subset of 1000 files.<br>
*RawFilesCheckpoint.ipynb: Jupyter notebook containing code for checking raw files.<br>
*Subset100decoycut.ipynb: Jupyter notebook containing code for generating a subset of files with a decoy cut.<br>
*Subset750decoycutcode.ipynb: Jupyter notebook containing code for generating a subset of files with a decoy PSAMids produced from Percolator output files<br>
*Subset750targetcutcode.ipynb: Jupyter notebook containing code for generating a subset of files with a target PSMids produced from Percolator output files<br>
*SubsetRawFilesforLargerFilesProteomeHD2.ipynb 

## Generation of Pep-to-Prot Mapping Text Template 

```
# Protein_inference from percolator files
library(data.table);library(stringr);library(tidyverse)
location_target_decoys ="~/Desktop/honours/microprot_group_picked_fdr/" #path to percolator targets+ decoys, include fasta folder in same folder

# fasta = "2023-02-21-decoys-contam-sp_isoforms_with_microproteins_.fasta.fas"
# reviewed_9606 <- seqinr::read.fasta(glue::glue("{location_target_decoys}{fasta}"),
#                                     seqtype = "AA",as.string = TRUE)
# complete_fasta = data.table::data.table(Annotation  = unlist(seqinr::getAnnot(reviewed_9606)),
#                             Seq = seqinr::getSequence(reviewed_9606, as.string = T) |>  unlist())
# complete_fasta[,Annotation:=fifelse(str_detect(Annotation,"^>(pI|sR|sU|nC)"),paste0(Annotation,"_HUMAN microprot"),Annotation)]
# seqinr::write.fasta(sequences =  as.list(complete_fasta$Seq),
#                     names  = as.list(complete_fasta$Annotation %>% str_remove_all("^>")),
#                     file.out = glue::glue("{location_target_decoys}micro_prot_decoy.fasta"))
fasta = "micro_prot_decoy.fasta"
targets = "combined_targets.tsv" 
decoys  ="combined_decoys.tsv"  #path to percolator decoys

concat_command = glue::glue("bash -c 'cat {location_target_decoys}{targets} <(tail -n+2 {location_target_decoys}{decoys}) > {location_target_decoys}combined_targets+decoys.tsv'")
system(concat_command)
# combined_target_decoy = data.table::fread(glue::glue("{location_target_decoys}/combined_targets+decoys.tsv"),sep = ",")
# combined_target_decoy_micro[,proteinIds:= stringr::str_replace_all(proteinIds,"\t","_;" )]
# data.table::fwrite(combined_target_decoy_micro, glue::glue("{location_target_decoys}/combined_targets+decoys.tsv"), sep = ",")
# combined_target_decoy_micro = data.table::fread(glue::glue("{location_target_decoys}/combined_targets+decoys.tsv"))
# combined_target_decoy_micro = combined_target_decoy_micro[stringr::str_detect(proteinIds,"^(sI|sR|sU|nC)")]
combined_target_decoy_micro = data.table::fread(glue::glue("{location_target_decoys}/combined_targets+decoys.tsv"), select = c("peptide","proteinIds"))
combined_target_decoy_micro[,`:=`(peptide_strip = str_remove_all(peptide,"\\[[:print:]*?\\]") |>
                                    str_remove("^.\\.") |>
                                    str_remove("\\..$"),
                                  prot_id = str_replace_all(proteinIds,"rev_(sp|pI|sR|nC)\\|","REV__") |>
                                    str_remove_all("(sp|pI|sR|nC)\\|") |>
                                    str_remove_all("\\|[:print:]*?(?=(;|$))"))]

pep_to_prot = combined_target_decoy_micro[,.(peptide_strip,prot_id)] |> unique()
fwrite(pep_to_prot, file = glue::glue("~/Desktop/honours/microprot_group_picked_fdr/pep_to_prot_mapping.txt"),
       col.names = F,sep = "\t")


```

### FDR evaluation
*3_Classic_TDS_Methods.ipynb: Jupyter notebook containing code for performing classic TDS-grouping strategy analyses -3 different methods <br>
*9_Methods_Results_1.ipynb: Jupyter notebook containing code for analyzing and comparing the results of 9 FDR STRATEGIES <br>
*9_Methods_Results_2.ipynb: Jupyter notebook containing additional code for analyzing and comparing Results of 9 FDR STRATEGIES  PART 2 <br>
*9_Methods_Target_Count_Bar_Plots.ipynb: Jupyter notebook containing code for generating bar plots of target protein counts AT 1% FDR AND 100 FDRN<br>
*InvestigatingMicroproteinCount3strategies.ipynb: Jupyter notebook containing code for investigating microprotein counts using three different strategies.<br>
*Isoform_Analysis.ipynb: Jupyter notebook containing code for analyzing isoforms in proteomics data.<br>
*Isoform_Count_6_Methods_22796_Files.ipynb: Jupyter notebook containing code for counting isoforms using six different methods in a dataset of 22796 files.<br>
*Isoform_Count_22796_Files.ipynb: Jupyter notebook containing code for counting isoforms in a dataset of 22796 files.<br>
*Isoform_Count_Count.ipynb: Jupyter notebook containing code for counting isoforms in a dataset.<br>
*Microprotein_Count_Dataframes.ipynb: Jupyter notebook containing code for counting microproteins using dataframes.<br>
*microprotein_large_datasets_count.ipynb: Jupyter notebook containing code for counting microproteins in large datasets.<br>
*protein_groups_analysis.ipynb: Jupyter notebook containing code for analyzing protein groups.<br>
*Protein_inference_analysis.ipynb: Jupyter notebook containing code for analyzing protein inference in proteomics data<br>
An in-depth paragraph about your project and overview of use.

### Example dataset analysis before running ProteomeHD2, this only has up to 2000 raw files. TROUBLESHOOTING THE METHODS:
THese can be seen executed in:
*Classic_no_grouping_no_remap_100.ipynb <br>
*Classic_no_grouping_no_remap_200.ipynb <br>
*Classic_no_grouping_no_remap_500.ipynb <br>
*Classic_no_grouping_no_remap_750.ipynb <br>
*Classic_no_grouping_no_remap_1000.ipynb <br>
*Classic_no_grouping_no_remap_2000.ipynb <br>
*Savitski_no_remap_100.ipynb <br>
*Savitski_no_remap_200.ipynb <br>
*Savitski_no_remap_500.ipynb <br>
*Savitski_no_remap_750.ipynb <br>
*Savitski_no_remap_1000.ipynb <br>
*Savitski_no_remap_2000.ipynb <br>
*Picked_protein_group_no_remap_100.ipynb
*Picked_protein_group_no_remap_200.ipynb
*Picked_protein_group_no_remap_500.ipynb
*Picked_protein_group_no_remap_750.ipynb
*Picked_protein_group_no_remap_1000.ipynb
*Picked_protein_group_no_remap_2000.ipynb


```

## TERMINAL COMMAND LINE FOR PGFDR TOOL in PYTHON BASED ENVIRONMENT

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

```
ommand to run if program contains helper info (CAN BE SEEN IN RUN_PICKED_GROUP_FDR_ANABEL.SH
Otherwise can be run on RStudio.
```
sec_command <- glue::glue("sed -i''-e 's/,/\\t/g;s/rev_/REV__/g;s/_;/_HUMAN\\t/g' {location_target_decoys}combined_targets+decoys.tsv")
sec_command = paste0("bash -c ",'"',sec_command,'"')
system(sec_command)



FDR_strategies <- c("picked_protein_group_no_remap", "classic_grouping_no_remap", "savitski_no_remap")
for(fdr in FDR_strategies[1]){
  Protein_group_command = 
    glue::glue("cd /Users/anabelyong/miniconda3/lib/python3.10/site-packages/picked_group_fdr/python -um picked_group_fdr --perc_evidence {location_target_decoys}combined_targets+decoys.tsv --protein_groups_out {location_target_decoys}proteinGroups_{fdr}.txt --method {fdr} --fasta {location_target_decoys}/{fasta} | tee proteinGroups.log")
  system(Protein_group_command)
```




## Authors

Contributors names and contact info:

ex.Anabel Yong (s1911568@ed.ac.uk)

## License

This project is licensed under University of Edinburgh.

## Acknowledgments
* Georg Kustatcher (georg.kustatscher@ed.ac.uk)
* Savvas Kourtis


