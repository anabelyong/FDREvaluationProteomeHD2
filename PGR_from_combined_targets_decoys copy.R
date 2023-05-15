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
sec_command <- glue::glue("sed -i''-e 's/,/\\t/g;s/rev_/REV__/g;s/_;/_HUMAN\\t/g' {location_target_decoys}combined_targets+decoys.tsv")
sec_command = paste0("bash -c ",'"',sec_command,'"')
system(sec_command)



FDR_strategies <- c("picked_protein_group_no_remap", "classic_no_grouping_no_remap", "savitski_no_remap")
for(fdr in FDR_strategies[1]){
  Protein_group_command = 
    glue::glue("cd /Users/anabelyong/miniconda3/lib/python3.10/site-packages/picked_group_fdr/python -um picked_group_fdr --perc_evidence {location_target_decoys}combined_targets+decoys.tsv --protein_groups_out {location_target_decoys}proteinGroups_{fdr}.txt --method {fdr} --fasta {location_target_decoys}/{fasta} | tee proteinGroups.log")
  system(Protein_group_command)
}







FDR_strategies <- c("picked_protein_group_no_remap", "classic_no_grouping_no_remap", "savitski_no_remap")
for(fdr in FDR_strategies[1]){
  Protein_group_command = 
    glue::glue("cd ~/miniconda3/lib/python3.10/site-packages/picked_group_fdr/ ;/Users/anabelyong/miniconda3/bin/python -um picked_group_fdr --perc_evidence {location_target_decoys}combined_targets+decoys.tsv --protein_groups_out {location_target_decoys}proteinGroups_{fdr}.txt --method {fdr} --fasta {location_target_decoys}/{fasta} | tee proteinGroups.log")
  system(Protein_group_command)
}






