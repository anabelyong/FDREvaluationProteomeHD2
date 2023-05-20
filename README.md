# Evaluation of Different FDR strategies with newly curated ProteomeHD2 dataset

This repository includes all the code generated in order to analyse the performance of different FDR methods in terms of target protein retrieval, isoform retrieval, microprotein retrieval, etc.
Code is as explained below.

## Code and their respective outputs
These are mostly executed on Jupyter Notebook where the figures and outputs can be visualised immediately after running it. When opening the ipynb.files the figures are presented. THey are saved in this format for data reproducibility and reliability.

## Getting Started

### Raw Files and Random Raw Files Generation(OPEN TO VISUALISE THE GENERATION OF THESE FILES!!!)
These include all the code required to subset all the target and decoy files which are produced before downstream processing in pgFDR, ppFDR and cFDR methods. Thee are preliminary sources required. These ipynb. files are advised to be opened to present the crediability of code and reliance of results.

*Random_5000_Files_Generator.ipynb: Jupyter notebook containing code for generating random subsets of files. <br>
*RANDOM1000Files3MethodsAnalysis.ipynb: Jupyter notebook containing code for analyzing proteomics data using three methods on a random subset of 1000 files.<br>
*RawFilesCheckpoint.ipynb: Jupyter notebook containing code for checking raw files. <br>
*Subset100decoycut.ipynb: Jupyter notebook containing code for generating a subset of files with a decoy cut.<br>
*Subset750decoycutcode.ipynb: Jupyter notebook containing code for generating a subset of files with a decoy PSAMids produced from Percolator output files<br>
*Subset750targetcutcode.ipynb: Jupyter notebook containing code for generating a subset of files with a target PSMids produced from Percolator output files<br>
*SubsetRawFilesforLargerFilesProteomeHD2.ipynb 

## Generation of Pep-to-Prot Mapping Text Template, otherwise, go to PGR_from_combined_targets_decoys copy.R

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


**Example of SeparateProteinIDs_Classic_no_grouping_no_remap.csv, How we analysed this to count the number of target proteins**
https://github.com/anabelyong/FDREvaluationProteomeHD2/blob/main/ExampleSeparateProteinIDs.png

<img src="https://github.com/anabelyong/FDREvaluationProteomeHD2/blob/main/ExampleSeparateProteinIDs.png" width="1000"/>

sp refers to SwissProt proteins, whereas REV_sp refers to their respective decoys. Target protein count were based on the number of sp detected. Decoy count is by number of REV__ proteins detected.

### Extension of FDR performance by extracting isoform weight and length through bioinformatics pipeline:
```
import pandas as pd
from multiprocessing import Pool
import pandas as pd
import requests
from bs4 import BeautifulSoup
from time import sleep

def get_fasta_isoform(isoform_id):
    # define the UniProt API endpoint for the protein isoform
    url = f"https://www.uniprot.org/uniprot/{isoform_id}.fasta"

    # send a GET request to the UniProt API endpoint to retrieve the protein sequence
    response = requests.get(url)
    
    lines = response.text.split('\n')
    fasta_isoform = ''.join(lines[1:])

    return fasta_isoform


def get_protein_info(uniprot_id):
    try:
        print(uniprot_id)
        # define the UniProt API endpoint for the protein isoform
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"

        # send a GET request to the UniProt API endpoint to retrieve the protein sequence
        response = requests.get(url)
        lines = response.text.split('\n')
        fasta_isoform = ''.join(lines[1:])

        # Define the URL for the ProtParam tool
        url = "https://web.expasy.org/cgi-bin/protparam/protparam"

        # Define the payload for the GET request
        payload = {"sequence": fasta_isoform}

        # Submit the GET request and get the response
        response = requests.get(url, params=payload)

        soup = BeautifulSoup(response.content, "html.parser")
        num_amino_acids_raw = soup.find(string="Number of amino acids:").next
        mol_weight_raw = soup.find(string="Molecular weight:").next

        num_amino_acids = int(num_amino_acids_raw)
        mol_weight = float(mol_weight_raw)

        return {'uniprot_id': uniprot_id, 'num_amino_acids': num_amino_acids, 'mol_weight': mol_weight}

    except Exception as e:
        print(f"Error processing {uniprot_id}: {e}")
        return None

if __name__ == '__main__':
    # Load the UniProt IDs from a CSV file
    # df = pd.read_csv('ProteinIsoformsHumanSavitskiNoRemap.csv')

    df = pd.DataFrame.from_dict({'Protein Uniprot': ['O94925-3', 'P14618-2', 'P50851-2', 'O94915-2', 'P56181-2']})

    # Use multiprocessing to speed up the retrieval of protein information
    with Pool(processes=4) as pool:
        results = pool.map(get_protein_info, df['Protein Uniprot'])

    # Combine the results into a single DataFrame
    df_results = pd.DataFrame.from_dict(results)

    # Save the results to a CSV file
    df_results.to_csv('isoform_info.csv', index=False)

```



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
Example WRONG Output from pgFDR tool, as shown in proteinGroups.log:
```
2023-03-13 11:42:52,091 - WARNING - Missing peptide: nMCPGSHVAPPAQQPSSPGSGCLAPETGPGPGRSDQAR ['REV__sp|O95866-2|G6B_HUMAN', 'REV__sp|O95866-4|G6B_HUMAN']
2023-03-13 11:42:52,092 - WARNING - Missing peptide: nMMTLEMIANVIDNSNK ['REV__sp|Q9BZK3|NACP4_HUMAN']
2023-03-13 11:42:52,092 - WARNING - Missing peptide: nMVDSEISPMK ['REV__sp|Q12879|NMDE1_HUMAN']
2023-03-13 11:42:52,092 - WARNING - Missing peptide: GPTGPFGRDGQPGPFGPSGPLGEAGPPGPLPVIKGPEGK ['REV__sp|P02462|CO4A1_HUMAN']
2023-03-13 11:42:52,093 - WARNING - Missing peptide: HVTRRGAGSMSAAEIR ['REV__sp|Q07954|LRP1_HUMAN']
2023-03-13 11:42:52,093 - WARNING - Missing peptide: FSSRFGGKQDQR ['REV__sp|A8MVX0-2|ARG33_HUMAN', 'REV__sp|A8MVX0|ARG33_HUMAN']
2023-03-13 11:42:52,093 - WARNING - Missing peptide: nARAGLESIR ['REV__sp|Q14849-2|STAR3_HUMAN', 'REV__sp|Q14849-3|STAR3_HUMAN', 'REV__sp|Q14849|STAR3_HUMAN']
2023-03-13 11:42:52,093 - WARNING - Missing peptide: HGPHHHPPGMEAPFDGQPYDFRHPYPPREFPHR ['REV__sp|Q8IWX8|CHERP_HUMAN']
2023-03-13 11:42:52,093 - WARNING - Missing peptide: YGPLGRAGQIGPSGPVGPTGAIGQYGPKGDQGPQGPYGPNGPLGPDGK ['REV__sp|Q96P44-3|COLA1_HUMAN']
2023-03-13 11:42:52,093 - WARNING - Missing peptide: PSDDTGRDTKR ['REV__sp|Q5H9L2|TCAL5_HUMAN', 'REV__sp|Q6IPX3-2|TCAL6_HUMAN', 'REV__sp|Q6IPX3|TCAL6_HUMAN', 'REV__sp|Q969E4|TCAL3_HUMAN']
2023-03-13 11:42:52,094 - WARNING - Missing peptide: PDCDAGHLHLSMRR ['REV__sp|O60522-2|TDRD6_HUMAN', 'REV__sp|O60522|TDRD6_HUMAN']
2023-03-13 11:42:52,094 - WARNING - Missing peptide: nMCDTGPTPSYQNEMARSEAREPVWSYSDAANQPK ['REV__sp|Q9NPY3|C1QR1_HUMAN']
2023-03-13 11:42:52,094 - WARNING - Missing peptide: INQKDMSVPAQGPFGSGEWGCFQHCKQQVEGLK ['REV__sp|O60942-4|MCE1_HUMAN']
2023-03-13 11:42:52,094 - WARNING - Missing peptide: ALPQLPAR ['REV__sp|Q9NRJ4|TULP4_HUMAN']
2023-03-13 11:42:52,095 - WARNING - Missing peptide: LRPDEEDPPPSR ['REV__sp|A6NKC0|F90A7_HUMAN', 'REV__sp|A8MXZ1|F90AN_HUMAN']
2023-03-13 11:42:52,095 - WARNING - Missing peptide: nMMFFQVILNAVGLIIIFVHFILK ['REV__sp|Q8NBW4-3|S38A9_HUMAN', 'REV__sp|Q8NBW4-4|S38A9_HUMAN', 'REV__sp|Q8NBW4|S38A9_HUMAN']
2023-03-13 11:42:52,095 - WARNING - Missing peptide: nMDETEEAETITGTVDDKDEDSK ['REV__sp|Q8N6K0|TEX29_HUMAN']
2023-03-13 11:42:52,096 - WARNING - Missing peptide: PCAPSHGTCMEPFDCEDK ['REV__sp|Q9H2U9|ADAM7_HUMAN']
2023-03-13 11:42:52,096 - WARNING - Missing peptide: VARFQKAR ['REV__sp|Q8WXH0-11|SYNE2_HUMAN', 'REV__sp|Q8WXH0-12|SYNE2_HUMAN', 'REV__sp|Q8WXH0-13|SYNE2_HUMAN', 'REV__sp|Q8WXH0-3|SYNE2_HUMAN', 'REV__sp|Q8WXH0-4|SYNE2_HUMAN', 'REV__sp|Q8WXH0-5|SYNE2_HUMAN', 'REV__sp|Q8WXH0-7|SYNE2_HUMAN', 'REV__sp|Q8WXH0|SYNE2_HUMAN']
2023-03-13 11:42:52,096 - WARNING - Missing peptide: nCCGGSGHGCGSGGGQQGSGRDCSNPR ['REV__sp|Q9BYE3|LCE3D_HUMAN']
2023-03-13 11:42:52,097 - WARNING - Missing peptide: STTDELGPPGRR ['REV__sp|O94966-2|UBP19_HUMAN', 'REV__sp|O94966-3|UBP19_HUMAN', 'REV__sp|O94966-4|UBP19_HUMAN', 'REV__sp|O94966-5|UBP19_HUMAN', 'REV__sp|O94966-6|UBP19_HUMAN', 'REV__sp|O94966-7|UBP19_HUMAN', 'REV__sp|O94966|UBP19_HUMAN']
2023-03-13 11:42:52,097 - WARNING - Missing peptide: LREEEMEVMDDEENK ['REV__sp|P80303-2|NUCB2_HUMAN', 'REV__sp|P80303|NUCB2_HUMAN']
2023-03-13 11:42:52,097 - WARNING - Missing peptide: RVRLLIFALSGGALLLLVLLTGGVAGWVLPGVDR ['REV__sp|Q92692-2|NECT2_HUMAN']
2023-03-13 11:42:52,887 - INFO - Assigning peptides to protein groups
2023-03-13 11:42:54,528 - INFO - Shared peptides: 262734; Unique peptides: 249117
2023-03-13 11:42:56,319 - INFO - Calculating protein group-level FDRs
2023-03-13 11:42:56,559 - INFO - Decoys: 1, Entrapments: 0, Pool: 50890
2023-03-13 11:42:56,604 - INFO - #Targets at 1% decoy FDR: 50890
2023-03-13 11:42:59,049 - INFO - Protein group results have been written to: proteinGroups.txt
2023-03-13 11:42:59,050 - INFO - PickedGroupFDR execution took 114.6 seconds wall clock time
```


Example Correct Output from the pgFDR tool in proteinGroups.log

```
2023-03-17 14:56:26,741 - INFO - Could not find pyproject.toml at /Users/anabelyong/miniconda3/lib/python3.10/site-packages/pyproject.toml containing the version number.
2023-03-17 14:56:26,741 - INFO - PickedGroupFDR version None
Copyright (c) 2020-2022 Matthew The. All rights reserved.
Written by Matthew The (matthew.the@tum.de) at the
Chair of Proteomics and Bioanalytics at the Technical University of Munich.
2023-03-17 14:56:26,741 - INFO - Issued command: picked_group_fdr.py --perc_evidence combined_targets+decoys.tsv --protein_groups_out proteinGroups.txt --method savitski_no_remap --peptide_protein_map pep_to_prot_mapping.txt
2023-03-17 14:56:26,745 - INFO - Protein group level estimation method: Savitski (best Percolator PEP, no protein grouping, discard shared peptides, picked target-decoy strategy)
2023-03-17 14:56:26,745 - INFO - Parsing Percolator output file
2023-03-17 14:56:26,746 - INFO -     Reading line 0
2023-03-17 14:56:31,941 - INFO -     Reading line 500000
2023-03-17 14:56:37,196 - INFO -     Reading line 1000000
2023-03-17 14:56:43,033 - INFO -     Reading line 1500000
2023-03-17 14:56:44,923 - INFO - Assigning peptides to protein groups
2023-03-17 14:56:47,068 - INFO - Shared peptides: 298092; Unique peptides: 296767
2023-03-17 14:56:49,487 - INFO - Calculating protein group-level FDRs
2023-03-17 14:56:49,637 - INFO - Decoys: 18176, Entrapments: 0, Pool: 20089
2023-03-17 14:56:49,681 - INFO - #Targets at 1% decoy FDR: 3766
2023-03-17 14:56:50,731 - INFO - Protein group results have been written to: proteinGroups.txt
2023-03-17 14:56:50,731 - INFO - PickedGroupFDR execution took 24.0 seconds wall clock time

```
ProteinGroups.output txt File for FDR Evaluation:
```
Protein IDs	Majority protein IDs	Peptide counts (unique)	Protein names	Gene names	Fasta headers	Best peptide	Number of proteins	Q-value	Score	Reverse	Potential contaminant
sp|P35527|K1C9_HUMAN	sp|P35527|K1C9_HUMAN	71				GGSGGSHGGGSGFGGESGGSYGGGEEASGSGGGYGGGSGK	1	0.0005350454788657035	18.66152557717679		
sp|P35908|K22E_HUMAN	sp|P35908|K22E_HUMAN	72				GGGFGGGSSFGGGSGFSGGGFGGGGFGGGR	1	0.0005350454788657035	17.13278951136139		
sp|P67809|YBOX1_HUMAN	sp|P67809|YBOX1_HUMAN	21				RPQYSNPPVQGEVMEGADNQGAGEQGR	1	0.0005350454788657035	16.793147149993302		
sp|P09972|ALDOC_HUMAN	sp|P09972|ALDOC_HUMAN	40				AEVNGLAAQGKYEGSGEDGGAAAQSLYIANHAY	1	0.0005350454788657035	16.784375653274928		
sp|P14550|AK1A1_HUMAN	sp|P14550|AK1A1_HUMAN	26				HIDCAAIYGNEPEIGEALKEDVGPGK	1	0.0005350454788657035	15.73462786131106		
sp|P04264|K2C1_HUMAN	sp|P04264|K2C1_HUMAN	84				FSSCGGGGGSFGAGGGFGSR	1	0.0005350454788657035	14.904022113958698		
sp|P13645|K1C10_HUMAN	sp|P13645|K1C10_HUMAN	53				GSSGGGCFGGSSGGYGGLGGFGGGSFR	1	0.0005350454788657035	14.389258946452074		
sp|P07737|PROF1_HUMAN	sp|P07737|PROF1_HUMAN	28				CSVIRDSLLQDGEFSMDLR	1	0.0005350454788657035	14.083014694909		
sp|P50395|GDIB_HUMAN	sp|P50395|GDIB_HUMAN	7				VPSTEAEALASSLMGLFEK	1	0.0005350454788657035	13.773804937633944		
sp|P14866|HNRPL_HUMAN	sp|P14866|HNRPL_HUMAN	17				TDNAGDQHGGGGGGGGGAGAAGGGGGGENYDDPHK	1	0.0005350454788657035	13.68738898894721		
sp|P36578|RL4_HUMAN	sp|P36578|RL4_HUMAN	63				NNRQPYAVSELAGHQTSAESWGTGR	1	0.0005350454788657035	13.428531086092361		
sp|P09651|ROA1_HUMAN	sp|P09651|ROA1_HUMAN	1				GGGGYGGSGDGYNGFGNDGGYGGGGPGYSGGSR	1	0.0005350454788657035	13.331491081692215		
sp|P08567|PLEK_HUMAN	sp|P08567|PLEK_HUMAN	27				KSEEENLFEIITADEVHYFLQAATPK	1	0.0005350454788657035	13.306237007888331		
sp|Q12905|ILF2_HUMAN	sp|Q12905|ILF2_HUMAN	25				NQDLAPNSAEQASILSLVTK	1	0.0005350454788657035	13.024205741353667		
sp|P60228|EIF3E_HUMAN	sp|P60228|EIF3E_HUMAN	33				LDLLSDTNMVDFAMDVYK	1	0.0005350454788657035	13.019927965292863		
sp|Q5D862|FILA2_HUMAN	sp|Q5D862|FILA2_HUMAN	56				HQEEESETEEDEEDTPGHK	1	0.0005350454788657035	13.014599126242835		
sp|O43815|STRN_HUMAN	sp|O43815|STRN_HUMAN	3				GLGPLAEAAAAGDGAAAAGAAR	1	0.0005350454788657035	12.937718930027357		
sp|P26196|DDX6_HUMAN	sp|P26196|DDX6_HUMAN	39				GPVKPTGGPGGGGTQTQQQMNQLK	1	0.0005350454788657035	12.90085478932859		
sp|Q15437|SC23B_HUMAN	sp|Q15437|SC23B_HUMAN	36				GPCVSENELGVGGTSQWK	1	0.0005350454788657035	12.866566484695655		
sp|Q06203|PUR1_HUMAN	sp|Q06203|PUR1_HUMAN	30				GQESAGIVTSDGSSVPTFK	1	0.0005350454788657035	12.755468068883696		
sp|Q9NYB0|TE2IP_HUMAN	sp|Q9NYB0|TE2IP_HUMAN	19				LGPASAADTGSEAKPGALAEGAAEPEPQR	1	0.0005350454788657035	12.7195513910777		
sp|Q9H7Z7|PGES2_HUMAN	sp|Q9H7Z7|PGES2_HUMAN	24				KVPILVAQEGESSQQLNDSSVIISALK	1	0.0005350454788657035	12.661898258666172		
sp|O00151|PDLI1_HUMAN	sp|O00151|PDLI1_HUMAN	31				VITNQYNNPAGLYSSENISNFNNALESK	1	0.0005350454788657035	12.642909800392198		

```
NEW METHODS GENERATED FOR THE PURPOSE OF THIS DISSERTATION: 
```
label = "Picked Protein Group FDR"

scoreType = "Perc bestPEP"
grouping = "no"
sharedPeptides = "discard"
pickedStrategy = "picked_group"


```
Format of example  .toml file I have generated for this dissertation. For reference, picked_no_grouping_no_remap represents Picked Group TDS + no grouping strategy <-- this was what I manipulated for evaluating the parameters./

## Authors

Contributors names and contact info:

ex.Anabel Yong (s1911568@ed.ac.uk)

## License

This project is licensed under University of Edinburgh.

## Acknowledgments
* Georg Kustatcher (georg.kustatscher@ed.ac.uk)
* Savvas Kourtis


