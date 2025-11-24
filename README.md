# iRLS_Metagenomic
## :page_facing_up: R Analysis Pipeline 
```
├── 0_Preprocessing_iRLS2_Metagenomics_OK.R
├── 1_iRLS2_Phyloseq_PhylogeneticTree_OK.R
├── 2_iRLS2_Metagenomics_Shrink_BigDataAnalysis_OK
│ ├── 3_iRLS_Metagenomica_TreatmentComparison_OK.R
│ ├── 4_20240724_iRLS_BDI_Analysis_OK.R
│ ├── 5_iRLS2_INSB_Metagenomics_Analysis_Family_OK.R
│ ├── 6_iRLS_Metagenomica_PHYLUM_OK.R
│ ├── 7_iRLS_Metagenomics_ALDEx2.R
│ ├── 8_iRLS2_picrust2_OK.R
├── 9_iRLS2_A_Diversity_OK.R
└── 10_iRLS2_B_Diversity_OK.R
```
# Gut Microbiota Analysis in Restless Legs Syndrome (RLS)

This repository contains code and documentation for the analysis of gut microbiota composition in patients with Restless Legs Syndrome (RLS), insomnia (INS), and healthy controls (CTRL). The project is part of a research study exploring microbiota alterations across different RLS phenotypes using 16S rRNA sequencing and downstream statistical tools.

---

## 📌 Study Objectives

To analyze the fecal microbiota composition in patients with **idiopathic Restless Legs Syndrome (RLS)** and investigate its relationship with different **RLS clinical phenotypes**.

---

## 🧪 Methods

- **Participants**:  
  - 37 RLS patients (28 females, mean age: 64.78 years)  
  - 31 INS patients (22 females, mean age: 60.64 years)  
  - 33 CTRL individuals (24 females, mean age: 62.54 years)  

- **Microbiota Analysis**:
  - **Sample collection**: Stool samples collected and processed.
  - **Sequencing**: 16S rRNA gene sequencing using **Illumina MiSeq** platform.
  - **Pipeline**:  
    - **DADA2** for sequence preprocessing and ASV table generation  
    - **DESeq2** and **ALDEx2** for differential abundance analysis  
    - Correction for **age**, **sex**, **BMI**, **sequencing run**, and **presence of mood disorders**

---
