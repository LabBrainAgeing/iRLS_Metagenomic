# iRLS_Metagenomic
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
