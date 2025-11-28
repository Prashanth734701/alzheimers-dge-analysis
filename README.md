
#  Differential Gene Expression Analysis in Alzheimer's Disease

##  Project Overview
This project investigates differential gene expression (DGE) patterns in **Alzheimer's disease (AD)** using RNA-seq data from **BioProject PRJNA1023207**. The analysis includes:

- AD vs healthy controls  
- Sex-specific differences (Male vs Female)

A fully automated **Nextflow** pipeline was used for preprocessing, followed by **DESeq2** for statistical modeling. Functional interpretation was performed using **GO**, **KEGG**, and phenotype enrichment analyses.

---

##  Methods

### **1. Data Acquisition**
- RNA-seq data downloaded from **NCBI SRA** (BioProject: *PRJNA1023207*).
- Samples were grouped into:
  - Male Control
  - Female Control
  - Male AD
  - Female AD

---

### **2. Preprocessing Pipeline (Nextflow)**

The automated pipeline includes:

- **FastQC** — Quality control  
- **Trimmomatic** — Adapter and quality trimming  
- **HISAT2** — Genome alignment  
- **featureCounts** — Gene-level read quantification  

The pipeline organizes all intermediate and final files for seamless downstream analysis.

---

### **3. Differential Expression Analysis (DESeq2)**

- Analysis performed in **R** using **DESeq2**
- Pairwise contrasts:
  - **AD vs Control**
  - **Male vs Female**
- Significance thresholds:
  - `adjusted p-value < 0.05`
  - `|log2FoldChange| > 1`

---

##  Key Findings

### **1. AD vs Control**

#### **MAP3K10 — Upregulated**
- Activates **MAPK/JNK stress pathway**
- Linked with inflammation, neuronal damage  
- KEGG: *Neurodegeneration pathways*  
- MGI phenotypes: altered cholesterol levels, reduced bone mineral content  

#### **IRAG1 — Downregulated**
- Involved in **cGMP-mediated signaling**
- KEGG: *cGMP-PKG signaling*, *Vascular smooth muscle contraction*  
- Suggests impaired neurovascular coupling and intracellular Ca²⁺ handling  

#### **KIF5C-AS (antisense)** — Downregulated
- May influence axonal transport  
- Potential downstream effects on neuronal connectivity and synaptic function  

---

### **2. Sex-Specific Differences (Male vs Female)**

#### **KDM5D — Male-biased**
- Y-linked chromatin modifier  
- Involved in:
  - H3K4 demethylation  
  - Androgen receptor signaling  
- Implicates sex-dependent epigenetic regulation in AD progression  

---

##  Biological Interpretation Summary
Functional annotations point to:

- Altered **second-messenger/vascular signaling** (IRAG1)  
- Activation of **MAPK/JNK stress and neurodegeneration pathways** (MAP3K10)  
- **Sex-specific chromatin remodeling** (KDM5D)

These mechanisms align with known AD biology and provide potential biomarkers or mechanistic targets.

---

##  Repository Structure
```
├── README.md
├── LICENSE
├── nextflow_pipeline/
├── scripts/
│   ├── deseq2_analysis.R
│   ├── plotting.R
├── results/
│   ├── QC/
│   ├── alignments/
│   ├── counts/
│   ├── DESeq2_results/
└── figures/
```

---

##  Technologies Used
- Nextflow  
- FastQC  
- Trimmomatic  
- HISAT2  
- featureCounts  
- R / DESeq2  
- KEGG, GO, phenotype enrichment tools  

---

##  Contact
**Name:**PRASHANTH E
**Email:** prashantheprashanth584@gmail.com 

---

## 📄 License
This project is released under the **MIT License**.

