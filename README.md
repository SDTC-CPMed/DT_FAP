
# A Familial Adenomatous Polyposis-Derived Signature May Predict Long-Term Sporadic Colorectal Cancer Risk and Reveal Candidates for Early Prevention
 

## Abstract

Key questions in cancer research are how to: 1) systematically characterize and prioritize mechanisms driving malignant transformation; 2) find biomarkers and drug targets for early diagnosis and prevention. Here, we address these challenges in colorectal cancer (CRC) by using familial adenomatous polyposis (FAP) – a Mendelian disease characterized by extensive colorectal polyposis and a high lifetime CRC risk – as a model system.

We constructed multicellular network models (MCNMs) of FAP based on cell-cell communication inferred using MultiNicheNet from scRNA-seq–derived intercellular molecular interactions and defined an FAP-derived signature (FdS) capturing the most influential cellular and molecular interactions. The mechanistic relevance of the FdS was evaluated using multi-omics analyses, including proteomics, spatial transcriptomics, and immunohistochemistry in colorectal tissues from FAP and sporadic CRC. Its potential as a predictive biomarker was assessed in plasma proteomic data from 50,000 UK Biobank participants using machine-learning models, and its therapeutic relevance was explored through drug enrichment analyses and in vitro studies in cancer cell lines.

A 19-gene FdS was identified. Supporting its relevance for malignant transformation, the FdS was enriched in early cancer hallmark pathways. Spatial transcriptomics and immunohistochemistry of more than 400 tissue areas demonstrated that the FdS was associated with malignant transformation in both FAP and CRC. In peripheral blood, plasma measurable FdS proteins were significantly associated with future CRC risk up to 14 years before diagnosis (AUC 0.73-0.78). Drug enrichment analyses across 19,156 drugs identified both clinically used CRC preventive and therapeutic agents, as well as novel candidate drugs or natural compounds with low cost and side effects. In vitro experiments with three of the later showed that cell lines representing early and late tumorigenesis were inhibited by single or combinatorial treatment, respectively.

This study introduces a novel framework based on Mendelian disease modeling to prioritize disease genes, biomarkers, natural compounds, and drugs for potential early diagnosis and prevention of malignant transformation in CRC. Natural compounds may merit further studies for long-term preventive treatments because of their low cost and limited side effects. 

<img width="895" height="726" alt="image" src="https://github.com/user-attachments/assets/00ae94ee-a799-4ac4-b5b8-240883515ad6" />


The code is written in R and Python. The code is organized into the following directories:

R code: R script for the MCNM construction, Immunohistochemistry (IHC), Spatial Transcriptomics, proteomics and machine learning analyses. R (version 4.2.2) 

Python code: Python script for pathway and drug enrichment. Python (version 3.12) 

