# FDA-HFP Software Projects

Welcome to the software portfolio of the FDA Human Foods Program Biostatistics group. This organization develops and maintains bioinformatics tools and pipelines for genomic surveillance, pathogen detection, and food safety applications.

## Overview

The Biostatistics team develops computational tools for analyzing whole genome sequencing (WGS) data from foodborne pathogens. Our software supports outbreak investigation, pathogen characterization, and regulatory decision-making to protect public health. The tools range from SNP analysis pipelines to serotyping tools, metagenomic classifiers, and workflow frameworks.

---

## Core SNP Analysis Tools

### [SNP Pipeline](https://github.com/CFSAN-Biostatistics/snp-pipeline)
The flagship CFSAN SNP Pipeline performs reference-based alignment to generate SNP matrices from NGS data for phylogenetic analysis of closely-related pathogenic organisms. Designed specifically for outbreak investigations, it automates alignment, SNP position identification, and SNP calling to produce high-quality matrices for food safety applications. Published in PeerJ Computer Science (2015), this Python-based pipeline remains one of the most widely adopted tools for foodborne pathogen surveillance.

### [CSP2](https://github.com/CFSAN-Biostatistics/CSP2)
CSP2 (CFSAN SNP Pipeline 2) is a Nextflow-based pipeline that provides fast and accurate SNP distance estimation from either WGS read data or genome assemblies. It represents the next generation of the original SNP Pipeline with improved performance and workflow management. The tool streamlines comparative genomics for outbreak cluster detection and strain relatedness analysis.

### [SNP Mutator](https://github.com/CFSAN-Biostatistics/snp-mutator)
SNP Mutator generates mutated sequence files from a reference genome, enabling creation of synthetic datasets for validation and benchmarking of variant calling pipelines. This utility tool supports quality control and method development by producing controlled test data with known mutations. It is particularly useful for evaluating the accuracy and sensitivity of SNP detection algorithms.

---

## Pathogen Typing and Serotyping

### [ShigaTyper](https://github.com/CFSAN-Biostatistics/shigatyper)
ShigaTyper performs rapid in silico serotyping of Shigella species from Illumina or Oxford Nanopore sequencing data with minimal computational requirements. The tool identifies serotypes and detects the ipaB virulence marker to assess invasion plasmid presence. Published in Applied and Environmental Microbiology (2019), it improves serotyping accuracy over traditional methods.

### [BetterCallSal](https://github.com/CFSAN-Biostatistics/bettercallsal)
BetterCallSal is a Nextflow workflow that assigns Salmonella serotypes based on genome similarity using multiple k-mer based methods (MASH, SOURMASH, and KMA). Designed for both metagenomic and quasi-metagenomic applications, it enables rapid serotype prediction from complex sample types. The pipeline supports high-throughput analysis in NCBI Pathogen Detection workflows.

### [SeroTools](https://github.com/CFSAN-Biostatistics/SeroTools)
SeroTools provides a comprehensive toolkit and repository for the White-Kauffmann-Le Minor scheme, the standard nomenclature system for Salmonella serotyping. This Python package facilitates consistent serotype assignment and interpretation across laboratory and bioinformatics workflows. It serves as a reference implementation for standardized Salmonella classification.

### [Cronology](https://github.com/CFSAN-Biostatistics/cronology)
Cronology is an automated Nextflow workflow for Cronobacter whole genome sequence assembly, subtyping, and isolate clustering based on the NCBI Pathogen Detection framework. The pipeline performs MLST typing and clusters isolates for outbreak detection in this emerging foodborne pathogen. It streamlines surveillance of Cronobacter species associated with infant formula contamination.

---

## Wastewater and Environmental Surveillance

### [C-WAP](https://github.com/CFSAN-Biostatistics/C-WAP)
The CFSAN Wastewater Analysis Pipeline (C-WAP) analyzes SARS-CoV-2 variants in wastewater samples using reference-based alignment and multiple variant detection methods including Kallisto, Freyja, and Kraken2/Bracken. This Nextflow pipeline processes both short-read (Illumina) and long-read (ONT/PacBio) sequencing data to estimate variant composition in environmental samples. Note: As of June 2023, C-WAP is no longer under active development; users are directed to the successor project Aquascope at CDC.

### [WW Simulations](https://github.com/CFSAN-Biostatistics/ww_simulations)
This repository contains simulation tools and analyses for wastewater surveillance applications. The tools support method development and validation for wastewater-based epidemiology. These simulations help evaluate the performance of variant detection algorithms under different sampling conditions.

### [WW-SC2-variant-estimations](https://github.com/CFSAN-Biostatistics/WW-SC2-variant-estimations)
This project provides methods for estimating SARS-CoV-2 variant proportions from wastewater sequencing data. It complements C-WAP with specialized statistical approaches for variant quantification. The tools support public health decision-making during the COVID-19 pandemic.

---

## Metagenomics and Taxonomic Classification

### [nowayout](https://github.com/CFSAN-Biostatistics/nowayout)
nowayout is an ultra-fast automated Nextflow pipeline for taxonomic classification of eukaryotic mitochondrial reads from metagenomic samples. Leveraging mitochondrial genome markers, it enables rapid species identification in complex food matrices and environmental samples. The tool addresses food fraud detection and allergen screening applications.

### [centriflaken](https://github.com/CFSAN-Biostatistics/centriflaken)
centriflaken is a Nextflow pipeline for precision metagenomics, focusing on high-resolution taxonomic profiling of complex microbial communities. The workflow integrates multiple classification approaches for accurate species-level identification in food and environmental samples. It supports both targeted and exploratory metagenomic analyses.

### [strainfish](https://github.com/CFSAN-Biostatistics/strainfish)
strainfish implements a weighted ensemble machine learning algorithm with multiple DNA sequence encoders specifically designed for classification of marker sequences. The tool combines various sequence representation methods to improve strain-level discrimination. It addresses the challenge of fine-scale taxonomic resolution in genomic surveillance.

---

## Data Processing Utilities

### [VCFtoolz](https://github.com/CFSAN-Biostatistics/vcftoolz)
VCFtoolz provides a collection of utilities for working with Variant Call Format (VCF) files, including filtering, merging, and format conversion operations. The Python package simplifies common VCF manipulation tasks in variant analysis pipelines. It supports quality control and data preparation for downstream phylogenetic analyses.

### [fastatools](https://github.com/CFSAN-Biostatistics/fastatools)
fastatools offers utilities for manipulating FASTA format sequence files, including extraction, filtering, and format conversion operations. These command-line tools streamline sequence data preparation and quality control. The package integrates seamlessly into bioinformatics workflows requiring FASTA file manipulation.

### [refchooser](https://github.com/CFSAN-Biostatistics/refchooser)
refchooser assists in selecting an optimal reference genome from a list of candidate assemblies for comparative genomics analyses. The tool evaluates assembly quality metrics and genomic similarity to identify the most suitable reference. This facilitates consistent and reproducible reference-based analyses across sample sets.

### [table-ops](https://github.com/CFSAN-Biostatistics/table-ops)
table-ops provides utilities for common tabular data operations in bioinformatics workflows. The Python package enables efficient manipulation, filtering, and transformation of delimited text files. It simplifies data wrangling tasks in analysis pipelines.

### [snps_in_genes](https://github.com/CFSAN-Biostatistics/snps_in_genes)
This tool identifies which genes contain SNPs detected by the CFSAN SNP Pipeline, enabling functional annotation of variants. By mapping SNP positions to gene coordinates, it facilitates interpretation of genomic variation in biological context. This supports investigation of virulence factors and antimicrobial resistance mechanisms.

---

## Machine Learning and Variant Calling

### [kyos](https://github.com/CFSAN-Biostatistics/kyos)
kyos provides tools for haploid variant calling using Deep Neural Networks, offering an alternative to traditional statistical variant callers. The machine learning approach leverages patterns in read alignment data to improve variant detection accuracy. This tool represents an exploration of AI-driven methods for genomic analysis.

---

## Oxford Nanopore Long-Read Tools

### [porerefiner](https://github.com/CFSAN-Biostatistics/porerefiner)
porerefiner is a management tool for Oxford Nanopore sequencing data, facilitating organization and quality control of long-read datasets. The Python package automates common tasks in ONT data handling and preprocessing. It streamlines workflows for laboratories adopting long-read sequencing technologies.

---

## High-Performance Computing Utilities

### [jobrunner](https://github.com/CFSAN-Biostatistics/jobrunner)
jobrunner provides an abstraction layer for executing computational jobs on HPC clusters using Grid Engine, Torque, or local execution. This Python package simplifies deployment of bioinformatics pipelines across different computing environments. It enables portable workflow implementations that adapt to available infrastructure.

### [qarrayrun](https://github.com/CFSAN-Biostatistics/qarrayrun)
qarrayrun is a helper tool for running array jobs on HPC computational nodes, facilitating parallel execution of embarrassingly parallel tasks. The utility simplifies batch processing of multiple samples in genomic analysis pipelines. It optimizes resource utilization on cluster computing systems.

---

## Workflow Development and Integration

### [wdl-commons](https://github.com/CFSAN-Biostatistics/wdl-commons)
wdl-commons provides a library of common WDL (Workflow Description Language) task components for bioinformatics workflows. These reusable modules standardize common operations across different pipelines. The repository promotes workflow portability and reproducibility in WDL-based analyses.

### [microrunqc-wdl](https://github.com/CFSAN-Biostatistics/microrunqc-wdl)
microrunqc-wdl implements quality control workflows for microbial sequencing runs in WDL format. The workflow automates assessment of sequencing quality metrics and generates standardized QC reports. It supports laboratory quality management and data release decisions.

---

## Galaxy Platform Tools

### [Seqsero2_galaxy](https://github.com/CFSAN-Biostatistics/Seqsero2_galaxy)
This repository provides a Galaxy wrapper for the SeqSero2 Salmonella serotyping tool developed by Zhang et al. The wrapper enables integration of SeqSero2 into Galaxy workflows for accessible, point-and-click Salmonella serotype prediction. It expands the Galaxy ecosystem for food safety microbiology.

### [amrfinder-galaxy](https://github.com/CFSAN-Biostatistics/amrfinder-galaxy)
amrfinder-galaxy wraps NCBI's AMRFinder Plus tool for antimicrobial resistance gene detection in Galaxy. This integration enables accessible AMR surveillance within the Galaxy platform. It supports regulatory monitoring of antimicrobial resistance in foodborne pathogens.

### [refchooser-galaxy](https://github.com/CFSAN-Biostatistics/refchooser-galaxy)
This Galaxy implementation of refchooser enables reference genome selection through the Galaxy interface. The tool makes reference selection accessible to users without command-line expertise. It streamlines reference-based analysis setup in Galaxy workflows.

### [galaxy-workflows](https://github.com/CFSAN-Biostatistics/galaxy-worflows)
This repository contains Galaxy workflows published by CFSAN for pathogen genomics analyses. The workflows integrate multiple tools into end-to-end pipelines for common analytical tasks. They provide validated, reproducible analysis protocols for the research community.

### [galaxy-pullthrough-cache](https://github.com/CFSAN-Biostatistics/galaxy-pullthrough-cache)
galaxy-pullthrough-cache implements a pull-through container cache for Galaxy, improving performance and reliability of containerized tool execution. This infrastructure component reduces dependency on external container registries. It enhances Galaxy deployment in restricted network environments.

---

## Data Resources and Infrastructure

### [One Health Enteric Package](https://github.com/CFSAN-Biostatistics/One_Health_Enteric_Package)
This repository defines a One Health-compatible metadata package for genomic surveillance of enteric microbes, enabling data integration across human, animal, food, and environmental sources. The standardized metadata schema facilitates multi-sector collaboration in foodborne disease surveillance. It promotes interoperability in global pathogen genomics networks.

### [CSP2_TestData](https://github.com/CFSAN-Biostatistics/CSP2_TestData)
CSP2_TestData provides test datasets for the CSP2 pipeline, enabling validation and benchmarking of installation and configuration. The repository includes example inputs and expected outputs for quality control. It supports reproducible testing of the CSP2 workflow.

### [data-commons](https://github.com/CFSAN-Biostatistics/data-commons)
This repository hosts shared data resources and common datasets used across CFSAN bioinformatics projects. It centralizes reference data and reduces redundancy across multiple pipelines. The commons approach promotes consistency in analytical workflows.

### [ncbi-pathogen-api-docs](https://github.com/CFSAN-Biostatistics/ncbi-pathogen-api-docs)
This repository documents NCBI's undocumented Pathogen Detection API, providing community knowledge about programmatic access to NCBI Pathogen Detection data. The documentation fills gaps in official API coverage and enables advanced data retrieval workflows. It supports integration of NCBI Pathogen Detection resources into custom pipelines.

---

## Method Validation and Quality Control

### [MultiLab-POD-LOD-ICC](https://github.com/CFSAN-Biostatistics/MultiLab-POD-LOD-ICC)
This R Shiny application facilitates interlaboratory microbiological method validation studies, computing probability of detection (POD), limit of detection (LOD), and intraclass correlation coefficients (ICC). The tool standardizes statistical analysis for method validation in food safety microbiology. It supports regulatory acceptance of new testing methods.

### [wgs_competency](https://github.com/CFSAN-Biostatistics/wgs_competency)
This repository provides templates for quantitative and graphical reports summarizing bacterial whole genome DNA sequencing competency assessment analyses. The standardized reporting framework supports laboratory proficiency testing and quality assurance programs. It ensures consistent evaluation of WGS capabilities across testing laboratories.

### [RIPS](https://github.com/CFSAN-Biostatistics/RIPS)
RIPS (Rapid Intuitive Pathogen Surveillance) is an R-based tool for streamlined pathogen surveillance data analysis and visualization. The package provides rapid exploratory analysis of genomic surveillance data for public health applications. It enables intuitive interpretation of complex surveillance datasets.

---

## About HFP Biostatistics

The Biostatistics group develops computational methods and software tools to support the Human Foods Program's mission of protecting public health through science-based regulation.

Our software development follows open-source principles, with most tools available under permissive licenses. We welcome community contributions, bug reports, and feature requests through GitHub issues. For questions about specific tools, please refer to individual repository documentation.

## Citation and Contact

When using these tools in published research, please cite the original publications listed in each repository's documentation. For general inquiries about HFP bioinformatics tools, please open an issue in the relevant repository or contact the FDA Human Foods Program.

---

*Last updated: February 2026*
