# FDA-HFP Software Projects

Welcome to the software portfolio of the FDA Human Foods Program Biostatistics group. This organization develops and maintains bioinformatics tools and pipelines for genomic surveillance, pathogen detection, and food safety applications.

## Overview

The Biostatistics team develops computational tools for analyzing whole genome sequencing (WGS) data from foodborne pathogens. Our software supports outbreak investigation, pathogen characterization, and regulatory decision-making to protect public health. The tools range from SNP analysis pipelines to serotyping tools, metagenomic classifiers, and workflow frameworks.

---

## Under Active Development

These tools are under current active development:

| Tool | Description |
|------|-------------|
| [phraya](https://github.com/CFSAN-Biostatistics/phraya) | Ultra-fast sequence aligner written in Rust |
| [torchbase](https://github.com/CFSAN-Biostatistics/torchbase) | Python framework for microbial typing by reference |
| [MAGGIC](https://github.com/CFSAN-Biostatistics/MAGGIC) | Automated workflow for Metagenome-Assembled Genome generation and refinement |
| [CSP2](https://github.com/CFSAN-Biostatistics/CSP2) | Fast and accurate SNP distance estimation from WGS reads or assemblies |
| [bettercallsal](https://github.com/CFSAN-Biostatistics/bettercallsal) | Nextflow workflow for Salmonella serotyping by genome similarity |
| [nowayout](https://github.com/CFSAN-Biostatistics/nowayout) | Ultra-fast taxonomic classification of eukaryotic mitochondrial reads |
| [data-commons](https://github.com/CFSAN-Biostatistics/data-commons) | Shared data resources used across CFSAN bioinformatics projects |

---

## Platforms

### [GalaxyTrakr](https://galaxytrakr.org)
GalaxyTrakr is an FDA-hosted Galaxy instance providing free access to bioinformatics tools for whole genome sequencing-based genomic surveillance of foodborne pathogens. Operated by the Human Foods Program in partnership with public health laboratories worldwide, it offers point-and-click access to validated workflows for pathogen typing, outbreak investigation, and AMR detection — without requiring command-line expertise. GalaxyTrakr serves as the primary public interface for deploying and using the tools developed by this organization, and is widely used by state, federal, and international public health agencies.

---

## Core SNP Analysis Tools

### [CSP2](https://github.com/CFSAN-Biostatistics/CSP2)
CSP2 (CFSAN SNP Pipeline 2) is a Nextflow-based pipeline that provides fast and accurate SNP distance estimation from either WGS read data or genome assemblies. It represents the next generation of the original SNP Pipeline with improved performance and workflow management. The tool streamlines comparative genomics for outbreak cluster detection and strain relatedness analysis.

### [SNP Pipeline](https://github.com/CFSAN-Biostatistics/snp-pipeline) *(legacy)*
The original CFSAN SNP Pipeline performs reference-based alignment to generate SNP matrices from NGS data for phylogenetic analysis of closely-related pathogenic organisms. Published in PeerJ Computer Science (2015), this Python-based pipeline was widely adopted for foodborne pathogen surveillance. **Note: SNP Pipeline has been superseded by CSP2 and is no longer under active development. New projects should use CSP2.**

### [SNP Mutator](https://github.com/CFSAN-Biostatistics/snp-mutator)
SNP Mutator generates mutated sequence files from a reference genome, enabling creation of synthetic datasets for validation and benchmarking of variant calling pipelines. This utility tool supports quality control and method development by producing controlled test data with known mutations.

---

## Sequence Alignment and Typing

### [phraya](https://github.com/CFSAN-Biostatistics/phraya)
phraya is an ultra-fast sequence aligner implemented in Rust, designed for high-throughput genomic analysis. It offers significant performance improvements for large-scale alignment tasks in pathogen genomics workflows.

### [torchbase](https://github.com/CFSAN-Biostatistics/torchbase)
torchbase is a Python framework for microbial typing by reference, providing a foundation for building reference-based strain characterization workflows. It supports flexible integration of alignment and typing approaches for diverse pathogens.

---

## Pathogen Typing and Serotyping

### [BetterCallSal](https://github.com/CFSAN-Biostatistics/bettercallsal)
BetterCallSal is a Nextflow workflow that assigns Salmonella serotypes based on genome similarity using multiple k-mer based methods (MASH, SOURMASH, and KMA). Designed for both metagenomic and quasi-metagenomic applications, it enables rapid serotype prediction from complex sample types. The pipeline supports high-throughput analysis in NCBI Pathogen Detection workflows.

### [ShigaTyper](https://github.com/CFSAN-Biostatistics/shigatyper)
ShigaTyper performs rapid in silico serotyping of Shigella species from Illumina or Oxford Nanopore sequencing data with minimal computational requirements. The tool identifies serotypes and detects the ipaB virulence marker to assess invasion plasmid presence. Published in Applied and Environmental Microbiology (2019), it improves serotyping accuracy over traditional methods.

### [SeroTools](https://github.com/CFSAN-Biostatistics/SeroTools)
SeroTools provides a comprehensive toolkit and repository for the White-Kauffmann-Le Minor scheme, the standard nomenclature system for Salmonella serotyping. This Python package facilitates consistent serotype assignment and interpretation across laboratory and bioinformatics workflows.

### [Cronology](https://github.com/CFSAN-Biostatistics/cronology)
Cronology is an automated Nextflow workflow for Cronobacter whole genome sequence assembly, subtyping, and isolate clustering based on the NCBI Pathogen Detection framework. It streamlines surveillance of Cronobacter species associated with infant formula contamination.

---

## Metagenomics

### [MAGGIC](https://github.com/CFSAN-Biostatistics/MAGGIC)
MAGGIC is an automated workflow for the generation and refinement of Metagenome-Assembled Genomes (MAGs) from metagenomic sequencing data. It provides a streamlined approach to recovering high-quality genomes from complex microbial communities in food and environmental samples.

### [nowayout](https://github.com/CFSAN-Biostatistics/nowayout)
nowayout is an ultra-fast automated Nextflow pipeline for taxonomic classification of eukaryotic mitochondrial reads from metagenomic samples. Leveraging mitochondrial genome markers, it enables rapid species identification in complex food matrices and environmental samples, supporting food fraud detection and allergen screening.

### [centriflaken](https://github.com/CFSAN-Biostatistics/centriflaken)
centriflaken is a Nextflow pipeline for precision metagenomics, focusing on high-resolution taxonomic profiling of complex microbial communities. The workflow integrates multiple classification approaches for accurate species-level identification in food and environmental samples.

### [strainfish](https://github.com/CFSAN-Biostatistics/strainfish)
strainfish implements a weighted ensemble machine learning algorithm with multiple DNA sequence encoders specifically designed for classification of marker sequences. The tool combines various sequence representation methods to improve strain-level discrimination.

---

## Wastewater and Environmental Surveillance

### [C-WAP](https://github.com/CFSAN-Biostatistics/C-WAP)
The CFSAN Wastewater Analysis Pipeline (C-WAP) analyzes SARS-CoV-2 variants in wastewater samples using reference-based alignment and multiple variant detection methods including Kallisto, Freyja, and Kraken2/Bracken. **Note: As of June 2023, C-WAP is no longer under active development; users are directed to the successor project Aquascope at CDC.**

### [WW Simulations](https://github.com/CFSAN-Biostatistics/ww_simulations)
Simulation tools and analyses for wastewater surveillance applications, supporting method development and validation for wastewater-based epidemiology.

### [WW-SC2-variant-estimations](https://github.com/CFSAN-Biostatistics/WW-SC2-variant-estimations)
Methods for estimating SARS-CoV-2 variant proportions from wastewater sequencing data, complementing C-WAP with specialized statistical approaches for variant quantification.

---

## Machine Learning and Variant Calling

### [kyos](https://github.com/CFSAN-Biostatistics/kyos)
kyos provides tools for haploid variant calling using Deep Neural Networks, offering a machine learning alternative to traditional statistical variant callers. It leverages patterns in read alignment data to improve variant detection accuracy.

---

## Data Processing Utilities

### [VCFtoolz](https://github.com/CFSAN-Biostatistics/vcftoolz)
VCFtoolz provides utilities for working with Variant Call Format (VCF) files, including filtering, merging, and format conversion. The Python package simplifies common VCF manipulation tasks in variant analysis pipelines.

### [fastatools](https://github.com/CFSAN-Biostatistics/fastatools)
fastatools offers utilities for manipulating FASTA format sequence files, including extraction, filtering, and format conversion. These command-line tools streamline sequence data preparation and quality control.

### [refchooser](https://github.com/CFSAN-Biostatistics/refchooser)
refchooser assists in selecting an optimal reference genome from a list of candidate assemblies. The tool evaluates assembly quality metrics and genomic similarity to identify the most suitable reference.

### [table-ops](https://github.com/CFSAN-Biostatistics/table-ops)
table-ops provides utilities for common tabular data operations in bioinformatics workflows, enabling efficient manipulation, filtering, and transformation of delimited text files.

---

## Oxford Nanopore Long-Read Tools

### [porerefiner](https://github.com/CFSAN-Biostatistics/porerefiner)
porerefiner is a management tool for Oxford Nanopore sequencing data, facilitating organization and quality control of long-read datasets. The Python package automates common tasks in ONT data handling and preprocessing.

---

## Workflow Development and Integration

### [wdl-commons](https://github.com/CFSAN-Biostatistics/wdl-commons)
wdl-commons provides a library of common WDL task components for bioinformatics workflows, promoting workflow portability and reproducibility in WDL-based analyses.

### [microrunqc-wdl](https://github.com/CFSAN-Biostatistics/microrunqc-wdl)
microrunqc-wdl implements quality control workflows for microbial sequencing runs in WDL format. The workflow automates assessment of sequencing quality metrics and generates standardized QC reports.

---

## Galaxy Platform Tools

### [galaxytrakr-tools](https://github.com/CFSAN-Biostatistics/galaxytrakr-tools)
Tool definitions for deployment on the GalaxyTrakr platform, maintained as the primary source of validated tools available through [galaxytrakr.org](https://galaxytrakr.org).

### [Seqsero2_galaxy](https://github.com/CFSAN-Biostatistics/Seqsero2_galaxy)
Galaxy wrapper for the SeqSero2 Salmonella serotyping tool, enabling integration into Galaxy workflows for point-and-click Salmonella serotype prediction.

### [amrfinder-galaxy](https://github.com/CFSAN-Biostatistics/amrfinder-galaxy)
Galaxy wrapper for NCBI's AMRFinder Plus tool for antimicrobial resistance gene detection.

### [refchooser-galaxy](https://github.com/CFSAN-Biostatistics/refchooser-galaxy)
Galaxy implementation of refchooser for reference genome selection through the Galaxy interface.

### [galaxy-workflows](https://github.com/CFSAN-Biostatistics/galaxy-worflows)
Galaxy workflows published by CFSAN for pathogen genomics analyses, providing validated end-to-end pipelines for common analytical tasks.

### [galaxy-pullthrough-cache](https://github.com/CFSAN-Biostatistics/galaxy-pullthrough-cache)
A pull-through container cache for Galaxy, improving performance and reliability of containerized tool execution in restricted network environments.

---

## Data Resources and Infrastructure

### [One Health Enteric Package](https://github.com/CFSAN-Biostatistics/One_Health_Enteric_Package)
A One Health-compatible metadata package for genomic surveillance of enteric microbes, enabling data integration across human, animal, food, and environmental sources.

### [data-commons](https://github.com/CFSAN-Biostatistics/data-commons)
Shared data resources and common datasets used across CFSAN bioinformatics projects, promoting consistency in analytical workflows.

### [CSP2_TestData](https://github.com/CFSAN-Biostatistics/CSP2_TestData)
Test datasets for the CSP2 pipeline, enabling validation and benchmarking of installation and configuration.

---

## Method Validation and Quality Control

### [MultiLab-POD-LOD-ICC](https://github.com/CFSAN-Biostatistics/MultiLab-POD-LOD-ICC)
An R Shiny application for interlaboratory microbiological method validation studies, computing probability of detection (POD), limit of detection (LOD), and intraclass correlation coefficients (ICC).

### [WGS-Econ-Benefit](https://github.com/CFSAN-Biostatistics/WGS-Econ-Benefit)
A Shiny app for quantifying the benefit versus cost of whole genome sequencing for three foodborne pathogens.

### [RIPS](https://github.com/CFSAN-Biostatistics/RIPS)
RIPS (Rapid Intuitive Pathogen Surveillance) is an R-based tool for streamlined pathogen surveillance data analysis and visualization.

---

## High-Performance Computing Utilities

### [jobrunner](https://github.com/CFSAN-Biostatistics/jobrunner)
An abstraction layer for executing jobs on HPC clusters using Grid Engine, Torque, or local execution. Enables portable workflow implementations across different computing environments.

### [qarrayrun](https://github.com/CFSAN-Biostatistics/qarrayrun)
A helper tool for running array jobs on HPC computational nodes, facilitating parallel execution of embarrassingly parallel tasks.

---

## About HFP Biostatistics

The Biostatistics group develops computational methods and software tools to support the Human Foods Program's mission of protecting public health through science-based regulation.

Our software development follows open-source principles, with most tools available under permissive licenses. We welcome community contributions, bug reports, and feature requests through GitHub issues.

## Citation and Contact

When using these tools in published research, please cite the original publications listed in each repository's documentation. For general inquiries about HFP bioinformatics tools, please open an issue in the relevant repository or contact the FDA Human Foods Program.

---

*Last updated: June 2026*
