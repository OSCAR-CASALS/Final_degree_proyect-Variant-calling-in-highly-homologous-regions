# Variant calling in highly homologous regions
This pipeline is designed to perform variant calling on regions with high homology to a pseudogene by masking all of them with the exeption of one that will represent the rest.

## Required Input
This pipeline requires:
  - A reference genome.
  - A gff3 annotation of the reference genome.
  - A fastq file with the reads we want to do variant calling on.
  
## How to use it?
All the steps the pipeline have been divided in 4 different python scripts that can be run separately or together with **Execute_all.py**. If each script is going to be run separately the next steps must be followed:

  - **STEP_01:** The first script to be run is _Init.sh_ since it will create the genes Blast database and the pseudogenes fasta file required tu run the rest of the pipeline.
  - **STEP_02:** Once the necessary files are created _Create_beds.py_ must be executed in order to make the Align_to and Realign bed files that will be used in the next step for variant calling.
  - **STEP_03:** With the bed files created _Find_variants.py_ must be run in order to create a vcf file with all the mutations detected in regions of high homology.
  - **STEP_04:** Once the last step has finished _Find_hom_zones.py_ can be run to add the 'High_Homology_Region' flag to any variant in a zone that presented high similarity with a pseudogene.

Either the pipeline is run by **Execute_all.py** or executing each step separately all the python and bash scripts in this github must be in the same directory.

## Dependencies
To run this pipeline the following libraries are required:
  - GATK
  - BEDTOOLS
  - BLAST+
  - ARGPARSE
 
 This github already provides a YML file called _Main_env.yml_ with the conda environment that was used to run this script when it was tested.
 
 ## REFERENCES
Pink, Ryan Charles, Kate Wicks, Daniel Paul Caley, Emma Kathleen Punch, Laura Jacobs, and David Raul Francisco Carter. “Pseudogenes: Pseudo-Functional or Key Regulators in Health and Disease?” RNA 17, no. 5 (May 2011): 792–98. https://doi.org/10.1261/rna.2658311.
 
Genome.gov. “Pseudogén,” September 14, 2022. https://www.genome.gov/es/genetics-glossary/Pseudogen. 

Ebbert, Mark T. W., Tanner D. Jensen, Karen Jansen-West, Jonathon P. Sens, Joseph S. Reddy, Perry G. Ridge, John S. K. Kauwe, et al. “Systematic Analysis of Dark and Camouflaged Genes Reveals Disease-Relevant Genes Hiding in Plain Sight.” Genome Biology 20, no. 1 (May 20, 2019): 97. https://doi.org/10.1186/s13059-019-1707-2. 

Ebbert, Mark T. W., Tanner D. Jensen, Karen Jansen-West, Jonathon P. Sens, Joseph S. Reddy, Perry G. Ridge, John S. K. Kauwe, et al. “Systematic Analysis of Dark and Camouflaged Genes Reveals Disease-Relevant Genes Hiding in Plain Sight.” Genome Biology 20, no. 1 (May 20, 2019): 97. https://doi.org/10.1186/s13059-019-1707-2. 

Poliseno, Laura, Andrea Marranci, and Pier Paolo Pandolfi. “Pseudogenes in Human Cancer.” Frontiers in Medicine 2 (September 25, 2015): 68. https://doi.org/10.3389/fmed.2015.00068.

Holley, Guillaume, Doruk Beyter, Helga Ingimundardottir, Peter L. Møller, Snædis Kristmundsdottir, Hannes P. Eggertsson, and Bjarni V. Halldorsson. “Ratatosk: Hybrid Error Correction of Long Reads Enables Accurate Variant Calling and Assembly.” Genome Biology 22, no. 1 (January 8, 2021): 28. https://doi.org/10.1186/s13059-020-02244-4.
Marwaha, Shruti, Joshua W. Knowles, and Euan A. Ashley. “A Guide for the Diagnosis of Rare and Undiagnosed Disease: Beyond the Exome.” Genome Medicine 14, no. 1 (February 28, 2022): 23. https://doi.org/10.1186/s13073-022-01026-w.

Duret, L. (2008) Neutral theory: The null hypothesis of molecular evolution. Nature Education 1(1):218
Ebler, Jana, Peter Ebert, Wayne E. Clarke, Tobias Rausch, Peter A. Audano, Torsten Houwaart, Yafei Mao, et al. “Pangenome-Based Genome Inference Allows Efficient and Accurate Genotyping across a Wide Spectrum of Variant Classes.” Nature Genetics 54, no. 4 (April 2022): 518–25. https://doi.org/10.1038/s41588-022-01043-w. 

Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li. Twelve years of SAMtools and BCFtools. GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008

Danecek, Petr, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, et al. “Twelve Years of SAMtools and BCFtools.” GigaScience 10, no. 2 (February 16, 2021): giab008. https://doi.org/10.1093/gigascience/giab008.

Stern, Adi, and Raul Andino. “Viral Evolution: It Is All About Mutations.” In Viral Pathogenesis: From Basics to Systems Biology: Third Edition, 233–40. Elsevier Inc., 2016. https://doi.org/10.1016/B978-0-12-800964-2.00017-3.

Tanner, Georgette, David R. Westhead, Alastair Droop, and Lucy F. Stead. “Simulation of Heterogeneous Tumour Genomes with HeteroGenesis and in Silico Whole Exome Sequencing.” Bioinformatics (Oxford, England) 35, no. 16 (August 15, 2019): 2850–52. https://doi.org/10.1093/bioinformatics/bty1063.

Blueprint Genetics. “Blueprint Genetics’ Approach to Pseudogenes and Other Duplicated Genomic Regions.” Accessed May 14, 2023. https://blueprintgenetics.com/pseudogene/.

Li, Heng. “Aligning Sequence Reads, Clone Sequences and Assembly Con*gs with BWA-MEM,” 2014, 0 Bytes. https://doi.org/10.6084/M9.FIGSHARE.963153.V1.

McKenna, Aaron, Matthew Hanna, Eric Banks, Andrey Sivachenko, Kristian Cibulskis, Andrew Kernytsky, Kiran Garimella, et al. “The Genome Analysis Toolkit: A MapReduce Framework for Analyzing next-Generation DNA Sequencing Data.” Genome Research 20, no. 9 (September 2010): 1297–1303. https://doi.org/10.1101/gr.107524.110.

Lapidus, A. L. “Genome Sequence Databases: Sequencing and Assembly.” In Encyclopedia of Microbiology (Third Edition), edited by Moselio Schaechter, 196–210. Oxford: Academic Press, 2009. https://doi.org/10.1016/B978-012373944-5.00028-6.

Quinlan, Aaron R., and Ira M. Hall. “BEDTools: A Flexible Suite of Utilities for Comparing Genomic Features.” Bioinformatics 26, no. 6 (March 15, 2010): 841–42. https://doi.org/10.1093/bioinformatics/btq033.

Cock, Peter J. A., Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg, et al. “Biopython: Freely Available Python Tools for Computational Molecular Biology and Bioinformatics.” Bioinformatics 25, no. 11 (June 1, 2009): 1422–23. https://doi.org/10.1093/bioinformatics/btp163.

“Picard Tools - By Broad Institute.” Accessed May 14, 2023. https://broadinstitute.github.io/picard/. 

“Homo Sapiens Chromosome 3, GRCh38 Reference Primary Assembly,” December 20, 2013. 568336021. NCBI Nucleotide Database. http://www.ncbi.nlm.nih.gov/nuccore/CM000665.2. 

“Nh13/DWGSIM: Whole Genome Simulator for Next-Generation Sequencing.” Accessed May 15, 2023. https://github.com/nh13/DWGSIM. 

Lai, Jianbo, Peifen Zhang, Jiajun Jiang, Tingting Mou, Yifan Li, Caixi Xi, Lingling Wu, et al. “New Evidence of Gut Microbiota Involvement in the Neuropathogenesis of Bipolar Depression by TRANK1 Modulation: Joint Clinical and Animal Data.” Frontiers in Immunology 12 (2021): 789647. https://doi.org/10.3389/fimmu.2021.789647. 

Camacho, Christiam, George Coulouris, Vahram Avagyan, Ning Ma, Jason Papadopoulos, Kevin Bealer, and Thomas L. Madden. “BLAST+: Architecture and Applications.” BMC Bioinformatics 10, no. 1 (December 15, 2009): 421. https://doi.org/10.1186/1471-2105-10-421. 

B, Harikrishnan N. “Confusion Matrix, Accuracy, Precision, Recall, F1 Score.” Analytics Vidhya (blog), June 1, 2020. https://medium.com/analytics-vidhya/confusion-matrix-accuracy-precision-recall-f1-score-ade299cf63cd. 

“Phred Quality Score - an Overview | ScienceDirect Topics.” Accessed May 15, 2023. https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/phred-quality-score. 
