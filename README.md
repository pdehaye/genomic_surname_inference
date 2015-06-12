Replication of "Identifying Personal Genomes by Surname Inference"
==================================================================

        Identifying Personal Genomes by Surname Inference
        Melissa Gymrek, Amy L. McGuire, David Golan, Eran Halperin, Yaniv Erlich
        Science 18 January 2013: 
        Vol. 339 no. 6117 pp. 321-324 
        DOI: 10.1126/science.1229566

How much does your genome say about you?
----------------------------------------

For a significant percentage of the male population, your Y chromosome together with genealogical databases pretty much give out your **last name**. This is an attempt at implementing the original attack of Gymrek et al. (or rather reuse the output of their software lobSTR). 

*WARNING WARNING WARNING*

This is an effort to teach myself some computational biology and associated data wrangling. I am neither a biologist nor a computational biologist. I know some of the markers are off, and have not gone through the complete validation tests described in the original paper. Any comment/correction/suggestion welcomed.

*WARNING WARNING WARNING*


Content of this repository
--------------------------

Main file: `analyse_str.py` (see doc there)

File `internet_archive_www_smgf_org_slash_ychromosome_slash_marker_standards.jspx` 
is an Internet Archive scrape of old website listing correspondance between genealogical and NIST Y-STR markers, located [here](https://web.archive.org/web/20130518082814/http://www.smgf.org/ychromosome/marker_standards.jspx)


`Science_Gymrek_et_al/` contains Science original paper and supplementary material


`lobSTR_ystr_hg19.bed` is  a .bed reference file for Y-STRs


`1000genomes_20130606_g1k.ped` and `1000genomes_20130606_g1k_CEU.ped` are pedigree files for 1000 genomes participants

Requires `1000Genomes_Phase_1.vcf` file obtained [from Gymrek's website](http://files.teamerlich.org/lobSTR/1000Genomes_Phase_1.vcf.gz)


