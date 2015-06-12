

Attempt at replicating the surname inference attack via genomics of Gymrek et al. 

WARNING WARNING WARNING
I am not a biologist, so this might be wrong. Most likely the corrections to the marker values, as there are many standards.
I did not calibrate against known values, as was done in the first part of the Gymrek paper. Still ongoing, would appreciate corrections.
WARNING WARNING WARNING


Main file: `analyse_str.py` (see doc there)


File `internet_archive_www_smgf_org_slash_ychromosome_slash_marker_standards.jspx` 
is an Internet Archive scrape of old website listing correspondance between genealogical and NIST Y-STR markers, located [here](https://web.archive.org/web/20130518082814/http://www.smgf.org/ychromosome/marker_standards.jspx)


`Science_Gymrek_et_al/` contains Science original paper and supplementary material


`lobSTR_ystr_hg19.bed` is  a .bed reference file for Y-STRs


`1000genomes_20130606_g1k.ped` and `1000genomes_20130606_g1k_CEU.ped` are pedigree files for 1000 genomes participants

Requires `1000Genomes_Phase_1.vcf` file obtained [from Gymrek's website](http://files.teamerlich.org/lobSTR/1000Genomes_Phase_1.vcf.gz)


