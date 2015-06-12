#!/usr/bin/python
"""
This file processes the output of a lobSTR run already intersected with the regular markers used in genealogy. 
It ultimately gives Y-STR counts for as many markers as possible, in the correct nomenclature.
This nomenclature is complicated, and the original Science Supplementary Materials link to a website that does not exist anymore. 
For this reason, I have decided to include here an archived copy of 
 http://www.smgf.org/ychromosome/marker_standards.jspx
at internet_archive_www_smgf_org_slash_ychromosome_slash_marker_standards.jspx

It gives links to corresponding ysearch (click and fill captcha)
It also gives entries to put in yhrd.org to perform similar search, and a row for the ybase table in the supplementary materials

I was unable to fully replicate the Gymrek et al Surname Inference attack, but one of the databases they used is now shutdown (SMGF)
The procedure here still works to generate profiles, it is just that they are not accessible publicly in genealogical databases.


1000Genomes_Phase_1.vcf comes from Melissa gymrek's website
lobSTR_ystr_hg19.bed comes from Melissa gymrek's website
(could also potentially use files such as ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf ALL.wgs.v2_lobSTR.20130502.microsat.integrated.genotypes.vcf
but couldn't fully understand the format. The x/y notation confused me for the Y chromosome.)

to prepare intersected bed files:   ./intersectBed -a 1000Genomes_Phase_1.vcf -b lobSTR_ystr_hg19.bed -wa -wb -header > intersect.gymrek.ystr
to run: python analyse_str.py intersect.gymrek.ystr
"""

import sys
from copy import deepcopy
from collections import defaultdict
from additional_data import ysearch_request_ordering, ysearch_form_data_parsed, ysearch_form_data, ysearch_request_initial, ysearch_request_final, yhrd_yfiler_plus, yhrd_minimal, ybase_markers, lds_markers, lds_volunteers


VERBOSE = False
SELECT_CEU = True     # Only select CEU individuals (Utah)

filename = sys.argv[1]

with open(filename, "r") as f:
    data = f.read()

lines = data.splitlines()

header = lines[19].split("\t") # Line 19 is the last of comments


def parse_format(x):
    """  Data is in format 
            GT:ALLREADS:AML:DP:GB:PL:Q:STITCH
         or when absent the entry is ./.
         We want to get to GB, itself of shape x/y. We want to return (x,y) as tuple of ints
    """
    x = x.split(":")
    try:
       x = x[4]
    except: 
       return ("?", "?")
    x = x.split("/")
    x = tuple(map(int,x))
    return x


marked_volunteers = defaultdict(list)  # A dict of markers for each volunteer

for line in lines[1:]:
    if line[0] == "#":
        continue
    data = line.split("\t")

    # each line is composed of the join of the original 1000genomes file and the lobSTR reference file
    # the original 1000 genomes file makes up the header, so len(header) helps delimit the two
    # value of len(header) is 1018 for file used, but written this way for ease of reuse

    chromosome = data[0:1]
    allotype = data[3:4]  # reference and alt
    # in between: quality information, format (i.e. GT:ALLREADS:AML:DP:GB:PL:Q:STITCH)

    genotype = data[9:len(header)]
    period = int(data[len(header)+3])
    ref = int(data[len(header) + 4])
    marker = data[len(header) + 5]   

    for i, d in enumerate(data):
        if i >= 9 and i < len(header):
            volunteer = header[i]
	    d_parsed = parse_format(d)
            #print i, d_parsed, period, ref
	    bp_diff = d_parsed[1]
            if type(bp_diff) == int:
                if bp_diff % period == 0:
                   marked_volunteers[volunteer].append((marker, bp_diff/period + ref))  
            # bp_diff counts difference of bp length from reference genome 
            # better be a multiple of the period, and need to add the repetition count already in ref genome

marked_volunteers = sorted(marked_volunteers.items(), key = lambda (volunteer, markers): len(markers), reverse = True) 
# Sorts them according to number of successful markers

def marker_lookup(marker, markers, bounds):
    """
         Checks that marker is in the marker bounds, returns the value of the marker there
    """
    try:
	tmp = str(markers[marker])
        assert int(tmp) in map(int, bounds)
        return tmp
    except:
        if VERBOSE:
            print "Can't find marker: ", marker
	return "0"

def prune_dict(d):
    """ Rejects some STR due to unreliability, according to Grymek supplementary materials in Science"""
    reject = ["GAAT1B07", "DYS726", "CDYa", "CDYb", "DYS425", "DYF371", "DXYS156-Y", "DYS19b", "DYS640", "DYS464a", "DYS464b", "DYS464c", "DYS464d"]
    for mark in reject:
	d.pop(mark, None)

# On her website, Gymrek has the following:
# Markers with known issues

# DYS640 is annotated as having a reference allele as 9. However, some standards may use a reference allele of 11. The full sequence of the STR in hg19 is TAAAAAATAAAAATAAATAAATAAATAAATAAATAAATAAATAAATAAA, and it is unclear whether the first two repeats with extra "A"'s count in some cases. Due to ambiguity, we usually discard this marker.
# DYS425 actually has multiple locations on chrY, and should generally be ignored when using lobSTR.
# The CODIS markers FGA and D21S11 are usually unreliable with 100bp reads since most alleles are too long to be completely spanned.
# Some long markers (e.g. DYS449 and DYS448) have multiple repeat stretches broken up by non repeat regions. These have been broken up into two STRs (e.g. DYS449.1 and DYS449.2). To genotype these, both of those sub-markers need to be called and added together.
# DYS389B.1 and DYS389B.2 are misannotated in the reference bed, and should have reference alleles listed as 17 and 12, repsectively.
# DYS511 is misannotated in the reference bed and should have a reference allele of 9.
# DYS442 has had changes in its nomenclature (http://www.hprg.com/hapest5/page2.html). FamilyTreeDNA reports a value 5 units shorter than NIST.

# THESE RECOMMENDATIONS HAVE NOT BEEN IMPLEMENTED!!!!!

def FamilyTreeDNA_correct(d):
    """ Corrects on some markers, following FaimlyTreeDNA conventions """
    corrections = [("DYS441", -1), ("DYS442", -5), ("GATA-A10", -2), ("GATA-H4", -1)] 
    for marker, correction in corrections:
        if marker in d:
           d[marker] += correction

def print_YHRD(markers_dict):
    print "-"*80
    print "YHRD Yfiler Plus:"
    print "-"*80
    for marker in yhrd_yfiler_plus:
        print marker, "    ", markers_dict.get(marker, " ")
    print "-"*80
    print "YHRD Minimal"
    print "-"*80
    for marker in yhrd_minimal:
        print marker, "    ", markers_dict.get(marker, " ")
    print "-"*80

def ysearch_request(markers_dict):
    ysearch_request_middle = "".join( \
                                      "&L"+str(ysearch_form_index) + "=" + marker_lookup(marker_name, markers_dict, marker_range)\
                                      for (ysearch_form_index,marker_name,marker_range) in ysearch_form_data)
    return ysearch_request_initial + ysearch_request_middle + ysearch_request_final

    
def ybase(volunteer, markers_dict):
    return volunteer + "|  "+"|".join([ marker                               for marker in ybase_markers]) , \
           volunteer + "|  "+"|".join([ str(markers_dict.get(marker, " ? ")) for marker in ybase_markers])
    
for volunteer, markers in marked_volunteers:
    markers_dict = dict(markers)
    if SELECT_CEU and not "NA" in volunteer or not int(volunteer[2:4]) <= 13:
        continue

    print "="*80
    print "="*80
    print
    print "-"*80
    print "1000 genome identifier: ", volunteer, "     ", len(markers_dict), " markers found"
    print "-"*80
    print 
    print "-"*80
    print "Original, uncorrected marker values, all (including unreliable):"
    print "-"*80
    print markers_dict
    print "-"*80
    print
    print 
    print "-"*80
    print "-"*80
    print "UNCORRECTED yhrd values (not sure needs to be corrected!), includes unreliable markers according to Gymrek"
    print_YHRD(markers_dict)    
    print "-"*80
    print "-"*80
    print "Ybase string"
    a, b = ybase(volunteer, markers_dict)
    print a
    print b
    print "-"*80
    print "-"*80

    FamilyTreeDNA_correct(markers_dict)
    print "SHIFTED yhrd values (not sure needs to be shifted!), includes unreliable markers according to Gymrek"
    print_YHRD(markers_dict)    
    print
    print "-"*80
    print "-"*80
    prune_dict(markers_dict)   
    print "Ysearch request, with markers shifted according to FamilyTreeDNA conventions, and unreliable markers removed according to Gymrek's supplementary materials"
    print "-"*80
    print ysearch_request(markers_dict)

    print "-"*80
    print "Final marker list"
    for marker in sorted(markers_dict.items()):
        print "        ", marker
    print
    
for volunteer in lds_volunteers:
    print volunteer, "     ", "|".join(markers_dict.get(marker, "?") for marker in lds_markers)
