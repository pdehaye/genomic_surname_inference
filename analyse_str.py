#!/usr/bin/python
"""
This file processes the output of a lobSTR run already intersected with the regular markers used in genealogy. 
It ultimately gives Y-STR counts for as many markers as possible, in the correct nomenclature.
This nomenclature is complicated, and the original Science Supplementary Materials link to a website that does not exist anymore. 
For this reason, I have decided to include here an archived copy of 
 http://www.smgf.org/ychromosome/marker_standards.jspx
at internet_archive_www_smgf_org_slash_ychromosome_slash_marker_standards.jspx

It gives links to corresponding ysearch (click and fill captcha)
It also gives entries to put in yhrd.org to perform similar search

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


VERBOSE = False
SELECT_CEU = False     # Only select CEU individuals (Utah)

ysearch_request_ordering = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 42, 64, 65, 66, 67, 
68, 69, 70, 71, 49, 72, 73, 51, 74, 75, 76, 77, 78, 79, 80, 43, 44, 45, 46, 47, 48, 50, 52, 53, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100]
# The URL paraemters ordering for ysearch (not used elsewhere, but useful to present in same order as on website)

ysearch_form_data_parsed = [('DYS459a', 15, [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]), ('DYS462', 44, [8, 9, 10, 11, 12, 13, 14]), ('DYS505', 86, [10, 11, 12, 13, 14]), ('DYS448', 21, [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]), ('DYS607', 35, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]), ('DYS464b', 24, [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), ('DYS531', 54, [7, 8, 9, 10, 11, 12, 13, 14]), ('DYS450', 71, [6, 7, 8, 9, 10]), ('DYS494', 84, [8, 9, 10, 11]), ('DYS485', 83, [12, 13, 14, 15, 16, 17, 18]), ('DXYS156-Y', 100, [6, 7, 8, 9, 10, 11, 12, 13]), ('DYS534', 70, [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]), ('DYS439', 10, [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]), ('DYS395S1a', 56, [12, 13, 14, 15, 16, 17, 18]), ('DYS617', 74, [9, 10, 11, 12, 13, 14, 15, 16]), ('DYS511', 63, [7, 8, 9, 10, 11, 12, 13]), ('DYS437', 20, [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]), ('YCAIIa', 32, [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]), ('DYS413b', 65, [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]), ('DYS438', 41, [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]), ('DYS464e', 27, [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), ('DYS464a', 23, [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), ('DYS717', 98, [16, 17, 18, 19, 20, 21, 22]), ('DYS640', 78, [8, 9, 10, 11, 12, 13, 14]), ('DYS392', 12, [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), ('DYS492', 79, [9, 10, 11, 12, 13, 14, 15]), ('DYS393', 1, [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]), ('DYS449', 22, [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41]), ('DYS464d', 26, [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), ('DYS495', 85, [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]), ('DYS452', 52, [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38]), ('DYS454', 18, [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]), ('DYS572', 77, [8, 9, 10, 11, 12, 13]), ('YCAIIb', 33, [17, 18, 19, 20, 21, 22, 23, 24, 25, 26]), ('DYS716', 97, [24, 25, 26, 27, 28, 29]), ('DYS472', 61, [5, 6, 7, 8, 9]), ('DYS636', 93, [10, 11, 12, 13]), ('DYS426', 8, [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]), ('DYS576', 36, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]), ('DYS556', 90, [9, 10, 11, 12, 13, 14]), ('DYS575', 91, [8, 9, 10, 11]), ('DYS445', 50, [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]), ('DYS19/DYS394', 3, [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), ('DYS455', 17, [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]), ('DYS391', 5, [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]), ('DYS19b', 4, [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]), ('DYS385a', 6, [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]), ('DYS490', 69, [8, 9, 10, 11, 12, 13, 14, 15, 16]), ('DYS442', 40, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]), ('DYS557', 66, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), ('DYS456', 34, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]), ('DYS590', 58, [6, 7, 8, 9]), ('DYS434', 81, [8, 9, 10, 11]), ('DYS635', 46, [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28]), ('DYS537', 59, [8, 9, 10, 11, 12, 13]), ('DYS425', 42, [8, 9, 10, 11, 12, 13, 14, 15]), ('DYS578', 55, [5, 6, 7, 8, 9, 10, 11]), ('DYS447', 19, [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]), ('DYS463', 53, [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]), ('DYS458', 14, [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]), ('DYS413a', 64, [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]), ('DYS444', 49, [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]), ('DYS464c', 25, [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), ('DYS385b', 7, [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]), ('DYS435', 82, [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]), ('DYS460', 30, [7, 8, 9, 10, 11, 12, 13]), ('DYS390', 2, [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]), ('DYS481', 72, [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]), ('DYS638', 94, [9, 10, 11, 12]), ('DYS395S1b', 57, [12, 13, 14, 15, 16, 17, 18, 19]), ('DYS406S1', 62, [7, 8, 9, 10, 11, 12, 13, 14]), ('DYS389II', 13, [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]), ('DYS461', 43, [9, 10, 11, 12, 13, 14]), ('DYS568', 75, [7, 8, 9, 10, 11, 12, 13, 14]), ('GATA-A10', 45, [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]), ('CDYb', 39, [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]), ('DYS726', 99, [11, 12, 13, 14, 15, 16]), ('DYS441', 48, [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]), ('DYS643', 95, [8, 9, 10, 11, 12, 13, 14]), ('DYS389I', 11, [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]), ('GATA-H4', 31, [7, 8, 9, 10, 11, 12, 13, 14]), ('DYS565', 80, [9, 10, 11, 12, 13, 14]), ('DYS446', 51, [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]), ('DYS549', 89, [10, 11, 12, 13, 14]), ('DYS464f', 28, [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), ('DYS594', 67, [8, 9, 10, 11, 12, 13]), ('DYS388', 9, [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]), ('DYS464g', 29, [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), ('DYS522', 87, [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]), ('DYS589', 92, [6, 7, 8, 9, 10, 11, 12, 13, 14, 15]), ('CDYa', 38, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]), ('DYS533', 88, [9, 10, 11, 12, 13, 14]), ('DYS520', 73, [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]), ('DYS487', 76, [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]), ('DYS641', 60, [6, 7, 8, 9, 10, 11, 12]), ('DYS714', 96, [21, 22, 23, 24, 25, 26, 27, 28, 29]), ('GAAT1B07', 47, [7, 8, 9, 10, 11, 12, 13]), ('DYS459b', 16, [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]), ('DYS570', 37, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]), ('DYS436', 68, [8, 9, 10, 11, 12, 13, 14, 15, 16])]
# Parsed from the ysearch form, w/ cosmetic changes to make marker names match lobSTR_ystr_hg19.bed

ysearch_form_data = [(1, 'DYS393', [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]), (2, 'DYS390', [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]), (3, 'DYS19/DYS394', [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), (4, 'DYS19b', [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]), (5, 'DYS391', [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]), (6, 'DYS385a', [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]), (7, 'DYS385b', [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]), (8, 'DYS426', [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]), (9, 'DYS388', [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]), (10, 'DYS439', [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]), (11, 'DYS389I', [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]), (12, 'DYS392', [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), (13, 'DYS389II', [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]), (14, 'DYS458', [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]), (15, 'DYS459a', [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]), (16, 'DYS459b', [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]), (17, 'DYS455', [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]), (18, 'DYS454', [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]), (19, 'DYS447', [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]), (20, 'DYS437', [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]), (21, 'DYS448', [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]), (22, 'DYS449', [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41]), (23, 'DYS464a', [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), (24, 'DYS464b', [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), (25, 'DYS464c', [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), (26, 'DYS464d', [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), (27, 'DYS464e', [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), (28, 'DYS464f', [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), (29, 'DYS464g', [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), (30, 'DYS460', [7, 8, 9, 10, 11, 12, 13]), (31, 'GATA-H4', [7, 8, 9, 10, 11, 12, 13, 14]), (32, 'YCAIIa', [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]), (33, 'YCAIIb', [17, 18, 19, 20, 21, 22, 23, 24, 25, 26]), (34, 'DYS456', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]), (35, 'DYS607', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]), (36, 'DYS576', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]), (37, 'DYS570', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]), (38, 'CDYa', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]), (39, 'CDYb', [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]), (40, 'DYS442', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]), (41, 'DYS438', [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]), (54, 'DYS531', [7, 8, 9, 10, 11, 12, 13, 14]), (55, 'DYS578', [5, 6, 7, 8, 9, 10, 11]), (56, 'DYS395S1a', [12, 13, 14, 15, 16, 17, 18]), (57, 'DYS395S1b', [12, 13, 14, 15, 16, 17, 18, 19]), (58, 'DYS590', [6, 7, 8, 9]), (59, 'DYS537', [8, 9, 10, 11, 12, 13]), (60, 'DYS641', [6, 7, 8, 9, 10, 11, 12]), (61, 'DYS472', [5, 6, 7, 8, 9]), (62, 'DYS406S1', [7, 8, 9, 10, 11, 12, 13, 14]), (63, 'DYS511', [7, 8, 9, 10, 11, 12, 13]), (42, 'DYS425', [8, 9, 10, 11, 12, 13, 14, 15]), (64, 'DYS413a', [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]), (65, 'DYS413b', [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]), (66, 'DYS557', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]), (67, 'DYS594', [8, 9, 10, 11, 12, 13]), (68, 'DYS436', [8, 9, 10, 11, 12, 13, 14, 15, 16]), (69, 'DYS490', [8, 9, 10, 11, 12, 13, 14, 15, 16]), (70, 'DYS534', [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]), (71, 'DYS450', [6, 7, 8, 9, 10]), (49, 'DYS444', [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]), (72, 'DYS481', [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]), (73, 'DYS520', [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]), (51, 'DYS446', [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]), (74, 'DYS617', [9, 10, 11, 12, 13, 14, 15, 16]), (75, 'DYS568', [7, 8, 9, 10, 11, 12, 13, 14]), (76, 'DYS487', [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]), (77, 'DYS572', [8, 9, 10, 11, 12, 13]), (78, 'DYS640', [8, 9, 10, 11, 12, 13, 14]), (79, 'DYS492', [9, 10, 11, 12, 13, 14, 15]), (80, 'DYS565', [9, 10, 11, 12, 13, 14]), (43, 'DYS461', [9, 10, 11, 12, 13, 14]), (44, 'DYS462', [8, 9, 10, 11, 12, 13, 14]), (45, 'GATA-A10', [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]), (46, 'DYS635', [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28]), (47, 'GAAT1B07', [7, 8, 9, 10, 11, 12, 13]), (48, 'DYS441', [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]), (50, 'DYS445', [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]), (52, 'DYS452', [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38]), (53, 'DYS463', [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]), (81, 'DYS434', [8, 9, 10, 11]), (82, 'DYS435', [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]), (83, 'DYS485', [12, 13, 14, 15, 16, 17, 18]), (84, 'DYS494', [8, 9, 10, 11]), (85, 'DYS495', [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]), (86, 'DYS505', [10, 11, 12, 13, 14]), (87, 'DYS522', [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]), (88, 'DYS533', [9, 10, 11, 12, 13, 14]), (89, 'DYS549', [10, 11, 12, 13, 14]), (90, 'DYS556', [9, 10, 11, 12, 13, 14]), (91, 'DYS575', [8, 9, 10, 11]), (92, 'DYS589', [6, 7, 8, 9, 10, 11, 12, 13, 14, 15]), (93, 'DYS636', [10, 11, 12, 13]), (94, 'DYS638', [9, 10, 11, 12]), (95, 'DYS643', [8, 9, 10, 11, 12, 13, 14]), (96, 'DYS714', [21, 22, 23, 24, 25, 26, 27, 28, 29]), (97, 'DYS716', [24, 25, 26, 27, 28, 29]), (98, 'DYS717', [16, 17, 18, 19, 20, 21, 22]), (99, 'DYS726', [11, 12, 13, 14, 15, 16]), (100, 'DXYS156-Y', [6, 7, 8, 9, 10, 11, 12, 13])]
# Sorted form data: HTML form index, marker name, marker range

ysearch_request_initial = "http://www.ysearch.org/search_search.asp?uid=&freeentry=true"
ysearch_request_final = "&min_markers=8&mismatches_max=0&mismatch_type=sliding&mismatches_sliding_starting_marker=8&recaptcha_challenge_field=03AHJ_Vus-aBXBSTqTsvtzTMXaNFeZpK7UVgSAyAvGVqgXPlQmG5RvMBAnkq3ulYz8su4RgQslzWCwVZHO--t5N4srYiZKDHCWSR2in9AlXvmPd25E-uGGn_Odje5Vli4CsYSRd2ibq4ZWu_G_cnuf6kPZVBHT3W7p0Ds_tX3jnBu4aLq5-Kq6kfRoUqZPhSOgnBgrQ6RrU12njumgyxd6D7nx_saYp2VYyRT_gOh6d_D4aWoAH1se52DjeBfv3ICdEvtnPxvP6ZaO&recaptcha_response_field=217&haplo=&region="
# together with ysearch _request middle, this is the form data
# final always includes an obsolete recaptcha, so this link only gets you to a new one, but with the form filled in


yhrd_yfiler_plus = ["DYS576", "DYS389I", "DYS635", "DYS389II", "DYS627", "DYS460", "DYS458", "DYS19", "GATA-H4", "DYS448", "DYS391", "DYS456", "DYS390", "DYS438", "DYS392", "DYS518", "DYS570", "DYS437", "DYS385a/b", "DYS449", "DYS393", "DYS439", "DYS481", "DYF38751", "DYS533"]
yhrd_minimal = ["DYS19", "DYS389I", "DYS389II", "DYS390", "DYS391", "DYS392", "DYS393", "DYS385a/b"] 
# two different ways to perform a search on yhrd

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

    print 
    print "-"*80
    print "Final marker list"
    for marker in sorted(markers_dict.items()):
        print "        ", marker
    print
