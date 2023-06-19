import os
import sys
import matplotlib.pyplot as plt

vcf_file = sys.argv[1] # First argument: VCF input file

VCF = open(vcf_file, 'r') # Open the file in read only mode
reg_size = int(sys.argv[2])# Second argument: size of the groups.
title = sys.argv[3]#Third argument: title of the graph.
line = VCF.readline()# Read first line of the file
variants = {}#Dictionary that will keep the position of the variants in each chromosome(in the tests done in the project there is only one chromosome).
score = {}#Dictionary with the score of each variant in each chromosome.

while line != '':#Read each line until we have reached the end of the file.
    if line[0] != '#':#Removing headers
        div_line = line.split('\t')#Separating each element in the line.
        if div_line[0] not in variants.keys():#Checking if the variant is the dictionaries, if it is not a list with the position of the first variant and score are created in the respective dictionaries.
            variants[div_line[0]] = [div_line[1]]
            score[div_line[0]] = [float(div_line[5])]
        else:
            variants[div_line[0]].append(div_line[1])
            score[div_line[0]].append(float(div_line[5]))
    line = VCF.readline()#Moving to the next line of the vcf file.
VCF.close()#Closing the file once we have finished examining it.

for chrom in variants.keys():#Looping thorugh each chromosome.
    x = variants[chrom]#The x axis of the boxplot will be the variant positions.
    y = score[chrom]#The y axis in the boxplot will be the score of each variant.
    scores = []#Since to create a boxplot the variants must be divided in groups, this list will contain the scores of each group.
    regions = []#This list will contain the region each group covers.
    count = 0#Variable that will count how many memebrs each group has till its reach the maximum allowed.
    region_score = []#This list will contain the scores of the region the next loop is currently iterating, once a group has reached its full capacity it is appended to scores and rest to an empty list.
    region_name = ''#Name of the region, it follows the next format: First position - Last position
    for i in range(0, len(y)):#This loop iterates through each variant dividing them into groups.
        count += 1#Add 1 to the counter
        region_score.append(y[i])#Append the score of the variant to the region list.
        if count == 1:#if it is the first iteration of the group, ad the first position to its name followed by an "-"-
            region_name = x[i] + '-'
        if count >= reg_size or i == (len(y) - 1):#Once a group is full we keep it in scores and region and reset the other variables.
            scores.append(region_score)
            regions.append(region_name + x[i])
            region_score = []
            region_name = ''
            count = 0
    plt.boxplot(scores, labels=regions)# The boxplot is created.
    plt.xticks(rotation=270)
    plt.xlabel('Position')
    plt.ylabel('Phred-scaled confidence')
    plt.title(title)
    plt.show()