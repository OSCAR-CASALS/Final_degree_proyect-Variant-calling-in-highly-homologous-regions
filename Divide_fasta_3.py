import os
from Bio import pairwise2
class DivideFasta:
    
    def __init__(self, seq_path):
        '''
        This function puts all sequences that can be found in a fasta file inside a dictionary so they can be easily extracted.
        '''
        self.seqs = {} #Dictionary that contains all sequences of a fasta
        file = open(seq_path, "r") # Opening fasta file
        #Iterating sequence by sequence and keeping them in the dictionary.
        line = file.readline()
        chr = ''
        while line != '':
            if line[0] == '>':
                chr = line[1:-1]
                self.seqs[chr] = ''
            else:
                self.seqs[chr] += line
            line = file.readline()
        file.close() #Close file
    
    def extract_region(self, reg_name):
        return self.seqs[reg_name][0:-1]
    
    def get_seq_names(self):
        return list(self.seqs.keys())


def Get_best_alignment_blastn(seq, fasta_div_obj, Temp_file, blast_database, max_percentatge, bed_output_real, to_mask_bed, max_ploidy = -1, Gene_reg = ''):
    '''
    This function extracts the sequence of the pseudogene it has received as input from a DivideFasta Object and blasts it agains a blast databse.
    Those matches with a percentage of sequence identity higher than 90% and lower than the threshold are kept in two bed files:
        - The bed_output_real file will contain all the matches found.
        - The to_mask_bed contains the gene that will represent the other high homology regions found in a posterior variant calling.
    '''
    # Create the file that will contain the results of performing a blast between the pseudogene sequence and the blast database.
    alignment_temp_fasta_file = Temp_file + '/' + seq + '.txt' 
    temporary_alignment = open(alignment_temp_fasta_file, 'w')
    #Extract the pseudogene sequence from the DivideFasta object.
    print('>' + seq, file = temporary_alignment)
    print(fasta_div_obj.extract_region(seq), file = temporary_alignment)
    temporary_alignment.close()
    #The Blast is performed.
    os.system('blastn -db ' + blast_database + ' -query ' + alignment_temp_fasta_file + ' -out ' + Temp_file + '/sample_al.txt -perc_identity 90.0 -outfmt "6 qseqid sseqid qstart qend sstart send evalue pident"')
    al_file = open(Temp_file + '/sample_al.txt', 'r')
    line = al_file.readline()
    ploidy = 2 # A value that will keep track of the ploidy is initialised.
    res = [] # List with all the regions that meet our conditions.
    al_region = [] # If the region specified is found, keep it here so it ca later be inputed in the Align_to file.
    particular_region_al_to = False # If the region specified is found this boolean will be true.
    #I filter the regions in the alignment that are below the percentage of identity threshold and errase compute the ploidy, if it is bigger than the established threshold and the region of interst is not in the alignment, the regions homologous to the pseudogene will not be taken into account.
    while line != '':
        divided_l = line.split('\t')
        identity = divided_l[-1]
        if float(identity) <= max_percentatge:
            pseudo_pos = divided_l[0].split(':')[1].split('-')
            name = divided_l[1].split(':')
            ploidy += 2 # For each high homology region found I add 2 to the ploidy
            pos = name[1].split('-')
            res.append([name[0], pos[0], pos[1], divided_l[0], int(pseudo_pos[0]) + int(divided_l[2]), int(pseudo_pos[0]) + int(divided_l[3]), int(pos[0]) + int(divided_l[4]), int(pos[0]) + int(divided_l[5]), divided_l[-1][0:-1]])
            if divided_l[1] == Gene_reg: # If the region of interest is in the alignment keep it in the al_region list so it represents the other regions homologous to it in a posterior variant calling.
                al_region = res[-1]
                particular_region_al_to = True
        line = al_file.readline()
        if ploidy > max_ploidy and ploidy != -1:
            if particular_region_al_to == True: # If the ploidy threshold has been surpased but the region of interest has been found we make an exeption for this pseudogene and take into account the higher ploidy bed.
                break
            else:
                al_file.close()
                return 0
    al_file.close()
    if ploidy > 2:
        #The Align_to file is created.
        output_file_align_to = to_mask_bed + '/align_to_ploidy_' + str(ploidy) + '.bed'
        Masked_file = open(output_file_align_to, 'a')
        r = []
        if al_region == []: #If the region specified has not been found we keep the last one found as the representative of the camouflaged group in a posterior variant calling, otherwise use the region found instead.
            r = res[-1]
        else:
            r = al_region
        print(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], sep = '\t', file = Masked_file)
        Masked_file.close()
        #The realign bed file is created and all the regions found are kept.
        re_to_f = bed_output_real + '/realign_ploidy_' + str(ploidy) + '.bed'
        re_to = open(re_to_f, 'a')
        for i in res:
            print(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], sep = '\t', file = re_to)
        re_to.close()
    return 0

def Flagg_high_homology_regions(vcf_file, hom_regs, output):
    '''
    This function locates the regions given as input in a vcf and flags them as regions with high homology.
    '''
    vcf = open(vcf_file, 'r')
    res = open(output, 'w')
    line = vcf.readline()
    while line != '':
        if line[0] != '#':
            l = line.split('\t')
            pos = int(l[1])
            for region in hom_regs:
                if pos >= region[0] and pos <= region[1]: #If the mutation in the vcf file is in one of the regions in the list, mark it as high homology zone.
                    line += ('\t' + 'High_Homology_Region')
        print(line, file = res)
        line = vcf.readline()
    res.close()
    vcf.close()
    

