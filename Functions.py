import os
class DivideFasta:
    
    def __init__(self, seq_path):
        '''
        This function puts all sequences that can be found in a fasta file inside a dictionary so they can be easily extracted.
        '''
        self.seqs = {}
        file = open(seq_path, "r")
        line = file.readline()
        chr = ''
        while line != '':
            if line[0] == '>':
                chr = line[1:-1]
                self.seqs[chr] = ''
            else:
                self.seqs[chr] += line
            line = file.readline()
        file.close()
    
    def extract_region(self, reg_name):
        return self.seqs[reg_name][0:-1]
    
    def get_seq_names(self):
        return list(self.seqs.keys())


def Make_beds(seq, fasta_div_obj, Temp_file, blast_database, max_percentatge, bed_output_real, to_mask_bed, max_ploidy = -1):
    '''
    This function is in charge of creating the Align_to and Realign bed files. This is done by blasting a pseudogene sequence against a blast database and 
    keeping the region with less e-value on the align_to file and the rest at the realing bed.
    '''
    alignment_temp_fasta_file = Temp_file + '/' + seq + '.txt'
    temporary_alignment = open(alignment_temp_fasta_file, 'w')
    print('>' + seq, file = temporary_alignment)
    print(fasta_div_obj.extract_region(seq), file = temporary_alignment)
    temporary_alignment.close()
    os.system('blastn -db ' + blast_database + ' -query ' + alignment_temp_fasta_file + ' -out ' + Temp_file + '/sample_al.txt -perc_identity 90.0 -outfmt "6 qseqid sseqid qstart qend sstart send evalue pident"')
    al_file = open(Temp_file + '/sample_al.txt', 'r')
    line = al_file.readline()
    ploidy = 2
    res = []
    while line != '':
        divided_l = line.split('\t')
        identity = divided_l[-1]
        if float(identity) <= max_percentatge:
            pseudo_pos = divided_l[0].split(':')[1].split('-')
            name = divided_l[1].split(':')
            ploidy += 2
            pos = name[1].split('-')
            res.append([name[0], pos[0], pos[1], divided_l[0], int(pseudo_pos[0]) + int(divided_l[2]), int(pseudo_pos[0]) + int(divided_l[3]), int(pos[0]) + int(divided_l[4]), int(pos[0]) + int(divided_l[5]), divided_l[-1][0:-1]])
        line = al_file.readline()
        if ploidy >= max_ploidy and ploidy != -1:
            break
    al_file.close()
    if ploidy > 2:
        #Align_to
        output_file_align_to = to_mask_bed + '/align_to_ploidy_' + str(ploidy) + '.bed'
        Masked_file = open(output_file_align_to, 'a')
        r = res[0]
        print(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], sep = '\t', file = Masked_file)
        Masked_file.close()
        #Realign vers
        re_to_f = bed_output_real + '/realign_ploidy_' + str(ploidy) + '.bed'
        re_to = open(re_to_f, 'a')
        for i in res:
            print(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], sep = '\t', file = re_to)
        re_to.close()
    return 0

def Flagg_high_homology_regions(vcf_file, hom_regs, out):
    '''
    This functtion locates the variants inside high homology regions(hom_regs argument) and adds a flag in them.
    '''
    vcf_int = open(vcf_file, 'r')
    res_int = open(out, 'w')
    line_int = vcf_int.readline()
    while line_int != '':
        if line_int[0] != '#':
            l_int = line_int.split('\t')
            pos_int = int(l_int[1])
            for region_int in hom_regs:
                print('Looking for homologous regions in region:', region_int[0], region_int[1])
                if pos_int >= region_int[0] and pos_int <= region_int[1]:
                    line_int = line_int[0:-1] + ('\t' + 'High_Homology_Region') + '\n'
        print(line_int[0:-1], file = res_int)
        line_int = vcf_int.readline()
    res_int.close()
    vcf_int.close()
    

