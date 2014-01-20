"""
This program parses a ProteinFunction file and convert it into CAFA protein prediction file
"""

import re

class protfun_parser:
    
    def __init__(self):
        # GO term names in Gene Ontology Category and their GO ids
        self.go_terms = {'Signal_transducer':'GO:0004871', 'Receptor':'GO:0004872', 'Hormone':'GO:0005179', 'Structural_protein':'GO:0008018', 
                         'Transporter':'GO:0005215', 'Ion_channel':'GO:0005216', 'Voltage-gated_ion_channel':'GO:0005244', 'Cation_channel':'GO:0005261', 'Transcription':'GO:0000988', 'Transcription_regulation':'GO:0006355', 'Stress_response':'GO:0006950', 
                         'Immune_response':'GO:0006955', 'Growth_factor':'GO:0008083', 'Metal_ion_transport':'GO:0030001'}
        
    def is_protein(self, line):
        """
        A protein name in ProtFun starts with >
        """
        if line[0] == '>':
            return True
        else:
            return False

    def protfun_to_cafa(self, protfun_file):
        prot_file = open(protfun_file)
        cafa_file = open('CAFA_prediction.txt','w')
        protfun_lines = prot_file.readlines()
        go_flag = False
        term_count = 0
        for line in protfun_lines:
            if self.is_protein(line): # Check if line is a protein name
                prot_name = line[1:-1]
            elif 'Gene Ontology category' in line: # Check if line is the start of GO category
                go_flag = True
            elif go_flag == True and term_count < 14: # parse the 14 GO terms and their probs
                go_line = re.split(r'[=> ]+', line)
                cafa_file.write(prot_name+' '+self.go_terms[go_line[1]]+' '+ "%.2f"%float(go_line[2])+'\n')
                term_count += 1
            elif go_flag == True and term_count >= 14: # stop parsing for the current protein
                go_flag = False
                term_count = 0
        
        prot_file.close()
        cafa_file.close()
