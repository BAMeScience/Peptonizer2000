import pandas as pd
import json
import os
import shutil



def ExportPeptidesForUnipeptCli(input,out):
#load dictionary from json
    with open(input) as json_file:
        peptide_to_psms = json.load(json_file)
        #print one peptide per line in new txt file as many times as it is in the dictionary sub entry 'key'
        with open(out, 'w') as f:
            for key in peptide_to_psms.keys():
                for i in range(peptide_to_psms[key]['psms']):
                    f.write("%s\n" % key)



ExportPeptidesForUnipeptCli('/home/tanja/Peptonizer2000/Peptonizer2000/results/CAMPI1_SIHUMIx_allbacteria_s8_root_new_rescore/CAMPI_SIHUMIx/UnipeptPeptides.json','/home/tanja/Peptonizer2000/Peptonizer2000/results/CAMPI1_SIHUMIx_allbacteria_s8_root_new_rescore/CAMPI_SIHUMIx/CAMPI1_SIHUMIx_allbacteria_s8_root_new_rescore_UnipeptPeptides.txt')