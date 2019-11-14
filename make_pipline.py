# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 11:52:11 2014

@author: mtinti-x
"""
import os

replace_list = eval(open('vars.txt').read())
template_file = open('_template_mapping.sh').read()
template_R = open('_template_count.R').read()
#print(template_file)
if __name__ == '__main__':
    run_all_content = ''
    for dictionary in replace_list:
        sh_script_name = dictionary['base_fastq']+'.sh'
        #run_all_content+='dos2unix '+sh_script_name+'\n'
        run_all_content+='chmod +x '+sh_script_name+'\n'
        run_all_content+='qsub '+sh_script_name+'\n'
        run_all_content+='mv '+sh_script_name+' '+dictionary['experiment']+'\n'
        

        sh_script_content = template_file.format(
                g_version=dictionary['g_version'],
                base_fastq=dictionary['base_fastq'],
                experiment=dictionary['experiment'],
                library=dictionary['library'])
        path_to_bam = os.path.join(dictionary['experiment'],'data',dictionary['base_fastq'])
        file_list = ['.sorted.bam','ff_barcode_sorted.bam','fr_barcode_sorted.bam','rf_barcode_sorted.bam','rr_barcode_sorted.bam']
        r_count = template_R.format(
        count_file = os.path.join(dictionary['experiment'], dictionary['base_fastq']+'count.txt'),
        gtf_file = os.path.join('genomes', dictionary['g_version'], dictionary['g_version']+'.gtf'),
        bam_files = ',\n'.join([ '\"'+os.path.join(path_to_bam, dictionary['base_fastq']+n)+'\"' for n in file_list ])
        )
        open('count_'+dictionary['base_fastq']+'.R','w').write(r_count)
        open(sh_script_name,'w').write(sh_script_content)
        
    run_all_content+='mv '+'run_all_'+dictionary['experiment']+'.sh'+' '+dictionary['experiment']+'\n'
    open('run_all_'+dictionary['experiment']+'.sh','w').write(run_all_content)
    

    
