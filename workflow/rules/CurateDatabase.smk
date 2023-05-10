#rules to produce files necessary for searching after filtering host spectra or to search all spectra but whithout host or crap DB added
rule RemoveDuplicates:
     input: DatabaseDirectory+ ReferenceDBName+'.fasta'
     output: DatabaseDirectory+ReferenceDBName+'_UNI.fasta'
     conda: 'envs/graphenv.yml'
     shell:  'seqkit rmdup -s {input} > {output}'

rule AddDecoysCrap:
     input: DatabaseDirectory+ReferenceDBName+'_UNI.fasta'
     output: 
          DatabaseDirectory+ReferenceDBName+'_UNI_concatenated_target_decoy.fasta',
          DatabaseDirectory+ReferenceDBName+'_UNI_decoy.fasta'
     conda: 'envs/decoy.yml'
     shell: 'decoypyrat {input} -o {output[1]} -d DECOY && cat {input} {output[1]} > {output[0]}'

# check if spectrum should be filtered or not
def SpectrumToUse(condition):
     if condition:
          return ResultsDir+SampleName+'/SpectraFilter/Filtered_'+HostName+SpectraFileType
     else:
          return DataDirectory+SampleName+'/'+SampleName+SpectraFileType
          
#if the spectra aren't beeing filtered, check whether host and crap should be added to the search DB
def DBToUse(condition):
     if condition:
        return DatabaseDirectory+HostName+'+crap+'+ReferenceDBName+'_UNI_concatenated_target_decoy.fasta'
     else:
        return DatabaseDirectory+ReferenceDBName+'_UNI_concatenated_target_decoy.fasta'

InputSpectrum = SpectrumToUse(FilterSpectra)
InputDB = DBToUse(AddHostandCrapToDB)

