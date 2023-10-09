###rule for parsing percolator file

rule UnipeptQuery:
    input: expand(ExperimentDir+'{spectrum_name}/ms2rescore/rescored_searchengine_ms2pip_rt_features.pout',spectrum_name = SpectraNames)
    params: 
          targetTaxa = targetTaxa,
          FDR = FDR
    log: ResultsDir + 'UnipeptResponse.log'
    output: 
          ResultsDir + 'UnipeptResponse.json',
          ResultsDir + 'UnipeptPeptides.json'
    conda: 'envs/Unipeptquery.yml'   
    shell: "python3 workflow/scripts/UnipeptGetTaxonomyfromPout.py --UnipeptResponseFile {output[0]} --pep_out {output[1]} --TaxonomyQuery {params.targetTaxa} --FDR {params.FDR} --PoutFile {input} --logfile {log}" 

rule ParseToUnipeptCSV:
    input: 
          ResultsDir + 'UnipeptResponse.json',
          ResultsDir + 'UnipeptPeptides.json',
          ResourcesDir + 'taxa_peptidome_size.tsv'
          
          
    params: 
      NumberofTaxa = TaxaNumber,
      TaxaRank = TaxaRank
           
    log: ResultsDir + 'ParsetoCSV.log'
    output: 
            ResultsDir + 'GraphDataframe.csv',
            ResultsDir +'TaxaWeights.csv'
    conda: 'envs/graphenv.yml' 
    shell: "python3 workflow/scripts/WeightTaxa.py --UnipeptResponseFile {input[0]} --UnipeptPeptides {input[1]} --out {output[0]} --TaxaWeightFile {output[1]} --NumberOfTaxa {params.NumberofTaxa} --PeptidomeSize {input[2]} --TaxaRank {params.TaxaRank}" 