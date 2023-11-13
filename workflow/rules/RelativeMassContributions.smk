rule CalculateMassContribution:
    input: 
        ResultsDir + 'PeptonizerResults.csv',
        ResultsDir + 'GraphDataframe.csv'
    params:
        PosteriorThreshold = PosteriorThreshold
    conda: 'envs/mass_contrib.yml'
    
    output:
        ResultsDir + 'MassContribution.csv'

    shell: 'python3 workflow/scripts/TaxaQuant.py --PeptonizerResults {input[0]} --GraphdataFrame {input[1]}  --PosteriorThreshold {params.PosteriorThreshold} --out {output}'

rule PlotMassContributions:
    input:
        ResultsDir + 'MassContribution.csv'
    output:
        ResultsDir + 'MassContribution.png'
    conda: 'envs/mass_contrib.yml'
    shell: 'python3 workflow/scripts/plot_mass_contributions.py --GroupedMassResults {input} --out {output}'