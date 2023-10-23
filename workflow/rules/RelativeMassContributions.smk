rule CalculateMassContribution:
    input: 
        Results + 'PeptonizerResults.csv',
        ResultsDir + 'GraphDataframe.csv'
    params:
        PosteriorThreshold = PosteriorThreshold
    output:
        ResultsDir + 'MassContribution.csv'
    shell: 'python3 workflow/scripts/taxa_clustering.py --PeptonizerResults {input[0]} --GraphDataFrame {input[1]} --out {output} --posterior_threshold {params.PosteriorThreshold}'

rule PlotMassContributions:
    input:
        ResultsDir + 'MassContribution.csv'
    output:
        ResultsDir + 'MassContribution.png'
    shell: 'python3 workflow/scripts/plot_mass_contribution.py --MassContribution {input} --out {output}'