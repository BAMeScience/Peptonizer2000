###rules that run ms2rescore

def InputMod(name, unimod_accession, mass_shift, amino_acid, n_term, c_term):
    return '{"name":"'+name+'", "unimod_accession":'+str(unimod_accession)+', "mass_shift":'+str(mass_shift)+', "amino_acid":'+str(amino_acid)+', "n_term":'+str(n_term)+', "c_term":'+str(c_term)+'},'

rule createMS2RescoreConfig:
    input: ExperimentDir+'{spectrum_name}/Xtandem/tandem.output.xml'
    output:
        ExperimentDir+'{spectrum_name}/ms2rescore/config.json'
    run:
        f_out = open(output[0], "w")

        lines = [' {\"$schema\":\"./config_schema.json\",','\"general\":{']

        lines.append('"pipeline":'+RescorePipeline+',')
        lines.append('"feature_sets":'+RescoreFeatures+',')
        lines.append('"run_percolator":'+RunPercolator+',')
        lines.append('"id_decoy_pattern": "DECOY",')
        lines.append('"num_cpu":'+'1'+',')
        lines.append('"config_file":null,')
        lines.append('"tmp_path":null,')
        lines.append('"mgf_path":null,')
        lines.append('"output_filename":null,')
        lines.append('"log_level":\"info\",')
        lines.append('"plotting":false')
        lines.append('},')
        lines.append('\"ms2pip\":{')
        lines.append('\"model\":'+FragModel+',')
        lines.append('\"modifications\":[')
        if mods:
           for mod in mods:
            lines.append(InputMod(mod[0],mod[1],mod[2],mod[3],mod[4],mod[5]))
        lines[-1] = lines[-1][:-1]

        lines.append(']')
        lines.append('},')
        lines.append('\"maxquant_to_rescore\":{},')
        lines.append('\"percolator\":{}')
        lines.append('}')

        f_out.writelines([line + "\n" for line in lines])
        f_out.close()


rule RunMS2Rescore:
    input:
        ExperimentDir+'{spectrum_name}/Xtandem/tandem.output.xml',
        ExperimentDir+'{spectrum_name}/ms2rescore/config.json',
        SpectraDir+'{spectrum_name}.mgf'
    conda: 'envs/graphenv.yml'
    log: ExperimentDir+'{spectrum_name}/ms2rescore/ms2rescore.log'
    params: OutputName = ExperimentDir+'{spectrum_name}/ms2rescore/rescored'
    output: ExperimentDir+'{spectrum_name}/ms2rescore/rescored_searchengine_ms2pip_rt_features.pout'
    shell: 'cp config/config.yaml '+ResultsDir +' && ms2rescore -c {input[1]} -m {input[2]} {input[0]} -o {params.OutputName}'



