import requests
import json
import logging
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('--UnipeptResponseFile', type = str, required = True, help = 'path to Unipept response .json file')
parser.add_argument('--TaxonomyQuery', required = True, help = 'taxa to query in Unipept. If querying all taxa, put [1]') 
parser.add_argument('--FDR', type = float, required = True, help = 'min peptide score for the peptide to be included in the search')
parser.add_argument('--PoutFile', type = str, nargs = '+', required = True, help = 'paths to percolator(ms2rescore) Pout files')
parser.add_argument('--UnipeptPeptides', type = str, required = True, help = 'path to file with peptides')
parser.add_argument('--TaxaRank', type = str, required = True, help = 'taxonomic rank at which you want the peptonizer results resolved')
parser.add_argument('--logfile', type = str, required = True, help = 'path to logfile')

args = parser.parse_args()


from typing import Dict, Tuple, List, Any

UNIPEPT_URL = "http://api.unipept.ugent.be"
UNIPEPT_PEPT2FILTERED_ENDPOINT = "/mpa/pept2filtered.json"
UNIPEPT_TAXONOMY_ENDPOINT = "/api/v2/taxonomy.json"


def MS2RescoreOutParser(pout_files, fdr_threshold, decoy_flag):
    '''
    Parses the ms2rescore pout file for peptides, psm numbers and peptide scores
    :param pout_file: str, path to pout file(s)
    :param fdr_threshold: float, fdr threshold below which psms are kept
    :param decoy_flag: str, can be emtpy string, decoy flag in pout file
    :return: dict, peptides:[score,#psms]
    '''

    pep_score = dict()
    pep_psm = dict()
    pep_score_psm = dict()
    
    for pout_file in pout_files:
        with open(pout_file, "r") as f:
            next(f)  # skip header
            for line in f:
                # skip empty lines
                if line.rstrip() == "":
                    continue
                splitted_line = line.rstrip().split("\t")[0:8]
                #assert len(splitted_line) >= 6, "Input file is wrongly formatted. Make sure that the input is a valid .pout file."
                peptide, psm_id, run, colelction, collection , score, q, pep = splitted_line
                if float(q) < fdr_threshold:
                    peptide = re.sub("\[.*?\]", "", peptide)
                    peptide = peptide.split("/")[0]
                    # update pep_psm
                    if peptide not in pep_psm.keys():
                        pep_psm[peptide] = set()
                        pep_psm[peptide].add(psm_id)
                    else:
                        pep_psm[peptide].add(psm_id)
                    # update pep_score
                    if peptide not in pep_score.keys():
                        if float(pep) <0.001:
                            pep_score[peptide] = '0.001'
                        else:
                            pep_score[peptide] = pep          #adjustement necessary to not have 0 and 1 fuck up probability calculations
                    else:
                        if float(pep) <0.001:
                            pep_score[peptide] = '0.001'
                        else:
                            pep_score[peptide] = min(pep,pep_score[peptide])
                    pep_score_psm[peptide] = [pep_score[peptide],len(pep_psm[peptide])]
                
    return pep_score_psm


def get_descendants_for_taxa(target_taxa: List[int], descendants_rank: str) -> List[int]:
    url = UNIPEPT_URL + UNIPEPT_TAXONOMY_ENDPOINT
    all_descendants = set()  # Using a set to avoid duplicates

    # Split the target taxa into batches of 15
    for i in range(0, len(target_taxa), 15):
        batch = target_taxa[i:i+15]

        # Prepare the request payload
        payload = {
            "input": batch,
            "descendants": 'true',
            "descendants_rank": descendants_rank
        }
        print(payload)

        # Perform the HTTP POST request
        response = requests.post(url, json=payload)

        # Check if the request was successful
        if response.status_code == 200:
            print(response)
            data = response.json()
            print('taxa in filter are' + str(data))

            # Extract descendants from each item in the response
            for item in data:
                all_descendants.update(item.get("descendants", []))
        else:
            logging.error(f"Failed to retrieve taxonomy data for batch {batch}. Status code: {response.status_code}")

    # Convert the set of descendants back to a list and return it
    return list(all_descendants)


# TODO check if we should still take into account the cutoff parameter? Maybe this is no longer an issue because
# of the new Unipept API...
def query_unipept_and_filter_taxa(peptides: List[str], taxa_filter: List[int]) -> List[Any]:
    url = UNIPEPT_URL + UNIPEPT_PEPT2FILTERED_ENDPOINT
    filtered_peptides = []

    # Convert to set in order to efficiently perform set intersection later on
    taxa_filter = set(taxa_filter)

    print("Querying unipept and filtering taxa...")

    # Split the peptides into batches of 100
    for i in range(0, len(peptides), 100):
        batch = peptides[i:i+100]

        # Prepare the request payload
        payload = {
            "peptides": batch
        }

        # Perform the HTTP POST request
        response = requests.post(url, json=payload)

        # Check if the request was successful
        if response.status_code == 200:
            data = response.json()

            # Process each peptide result in the response
            for peptide_data in data.get("peptides", []):
                original_taxa = peptide_data.get("taxa", [])

                # Find the intersection of original taxa and taxa_filter
                filtered_taxa = list(set(original_taxa) & taxa_filter)

                # Replace the taxa with the filtered taxa
                peptide_data["taxa"] = filtered_taxa

                # Append the modified peptide data to the result list
                filtered_peptides.append(peptide_data)
        else:
            logging.error(f"Failed to retrieve peptide data for batch {batch}. Status code: {response.status_code}")

    return filtered_peptides


def fetch_unipept_taxon_information(
    pep_score_psm: Dict[str, Tuple[float, int]],
    unipept_peptide_counts_file: str,
    unipept_response_file: str,
    taxonomy_query: str,
    rank: str,
    log_file: str
):
    unipept_peptides = dict()

    for pep in pep_score_psm.keys():
        unipept_peptides[pep] = {
            "score": pep_score_psm[pep][0],
            "psms": pep_score_psm[pep][1],
        }

    with open(unipept_peptide_counts_file, "a") as f_out:
        f_out.write(json.dumps(unipept_peptides))

    # Set up the logger such that we can write errors to the provided file
    logging.basicConfig(filename=log_file, level=logging.INFO)

    taxa_filter = get_descendants_for_taxa([int(item) for item in taxonomy_query.split(",")], rank)
    unipept_responses = query_unipept_and_filter_taxa(list(unipept_peptides.keys()), taxa_filter)

    with open(unipept_response_file, "w") as f_out:
        json.dump(unipept_responses, f_out)

peptides_to_query = MS2RescoreOutParser(args.PoutFile, args.FDR, '')
results = fetch_unipept_taxon_information(peptides_to_query, args.UnipeptPeptides, args.UnipeptResponseFile, args.TaxonomyQuery, args.TaxaRank, args.logfile)


