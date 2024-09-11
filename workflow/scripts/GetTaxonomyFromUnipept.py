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
parser.add_argument('--TaxaRank', type = str, nargs = '+', required = True, help = 'taxonomic rank at which you want the peptonizer results resolved')
parser.add_argument('--logfile', type = str, required = True, help = 'path to logfile')

args = parser.parse_args()


from typing import Dict, Tuple, List, Any


from taxon_manager import TaxonManager

UNIPEPT_URL = "https://api.unipept.ugent.be"
UNIPEPT_PEPT2FILTERED_ENDPOINT = "/mpa/pept2filtered.json"

UNIPEPT_PEPTIDES_BATCH_SIZE = 500


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
                    pep_score_psm[peptide] = {
                    "score": pep_score[peptide],
                    "psms": len(pep_psm[peptide])
                    }   
                
    return pep_score_psm

# TODO check if we should still take into account the cutoff parameter? Maybe this is no longer an issue because
# of the new Unipept API...
def query_unipept_and_filter_taxa(peptides: List[str], taxa_filter: List[int]) -> List[Any]:
    url = UNIPEPT_URL + UNIPEPT_PEPT2FILTERED_ENDPOINT
    filtered_peptides = []

    # Convert to set in order to efficiently perform set intersection later on
    taxa_filter = set(taxa_filter)

    print("Querying unipept and filtering taxa...")

    batches = len(peptides) // UNIPEPT_PEPTIDES_BATCH_SIZE

    # Split the peptides into batches of 100
    for i in range(0, len(peptides), UNIPEPT_PEPTIDES_BATCH_SIZE):
        print(f"Now querying Unipept with batch {(i // UNIPEPT_PEPTIDES_BATCH_SIZE) + 1} out of {batches}.")
        batch = peptides[i:i+UNIPEPT_PEPTIDES_BATCH_SIZE]

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
    peptide_scores: Dict[str, Dict[str, float | int]],
    taxonomy_query: str,
    rank: str,
    log_file: str
) -> List[Any]:
    # Set up the logger such that we can write errors to the provided file
    logging.basicConfig(filename=log_file, level=logging.INFO)

    taxa_filter = TaxonManager.get_descendants_for_taxa([int(item) for item in taxonomy_query.split(",")], rank)
    return query_unipept_and_filter_taxa(list(peptide_scores.keys()), taxa_filter)

    


peptides_to_query = MS2RescoreOutParser(args.PoutFile, args.FDR, '')
with open(args.UnipeptPeptides, "w") as f:
    json.dump(peptides_to_query, f)
results = fetch_unipept_taxon_information(peptides_to_query, args.TaxonomyQuery, args.TaxaRank, args.logfile)
with open(args.UnipeptResponseFile, 'w') as f:
    json.dump(results, f)
