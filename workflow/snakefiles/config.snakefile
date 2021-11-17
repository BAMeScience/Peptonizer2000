configfile: 'config/config.yaml'

DataDirectory = config['DataDir']
DatabaseDirectory = config ['DatabaseDir']
ResultsDir = config['ResultsDir']

PeptideShakerDir = config['PeptideShakerDir']
SearchGUIDir = config['SearchGUIDir']

searchengines = config['searchengines']
peptideFDR = config['peptideFDR']
proteinFDR = config['proteinFDR']
psmFDR = config['psmFDR']

SpectraFileType = config['SpectraFileType']
SampleName = config['SampleName']
HostName = config['HostName']
ReferenceDBName = config['ReferenceDBName']

TaxaInPlot = config['TaxaInPlot']
TaxaInProteinCount = config['TaxaInProteinCount']

TargetTaxa = config['TargetTaxa']
firstTarget = config['FirstTarget']

AlphaRange = config['Alpha']
BetaRange = config['Beta']
Prior = config['prior']