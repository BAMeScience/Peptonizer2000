<h2> PGM for Metaproteomic and virus strain identification <\h2>
<p> a- use PPI to improve protein inference <br>
b- use pgm to propagate probablities to taxonomic & other identifcation <\p>


<p>structured as a mutilayered graph <br>
One graph that enables various graph analysis steps & pgm <\p>
<br>
_____________________________________________________________________________________________________________
<p> Input: peptide and protein IDs and scores from previous DB searches using viral samples = fileID of my results, other parameters to be defined <br>
Output: Organism/ taxonomic classification with a confidence score, visualization of input/output data, other outputs to be defined <\p>

<br>
<p> 0- load search result data (in the future) <br>
1- query data for proteins found by search engine from viruses.String (in the future) <br>
2- build network for my proteins from viruses.STRING DB  done by String_Load.py     <br>           
3- add peptides & scores generated by cp-dt using Peptide_load.py <br>
4- performs residual message passing using the convolution tree, on separate subgraphs<br>
5- visualization of the factorgraph: save as gml and visualize with graphia </p>


<p> I am using a model similar too https://www.openms.de/comp/epifany/ <\p>

<p> next steps: <br>
- restructure git repo <br>
- write snakemake pipeline to make my life easier<br>
- add funtions to infer target taxa with possible restrictions to certains hosts by the user -> possibly use taxadb package already available?<br>
- add possibility to map peptides only to higher taxonomic layer/ if we have an an unknown virus, output viral families -> maybe just infer them from the virus strains identified?<br>
- find a way to download proteins for metaproteomic samples / generate a DB locally that stores taxon-peptide relationships? (e.g. was also done for bacterial samples for taxit)<br>
- try virus strain identification including non-curated Proteins<br>
- fix the problem of to little valid models identified in my xtandem searches...how is that possible?<br>
- add grid search for best chosen perameters: how to evaluate this? <br>
- strains not int the sample do not seem to get their probability reduced <br> how to (maybe) fix this?




- add GO, EC layers 
<\p>

<p> General mixed TODOs: <br>
- change way of acquiring the peptide-protein graph to allow general user inputs. includes: using scores from DBsearch, downloading proteins for given species/user defined target DB<br>
- add external loop evaluating the parameter settings (do this in extra python script)<br>
- organise my big python script into smaller scripts to be able to execute them serially with something like snakemake<br>
- think about what proteins to fetch from ncbi (only swissprot?)/ whether to construct local databases for each experiment (questions for later, best to talk with others)<br>
- contact Marten danens(?spelling) to see about sars cov 2 strain specific data<br>
- add options/ automated recognition to map peptides only to higher taxa if strain info is not available<br>
- detection & dampening of oscillations, (less urgent) <br>


<p> Ideas sidelined for now/for later:<br>
- pgmpy python library to implement simple version of inference with a noisy Or CPD at each peptide: this doesn't work with pgmpy because there is no way to "add" the proteins together to create the noisyOr and you can't combine the noisy OR and the Markovmodel in pgmpy <br>
- incorporate the connections between co-occuring proteins : use the algorithm proposed in https://advances.sciencemag.org/content/7/17/eabf1211<br>

<p> solved steps log <br>
- restrict update of messages to the parts of the graph that received new info in the last iterations step, until 17.05<br>
- save ad gml and <br>
- generation of taxonomic layer and passing message between taxonomic and petide layer **till 13.06 **. Is there already a tool that assigns proteins to taxa /peptides to taxa? possibility to get taxonomic tree already as  graph (for later visualization) optional: include LCA, set priors to only account for likely present species (?)<br>