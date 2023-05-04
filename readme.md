<div id="top"></div>


<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href=https://git.bam.de/tholstei/pepgm/>
    <img src="images/peptonizer.jpg" alt="Logo"  height="150">
  </a>

<h3 align="center">The Peptonizer 2000</h3>

  <p align="center">
    Integrating PepGM and Unipept for probability-based taxonomic inference of metaproteomic samples
    <br />
  </p>
</div>


<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
      </ul>
    </li>
    <li><a href="#input">Input</a></li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
        <li><a href="#preparation">Preparation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

Introducing the Peptonizer2000 - a tool that combines the capabilities of Unipept and PepGM to analyze
metaproteomic mass spectrometry-based samples. Originally designed for taxonomic inference of viral
mass spectrometry-based samples, we've extended PepGM's functionality to analyze metaproteomic samples by
retrieving taxonomic information from the Unipept database.

PepGM is a probabilistic graphical model developed by the eScience group at BAM (Federal Institute for Materials
Research and Testing) that uses belief propagation to infer the taxonomic origin of peptides and taxa in viral samples.
You can learn more about PepGM on our eScience group at BAM (Federal Institute for Materials Research and Testing).
Please refer to our [GitHub](https://github.com/BAMeScience/PepGM) page.

Unipept, on the other hand, is a web-based metaproteomics analysis tool that provides taxonomic information for
identified peptides. To make it work seamlessly with PepGM, we've extended Unipept with new functionalities that
restrict the taxa queried and provide all potential taxonomic origins of the peptides queried. Check out more
information about Unipept [here](https://unipept.ugent.be/).

With the Peptonizer2000, you can look forward to a comprehensive and streamlined workflow that simplifies
the process of identifying peptides and their taxonomic origins in metaproteomic samples.

The Peptonizer2000 workflow is comprised of the following steps:

1. Start by conducting a database search using X!Tandem and [MS2Rescore](https://github.com/compomics/ms2rescore).
   The reference database has to be provided by the user.
2. Query all identified peptides in the Unipept API,
   and restrict the taxonomic range queried based on any prior knowledge of the sample.
3. Assemble the peptide-taxon associations provided by Unipept into a bipartite graph,
   where peptides and taxa are represented by different nodes, and an edge is drawn between a peptide and a taxon
   if the peptide is part of the taxon's proteome.
4. Transform the bipartite graph into a factor graph using convolution trees and conditional probability table
   factors (CPD).
5. Run the belief propagation algorithm multiple times with different sets of CPD parameters until convergence,
   to obtain posterior probabilities of candidate taxa.
6. Use an empirically deduced metric to determine the ideal graph parameter set.
7. Output the top 15 scoring taxa as a results barchart. The results are also available as comma-separated files
   for further downstream analysis or visualizations.


<div align="center">
    <img src="images/workflow.png" alt="worklfow scheme" width="500">
</div>

<br>



<p align="right">(<a href="#top">back to top</a>)</p>

<!-- INPUT -->

## Input

* Your raw spectrum file in mgf format
* A reference database in fasta format <br>

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

Install PepGM. Find installation instructions [here](https://github.com/BAMeScience/PepGM).

### Configuration file 

The Peptonizer2000 relies on a configuration file in `yaml` format to set up the workflow.
An example configuration file is provided in `config/config.yaml`. <br>
Do not change the config file location.

<details>
  <summary>Details on the configuration parameters </summary> <br>
    <details> <summary>Run panel <br> </summary> 
    Set up the workflow of your Peptonizer run by providing parameters that fill wildcards to locate input files
    such as raw spectra or reference database files. Thus, use file basenames i.e., without file 
    suffix, that your files already have or rename them accordingly. <br><br>
    Run: Name of your run that is used to create a subfolder in the results directory. <br>
    Sample: Name of your sample that is used to create a subfolder in the run directory. <br> 
    Reference: Name of reference database (e.g. human). <br>

    <details> <summary>Input panel <br> </summary> 
    Specify input file and directory paths. <br><br>
    Sample spectra: Path to raw spectra file. <br>   
    Parameter: Path to SearchGUI parameter file. <br>
    Sample data: Path to directory that contains sample raw spectra files. <br>
    Database: Path to directory that contains the reference database. <br>
    Resources: Relative path to resources folder <br>
    Results: Relative path to results folder <br>
    TaxID mapping: Relative path to folder that contains mapped taxIDs. <br> <br>
    </details>
    <details> <summary>X!Tandem panel <br> </summary> 
    Fill in the parameters necessary for the X!Tandem search. Please refer to the [X!Tandem documentation](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0091507) regarding the parameters. <br><br>
    </details>
    <details> <summary>MS2Rescore panel <br> </summary> 
    Fill in the parameters necessary for the MS2Rescore rescoring. Please refer to the MS2Rescore documentation [X!Tandem documentation](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0091507) regarding the parameters. <br><br>
    </details>
    <details> <summary>PepGM panel <br></summary>
    Grid search: Choose increments for alpha, beta and prior that are to be included in the grid search to tune
    graphical model parameters. Do not put a comma between values. For metaproteomics analysis, you can significantly speed up the Peptonizer by reducing the amount of parameters searched. For samples where you expect good peptide coverage per species, we recommmend alpha between 0.8 and 0.9, beta between 0.2 and 0.4 and gamma of 0.5. <br>
    Results plotting: Number of taxa in the final species identification barplot. <br><br> 
    </details>

</details>


### Through the command line

The Peptonizer2000 needs to be run from the command line. To run the snakemake workflow, 
you need to be in your PepGM repository and have the Snakemake conda environment activated. 
Run the following command 
```sh
  snakemake --use-conda --conda-frontend conda --cores <n_cores> 
  ```
Where `n_cores` is the number of cores you want snakemake to use. 
<p align="right">(<a href="#top">back to top</a>)</p>

### Output files

All Peptonizer2000 output files are saved into the results folder and include the following: <br>

Main results: <br>
- PepGM_Results.csv: Table with values ID, score, type (contains all taxids under 'ID' and all probabilities under 'score' tosterior probabilities of n (default: 15) highest scoring taxa <br>
 <br>

Additional (intermediate): <br>
- Intermediate results folder sorted by their prior value for all possible grid search parameter combinations
- mapped_taxids_weights.csv: csv file of all taxids that had at least one protein map to them and their weight 
- PepGM_graph.graphml: graphml file of the graphical model (without convolution tree factors). Useful to visualize the graph structure and peptide-taxon connections <br>
- paramcheck.png: barplot of the metric used to determine the graphical model parameters for n (default: 15) best performing parameter combinations <br>
- log files for bug fixing

## Toy example
Will be provided.



<p align="right">(<a href="#top">back to top</a>)</p>

<!-- ROADMAP -->
## Roadmap

- [ ] Extensive testing
- [ ] Selection of taxonomic resolution for the results


See the [open issues](https://github.com/BAMeScience/repo_name/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Tanja Holstein - [@HolsteinTanja](https://twitter.com/HolsteinTanja) - tanja.holstein@bam.de <br>
Franziska Kistner - [LinkedIn](https://www.linkedin.com/in/franziska-kistner-58a57b18b) - franziska.kistner@bam.de
Pieter Verschaffelt - pieter.verschaffelt@ugent.be

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/BAMeScience/repo_name.svg?style=for-the-badge
[contributors-url]: https://github.com/BAMeScience/repo_name/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/BAMeScience/repo_name.svg?style=for-the-badge
[forks-url]: https://github.com/BAMeScience/repo_name/network/members
[stars-shield]: https://img.shields.io/github/stars/BAMeScience/repo_name.svg?style=for-the-badge
[stars-url]: https://github.com/BAMeScience/repo_name/stargazers
[issues-shield]: https://img.shields.io/github/issues/BAMeScience/repo_name.svg?style=for-the-badge
[issues-url]: https://github.com/BAMeScience/repo_name/issues
[license-shield]: https://img.shields.io/github/license/BAMeScience/repo_name.svg?style=for-the-badge
[license-url]: https://github.com/BAMeScience/repo_name/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 