# Primer_Design_SCARECROW


## Description
This project is the examination project for the Linux and Python class for biologists.

### Files
* ```description.txt```	contains the project description and the criteria that the programm respectively the primers must fulfill.
* ```primer_design.py``` contains the script.

There is an example provided to see how the code works:
* ```SCR_SCL_genes_ncbi_ex.txt``` contains a list of SCARECROW and SCARECROW-Like genes in the plant Arabidopsis thaliana and their NCBI-IDs that the program asks to enter via the console as theiy are required to access the correct FASTA files.
* ```SCR_SCL_primer_ex.txt``` contains the saved output of the program when all listed NCBI-IDs of the SCARECROW and SCARECROW-Like genes were entered.


---

## Installation
Clone this GitHub repository into your local compiler (PyCharm).

### Libraries
This project uses Bio to download FASTA files from NCBI. Make sure this library are installed correctly in PyCharm to run the script.

### Run
* To run the code via the terminal, change the current directory to the one where the script is saved and run it by entering ```python primer_design.py```. You will then be asked to enter information via the terminal.
* To run the code via the PyCharm Console, click the "Run" button and enter the information via the console.

---

## Application
* First, you will be asked to enter an email address which is required to access NCBI.

```
Please enter an email address: 
```

* Second, you will be asked to enter one or multiple NCBI Reference Sequence for the gene(s) the primers should be designed for. If not already downloaded, this program will download the FASTA files for the given NCBI references and save it in the current directory automatically (... this might take a while).

```
Please enter one or multiple (comma-seperated) NCBI Reference Sequence (ex. NM_115282.4, NM_103925.6): 
```

* Third, you will be asked for the design criteria for the primer design. The default values are given in parenthesis.

```
Enter min. primer length (ex. 20): 
Enter max. primer length (ex. 23): 
Enter min. annealing the temperature in °C (ex. 55.0): 
Enter max. annealing the temperature in °C (ex. 62.0): 
Enter max. difference in annealing temperatures between forward and reverse primers in °C (ex. 4.0): 
```

### Results
The respective result for each given NCBI reference consits of:
* Information on the gene
* Coding Region
* Comment on whether the designed primers fulfill the predefined criteria
* Forward and reverse primer in 5'-3' orientation, each with their melting temperature and GC-content
* The difference in annealing temperatures between forward and reverse primer

Example:
```
>NM_115282.4 Arabidopsis thaliana GRAS family transcription factor (SCR), mRNA
Coding region: 5'-ATG GAT CGG CAA ACG GAA GAC GTC AAA CAC ACA ACG ACG AAC ATT TTC CGA TCA CCC ACC TAA-3'

This primer pair for >NM_115282.4 checks all predefined criteria:
FORWARD primer: 5'-GGTCAAGAAGAAACGAAATG-3' 	(56.00°C, GC-content: 40.0 %)
REVERSE primer: 5'-AAATGGGAAGAGATTAGGTG-3' 	(56.00°C, GC-content: 40.0 %)
Difference in annealing Temperature: 	0°C
```
### Unmatched Design Criteria
If the predefined criteria cannot be matched, the output will contain a comment on what criteria is not fulfilled and the best found primer pair is printed.
For example:
```
The annealing temperatures of both primers for >NM_115282.4 are not within the predefined range.
``` 
or
```
The annealing temperatures of the FORWARD primer for >NM_115282.4 is not within the predefined range.
``` 
or
```
The annealing temperatures of the REVERSE primer for >NM_115282.4 is not within the predefined range.
``` 
or
```
The difference in annealing temperatures of the primers for >NM_115282.4 exceeds the predefined maximum.
``` 
Followed by:
```
This is the best match for the primer pair:
```
---
## Export Data
In the end the user is asked whether they want to save the generated output. If not (```n```) the program ends. If yes (```y```) the output is written and saved in a newly generated .txt file, which name contains the current date and time, in the current directory.

---
	
## Ideas to extend the project
* include the GC-content in the predefined criteria
* allow aimed point mutations in primer sequences via the console to increase or decrease the annealing temperature

---


Göttingen, 23.04.2023

Fiedler, Carlotta 
