# Primer_Design_SCARECROW
Examination project for the Linux and Python class for biologists

# Description
This project offers a web based seat reservation system for a plane with a database.
The project is part of the course "Programmieren für Data Scientists: Python" and is mandatory to pass this course.

---

# Installation
Clone this GitHub repository into your local compiler (PyCharm).

## Libraries
This project uses Bio. Make sure this library are correctly installed in PyCharm in order to run the script.

## Run
*To run the code via the terminal, change the current directory to the one where the script is saved and run it by entering ```python primer_design.py```. You will then be asked to enter information via the terminal.
*To run the code via the PyCharm Console, click the "Run" button and emter the information via the console.

---

# Application
*First, you will be asked to enter an email address which is required to access NCBI.

```
Please enter an email address: 
```

*Second, you will be asked to enter one or multiple NCBI Reference Sequence for the gene(s) the primers should be designed for. If not already downloaded, this program will download the FASTA files for the given NCBI references and save it in the current directory automatically.

```
Please enter one or multiple (comma-seperated) NCBI Reference Sequence (ex. NM_115282.4, NM_103925.6): 
```

*Third, you will be asked for the design criteria for the primer design. The default values are given in parenthesis.

```
Enter min. primer length (ex. 20): 
Enter max. primer length (ex. 23): 
Enter min. annealing the temperature in °C (ex. 55.0): 
Enter max. annealing the temperature in °C (ex. 62.0): 
Enter max. difference in annealing temperatures between forward and reverse primers in °C (ex. 4.0): 
```

## Results
The respective result for each given NCBI reference consits of:
*Information on the gene
*Coding Region
*Comment on whether the designed primers fulfill the predefined criteria
*Forward and reverse primer in 5'-3' orientation, each with their melting temperature and GC-content
*The difference in annealing temperatures between forward and reverse primer

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


---
	
# Ideas to extend the project
*include the GC-content in the predefined criteria
*allow aimed point mutations in primer sequences via the console to increase or decrease the annealing temperature

---


Göttingen, 23.04.2023

Fiedler, Carlotta 
