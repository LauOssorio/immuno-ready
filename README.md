# Disclaimer ⚠️
This tool is being developed primarily for educational purposes—to explore how to build tools like this and to create a pharma-oriented coding portfolio. It stands on the project developed at Le Wagon Data Science bootcamp by Alice Lemery, Jean-Baptiste Fallacher, Benjamin Galet and Laura Ossorio Carballo. While we strive for scientific accuracy, any conclusions drawn from its results should not be considered definitive or used as the sole basis for decision-making.

That said, if you run or collaborate with a wet lab specializing in immuno-oncology and are interested in validating this tool through experiments, we’d love to hear from you! Feel free to reach out—we're always open to meaningful collaborations.

# Hello, Immuno-World! 🌍
We are currently in the early stages of setting up the ImmunoReady tool. This README is a work in progress and will be updated frequently as development continues.

## Initial Setup
To get started, clone or download this GitHub repository to your local machine:

`git clone https://github.com/LauOssorio/immuno-ready.git`

## Git commit labels

Labels to include at the beginning of commit messages for improved traceability.

- feat:	New feature
- fix:	Bug fix
- docs:	Documentation only	docs
- refactor:	Code changes that don’t change behavior	refactor
- test:	Adding or updating tests	test

# Training raw data specifications:

## Peptides from HLA ligand atlas
This data base will be use to construct the category of peptides  would not trigger an immune response as their presence in normal tissues would indicate immunity tolerance.

The list of peptides (MHC Class I and II) was dowloaded the 29th of June 2025 from: https://hla-ligand-atlas.org/data

## Peptides from IEDB

![Alt text](doc/img/doc/img/IEDB_search_filters_for_training_set2025-06-30_at_19.50.37.png)


The API request is WIP, still not fucntional at the moment.

# Feature engineering
- MHC I and II shared epitope: sequence of the MHC class I epitope contained in MHC class II presented epitope.


https://github.com/IEDB/IQ-API-use-cases
## Acnowledgement of published work on the matter:
TO DO
