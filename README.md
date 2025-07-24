### Disclaimer ⚠️
This tool is being developed primarily for educational purposes—to explore how to build tools like this and to create a pharma-oriented coding portfolio. It stands on the project developed at Le Wagon Data Science bootcamp by Alice Lemery, Jean-Baptiste Fallacher, Benjamin Galet PhD and Laura Ossorio Carballo PhD. While we strive for scientific accuracy, any conclusions drawn from its results should not be considered definitive or used as the sole basis for decision-making.

That said, if you run or collaborate with a wet lab specializing in immuno-oncology and are interested in validating this tool through experiments, we’d love to hear from you! Feel free to reach out—we're always open to meaningful collaborations.



# 🛡️ ImmunoReady
Deep Learning for Immune Readiness & Personalized Cancer Vaccine Design


## 🧬 What is ImmunoReady?
ImmunoReady is a deep learning tool trained with a comprehensive immunological dataset to predict whether a MHC class I peptide can trigger an immune response.
It is being developed with the intention of detecting single amino acid changes in immunopeptides — a critical capability for assessing the potential effectiveness of personalized vaccines based on point mutations.

🧠 AI-Powered
🔬 Mutation-Sensitive
💉 Cancer Vaccine-Oriented

## 🎯 Key Features
🧠 Deep Learning Model trained on curated immunological data

🎯 Point Mutation Sensitivity to detect single-residue impact on immunogenicity

🧪 MHC Class I Focus for CD8+ T-cell epitope prediction

🧬 Personalized Vaccine Design Support

📊 Immunogenicity Scoring and interpretability tools

⚡ Fast Inference with minimal input

## Git commit labels

Labels to include at the beginning of commit messages for improved traceability.

- feat:	New feature
- fix:	Bug fix
- dev: development in progress
- docs:	Documentation only
- test:	Adding or updating tests

# Training raw data specifications:

## Peptides from HLA ligand atlas
This data base will be use to construct the category of peptides  would not trigger an immune response as their presence in normal tissues would indicate immunity tolerance.

The list of peptides (MHC Class I and II) was dowloaded the 29th of June 2025 from: https://hla-ligand-atlas.org/data

## Peptides from IEDB
Peptides that are predicted to bind
The API request is WIP, still not fucntional at the moment.

# Feature engineering
- MHC I and II shared epitope: sequence of the MHC class I epitope contained in MHC class II presented epitope.


https://github.com/IEDB/IQ-API-use-cases


## Acnowledgement of published work on the matter:
Li et al., 2021 – DeepImmuno: Deep learning-empowered prediction and generation of immunogenic peptides for T-cell immunity
👉 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7781330/

Wu et al., 2019 – DeepHLApan: A deep learning approach for neoantigen prediction considering both HLA-peptide binding and immunogenicity
👉 https://pubmed.ncbi.nlm.nih.gov/31736974/

Sidhom et al., 2018 – Deep learning for class I and class II HLA binding prediction (AI-MHC)
👉 https://pubmed.ncbi.nlm.nih.gov/33398286/

Diao et al., 2023 – Seq2Neo-CNN: Predicting immunogenicity of tumor neoantigens using convolutional neural networks
👉 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11581883/

Wang et al., 2023 – INeo-Epp: An ensemble learning framework for predicting immunogenic epitopes
👉 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11581883/

Li et al., 2024 – ProVaccine: Predicting immunogenic antigens using dual-attention neural networks
👉 https://arxiv.org/abs/2410.02647

Yu et al., 2024 – UnifyImmun: Cross-attention deep learning for joint prediction of peptide-HLA and peptide-TCR interactions
👉 https://arxiv.org/abs/2405.06653

Ma et al., 2025 – Triad sequence fusion model for predicting peptide–MHC–TCR binding and immunogenicity
👉 https://arxiv.org/abs/2501.01768

Goffinet et al., 2023 – MATE-Pred: A multimodal attention-based predictor of TCR–epitope binding affinity
👉 https://arxiv.org/abs/2401.08619
