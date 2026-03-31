# HyperbolIDRInteractome

A modular analysis pipeline for studying **intrinsically disordered regions (IDRs)** in the human proteome through the lens of **hyperbolic network geometry**. The pipeline embeds the STRING human protein–protein interaction (PPI) network in two-dimensional hyperbolic space (H²) and integrates IDR content, liquid–liquid phase separation (LLPS) propensity, network topology metrics, and community structure to test whether the position of a protein in hyperbolic space is a predictor of its disorder properties and interaction behaviour.

---

## Citation

If you use this pipeline, please cite the associated manuscript:

> [Authors]. [Title]. *In preparation* (2026).

Please also cite the underlying tools and databases:

- **AIUPred** – Erdős, G. & Dosztányi, Z. AIUPred: combining energy estimation with deep learning for the enhanced prediction of protein disorder. *Nucleic Acids Res.* **52**, W176–W181 (2024).
- **AlphaFold2** – Jumper, J. et al. Highly accurate protein structure prediction with AlphaFold. *Nature* **596**, 583–589 (2021). / Fleming, J. et al. AlphaFold protein structure database and 3D-Beacons: new data and capabilities. *J. Mol. Biol.* **437**, 168967 (2025).
- **FuzDrop** – Miskei, M.; Horvath, A.; Vendruscolo, M. & Fuxreiter, M. Sequence-based prediction of fuzzy protein interactions. *J. Mol. Biol.* **432**, 2289–2303 (2020). / Horvath, A.; Miskei, M.; Ambrus, V.; Vendruscolo, M. & Fuxreiter, M. Sequence-based prediction of protein binding mode landscapes. *PLoS Comput. Biol.* **16**, e1007864 (2020). / Vendruscolo, M. & Fuxreiter, M. FuzDrop: sequence-based prediction of the propensity of proteins for liquid–liquid phase separation and aggregation. *Nature Protocols*, 1–27 (2026).
- **STRING** – von Mering, C. et al. STRING: a database of predicted functional associations between proteins. *Nucleic Acids Res.* **31**, 258–261 (2003). / Szklarczyk, D. et al. The STRING database in 2025: protein networks with directionality of regulation. *Nucleic Acids Res.* **53**, D730–D737 (2025).
- **Gene Ontology** – Ashburner, M. et al. Gene ontology: tool for the unification of biology. *Nat. Genet.* **25**, 25–29 (2000). / The Gene Ontology Consortium. The Gene Ontology knowledgebase in 2023. *Genetics* **224**, iyad031 (2023).
- **Reactome** – Milacic, M. et al. The Reactome pathway knowledgebase 2024. *Nucleic Acids Res.* **52**, D672–D678 (2024).
- **UniProt** – The UniProt Consortium. UniProt: the universal protein knowledgebase. *Nucleic Acids Res.* **46**, 2699 (2018).
- **InterPro** – Apweiler, R. et al. The InterPro database, an integrated documentation resource for protein families, domains and functional sites. *Nucleic Acids Res.* **29**, 37–40 (2001). / Blum, M. et al. InterPro: the protein sequence classification resource in 2025. *Nucleic Acids Res.* **53**, D444–D456 (2025).
- **ELM** – Dinkel, H. et al. ELM—the database of eukaryotic linear motifs. *Nucleic Acids Res.* **40**, D242–D251 (2012). / Kumar, M. et al. ELM—the eukaryotic linear motif resource—2024 update. *Nucleic Acids Res.* **52**, D442–D455 (2024).
- **PROSITE** – Sigrist, C. J. et al. The PROSITE database for protein families, domains and sites. *Nucleic Acids Res.* **54**, D451–D458 (2026).
- **CD-CODE** – Kuznetsova, K. et al. CD-CODE 2.0: an enhanced condensate knowledgebase integrating pathobiology, condensate modulating drugs and host–pathogen interactions. *Nucleic Acids Res.* **54**, D375–D382 (2026).
- **Mercator** – García-Pérez, G.; Allard, A.; Serrano, M. Á. & Boguñá, M. Mercator: uncovering faithful hyperbolic embeddings of complex networks. *New J. Phys.* **21**, 123033 (2019).
- **iGraph** – Csardi, G. & Nepusz, T. The igraph software package for complex network research. *InterJournal Complex Systems* 1695 (2006).
- **ReactomePA** – Yu, G. & He, Q.-Y. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. *Mol. BioSyst.* **12**, 477–479 (2016).

---

## Contact

Frank Hause  
Center for Structural Mass Spectrometry  
Martin Luther University Halle-Wittenberg  
frank.hause@pharmazie.uni-halle.de
