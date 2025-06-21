# SimHumanity
_SLiM scripts for full-human-genome evolutionary simulations_

---

This repository contains scripts and data associated with the publication:

> Haller, B.C., Nelson, C.W., and Messer, P.W. (2025). SimHumanity: Using SLiM 5.0 to run full-genome simulations of human evolution. [Journal info TBD]

In particular, it contains:

- `README.md` : this file

- `LICENSE` : specifying the license for this repository

- `model1_foundation.slim` : the first version of the SimHumanity model presented in the paper (Code Sample 1)

- `model1_foundation_FIGURE.slim` : an extended version of `model1_foundation.slim` that adds plotting code that produces Fig. 1

- `model2_demography.slim` : the second version of the SimHumanity model presented in the paper (Code Sample 2)

- `model2_demography_FIGURE.slim` : an extended version of `model2_demography.slim` that adds plotting code that produces Fig. 2

- `model3_treeseq.slim` : the third version of the SimHumanity model presented in the paper (Code Sample 3)

- `model3_treeseq_FIGURE.slim` : an extended version of `model3_treeseq.slim` that sets the seed value used to produce Fig. 3

- `model3_analysis.py` : a Python script that performs recapitation, neutral mutation overlay, and analysis for the trees archive produced by `model3_treeseq_FIGURE.slim`, leading to Fig. 3

- `stdpopsim extraction` : a folder containing genetic structure data for all chromosomes except the mtDNA, together with scripts and materials for extracting that genetic structure data from `stdpopsim`

- `mtDNA info` : a folder containing genetic structure data for the mtDNA chromosome, together with info about how we obtained that data from the NCBI Reference Sequence we used

If you have questions, please contact the corresponding author(s) for the paper.
