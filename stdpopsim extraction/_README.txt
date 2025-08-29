This "stdpopsim extraction" folder contains recombination maps, genetic maps, etc. that were extracted from stdpopsim, together with details (below) regarding how this extraction was done.  Probably I could've just mined the data out of the stdpopsim repository, but that seems maybe error-prone and dependent on the structure of their repository; this way seems more reliable.

This folder contains:

- this file
- generate_data.txt : a Bash script that generates SLiM scripts using stdpopsim
- extract_data.txt : a Bash script that performs the data extraction once generate_data.txt has completed
- raw : a directory containing raw output from stdpopsim in the form of SLiM scripts and notes
- extracted : a directory containing the extracted recombination rate maps, genomic maps, etc.
- out_event.txt : a SLiM script snippet (not a full script) used by extract_data.txt to do its work

----------------------------------------------------------

The data generation provided here was done with stdpopsim version 0.3.0.  NOTE that in stdpopsim 0.3.0 there was a bug in DFE PosNeg_R24, such that the beneficial and deleterious categories were swapped.  I have compensated for that error downstream, specifying the correct DFE in the SLiM scripts provided in this publication.  The stdpopsim bug has no effect on the data that we automatically extract from the scripts, since we do not extract DFE information.

Choices made:

- For the recombination map I chose HapMapII_GRCh38 somewhat arbitrarily.  Maybe DeCodeSexAveraged_GRCh38 would be a good choice too?

- For the DFE I've chosen PosNeg_R24 as one that sounds quite up to date, without getting into the complexity of different dominance coefficients for different selection coefficient categories.  It also includes beneficials, which is nice for demo purposes.  I use the overall mutation rate of 2.0e-8 from this DFE as well, NOT the rate from the demographic model, since the DFE was fitted with 2.0e-8.

- For the coding vs. non-coding regions I've used ensembl_havana_104_CDS; could use ensembl_havana_104_exons instead, but showing coding vs. non-coding seems natural.  Anyhow, this is a somewhat arbitrary choice.

- For the demography I specified OutOfAfrica_2T12, but that aspect of the SLiM code is not used by the data extraction process anyway, so it doesn't matter at all.  We use a different demographic model in our SLiM models.


DATA GENERATION USING STDPOPSIM:

cd ".../stdpopsim extraction"
bash generate_data.txt


DATA EXTRACTION FROM THE PRODUCED SLiM SCRIPTS:

cd ".../stdpopsim extraction"
bash extract_data.txt

