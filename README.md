# Observation Mutual Hazard Networks
This is the repository for the manuscript "Overcoming Observation Bias for Cancer Progression Modeling". It contains the data and resulting models. 
The script ReproduceResults.py learns the models from the data. It requires the official MHN-package which can be installed with 
```bash
pip install mhn
```
The package is currently compatible with Python versions 3.8 to 3.11 but not yet with the most recent version 3.12.
Readers interested in implementation details can access the code at the (package repository)[https://github.com/spang-lab/LearnMHN/]. The learning algorithm for the new observation MHN models exploits their likelihood-equivalence to certain classical MHNs (section 2.3 in the manuscript) and can thus reuse much of the existing numerical machinery.
