## Explorations on class-specificity, early mitigation, and ambiguous peptides on meaningful batch effect correction in mass-spectrometry-based proteomics --- A commentary on proBatch

Codes used to analyse and generate figures for the manuscript. 

### /codes/ folder

Python and R codes used for the analyses.

GPCA_R.ipynb - gPCA and class-specific adaptation of gPCA to extract delta

pca_kw.py - Principal component analysis association test.

### /matrices/ folder

Matrices used to generate plots

gr.csv: Uncorrected protein data

grp_jelena.csv: Williams et al (2021) median centering

rcl.csv: Class-specific protein level

rgl.csv: Global protein level

ecl.csv: Class-specific peptide level

egl.csv: Global peptide level

gcl.csv: Class-specific ambiguous peptide guided

ggl.csv: Global ambiguous peptide guided


## Acknowledgements

 - [proBatch resource by ‪Čuklina et al (2021)](https://github.com/symbioticMe/proBatch)
 

Williams, E.G., Pfister, N., Roy, S., Statzer, C., Haverty, J., Ingels, J., Bohl, C., Hasan, M., Čuklina, J., Bühlmann, P., Zamboni, N., Lu, L., Ewald, C.Y., Williams, R.W., Aebersold, R., 2021. Multiomic profiling of the liver across diets and age in a diverse mouse population. Cell Systems S2405471221003446. https://doi.org/10.1016/j.cels.2021.09.005 *Source of data*

Čuklina, J., Lee, C.H., Williams, E.G., Sajic, T., Collins, B.C., Rodríguez Martínez, M., Sharma, V.S., Wendt, F., Goetze, S., Keele, G.R., Wollscheid, B., Aebersold, R., Pedrioli, P.G.A., 2021. Diagnostics and correction of batch effects in large‐scale proteomic studies: a tutorial. Mol Syst Biol 17. https://doi.org/10.15252/msb.202110240 *proBatch resource*
