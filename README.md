# ZZ-exchange-analysis
To fit peak intensities of ZZ-exchange diagonal and cross peaks to solutions of Bloch-McConnell equations for a simple two state exchange process
## Usage and Description
Take the peak intensity data from PINT lineshape analysis software, for both the diagonal peaks and cross peaks. The code takes a list of residues in the format resno-res (eg. 2M, 54G etc.). For example, for residue metheonine at position 2 in the proetin sequence, the labeling of the diagonal peaks should be 2Ma, 2Mb for major and minor conformation respectively. The labeling of the cross peaks should be 2Mab, 2Mba with the conformation correlation label (ab or ba) is done mathcing the proton dimensions of the major and minor states.

The code fits the two diagonal peaks and two cross peaks simulataneously and give the fit plots as svg and pdf. The fit parameters are also written to a tab separated file.
