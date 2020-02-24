**JASPrecipitationExtremes**

The data and code in this repository can be used to reproduce "Convective dynamics and the response of precipitation extremes to warming in radiative-convective equilibrium" published in *Journal of the Atmospheric Sciences* by myself, Tim Cronin, and Tom Beucler.

Please note that the code in this repository comes without any warranty, and without even the implied warranty of metchantibility or fitness for a particular purpose.

This respository contains several folders:

- `data` contains postprocessed simulation data and an optical rain gauge precipitation timeseries from the Nauru DOE ARM facility. This repository does *not* contain raw simulation output, but it's available on request (email me at thabbott at mit dot edu).
- `figures` contains matlab scripts for recreating all non-appendix figures.
- `figures_rcemip*` contain matlab scripts for recreating figures from the appendices that use data from small-domain simulations.
- `manuscript` contains TeX source files for compiling multi-panel figures and the manuscript itself.

The repository also contains three bash scripts to automatically regenerate all figures and recompile the paper. They require `matlab`, `pdflatex`, and `bibtex` to be available on the search path. To use the scripts, run
```bash
$ bash reproduce_figures.sh	# Generate all figures---may take some time
$ bash reproduce_tables.sh	# Calculate the values in Table A1
$ bash recompile_manuscript.sh	# Compile multi-panel figures and the manuscript
```
If they succeed, they should produce the [arXiv preprint of this paper](https://arxiv.org/abs/1909.01941) as a PDF at `figures/JASPrecipitationExtremes.pdf`.
