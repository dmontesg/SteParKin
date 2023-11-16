---
title: 'SteParKin: a Python code to determine Galactic-space velocity components, evaluate membership to stellar kinematic groups and assign stellar populations'
tags:
  - Python
  - astronomy
authors: 
  - name: David Montes
    orcid: 0000-0002-7779-238X
    equal-contrib: true
    affiliation: 1 
 - name: Hugo M. Tabernero
    orcid: 0000-0002-8087-4298
    affiliation: 1 
affiliations:
 - name: Departamento de Física de la Tierra y Astrofísica and IPARCOS-UCM (Instituto de Física de Partículas y del Cosmos de la UCM), Facultad de Ciencias Físicas, Universidad Complutense de Madrid, E-28040, Madrid, Spain
   index: 1
date: 16 November 2023
bibliography: paper.bib
---

# Summary

Steparkin is a code that uses the Galactic-space velocity components (U, V, W) of the stars to evaluate their membership to young (< 1 Gyr) kinematic moving groups and associations (5 by default) and to assign their stellar populations as proposed by Bensby et al. 2003, 2005. The code takes as input a file that contains the equatorial coordinates (in J2000 and epoch 2000), proper motions, parallax/distance, and radial velocity (along with their corresponding errors) of a star and retrieves a CSV file and a series of plots (UV and UW planes and the Toomre diagram) as outputs. The CSV file contains the Galactic-space velocity components (U, V, W), their corresponding errors (computed following the approach described in `@JohnsonSoderblom1987`), the associations and/or groups to which the star is candidate, and its most probable stellar population along with the name, the value, and error in radial velocity of the star from the input file.

# Statement of need

bla bla

# Scheme

The code is constituted by 3 main functions:

<div><li style="text-align: justify"><b>spk_groups</b>: This function associates stars to a young moving groups or associations if its Galactic-space velocity components are within the 3D ellipsoid that defines each moving group or association. The ellipsoids are defined in the UVW space and their parameters are recorded in the association_parameters.csv file. This function also determines to which stellar population (thin disk (D), thick disk (TD), in between them (TD-D), or halo (H)) the star belongs following the probabilistic approach described in Bensby et al. 2003, 2005. The values of the parameters that define each population (the observed fraction in the solar neighborhood (X_ns), the characteristic velocity dispersions (sig_u, sig_v, sig_w), and the asymmetric drift (v_asym)) are stored in the param_prob_populations.csv file. Note that CSV files must be in the same directory as the function (spk_groups.py) for this to work. The optional parameter file_name takes a string to be included in the names of the output files. The default value is empty.

The input file (spk_wrapper, see below) or input DataFrame (spk_groups) must have the following columns: </li></div>

| Column | Description |
| :------:|:-----------|
| NAME |  name of the star (numbers are allowed) | 
| RA |  right ascension in degrees | 
| DEC|  declination in degrees | 
| PMRA|  proper motion in right ascension (mu<sub>alpha</sub>cos(delta)) in mas yr-1 | 
| EPMRA |  error of proper motion in right ascension in mas yr-1 | 
| PMDEC|  proper motion in declination in mas yr-1 | 
| EPMDEC|  error of proper motion in declination in mas yr-1 | 
| RV |  radial velocity in km s-1 | 
| ERV |  error of radial velocity in km s-1 | 
| PLX (optional<sup>*</sup>) |  parallax in mas | 
| EPLX (optional<sup>*</sup>) |  error of parallax in mas | 
| d (optional<sup>*</sup>)  |  distance in pc | 
| ed (optional<sup>*</sup>)  |  error of distance in pc | 

<blockquote>
  
<h6> <sup>*</sup> The user use spk_wrapper can choose between giving the parallax or the distance (and the corresponding error) of the stars (always one of them for each set analyzed at a time). If only distances and their error are given, the parallaxes and their errors are calculated from them. If both are provided, distance is ignored. However, if spk_groups is used independently, only PLX and EPLX are allowed. </h6>

</blockquote>


# Acknowledgements

We acknowledge contributions from ...
We acknowledge financial support from the Agencia Estatal de Investigación of the Spanish Ministerio de Ciencia e Innovación through projects PID2019-109522GB-C5[1:4]/AEI/10.13039/501100011033



# References

Johnson & Soderblom 1987
