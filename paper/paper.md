title: 'SteParKin: a Python code to determinae Galactic-space velocity components and evaluate membership to stellar kinematic groups and associations and assign stellar populations'
tags:
  - Python
  - astronomy
authors:
  - name: David Montes
    orcid: 0000-0002-7779-238X
    equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)

affiliations:
 - name: Departamento de F\'{i}sica de la Tierra y Astrof\'{i}sica 
and IPARCOS-UCM (Instituto de F\'{i}sica de Part\'{i}culas y del Cosmos de la UCM), 
Facultad de Ciencias F\'{i}sicas, Universidad Complutense de Madrid, E-28040, Madrid, Spain
   index: 1

# Summary

Steparkin is a code that uses the Galactic-space velocity components (U, V, W) of the stars to evaluate their membership to young (< 1 Gyr) kinematic moving groups and associations (5 by default) and to assign their stellar populations as proposed by Bensby et al. 2003, 2005. The code takes as input a file that contains the equatorial coordinates (in J2000 and epoch 2000), proper motions, parallax/distance, and radial velocity (along with their corresponding errors) of a star and retrieves a CSV file and a series of plots (UV and UW planes and the Toomre diagram) as outputs. The CSV file contains the Galactic-space velocity components (U, V, W), their corresponding errors (computed following the approach described in Johnson & Soderblom (1987)), the associations and/or groups to which the star is candidate, and its most probable stellar population along with the name, the value, and error in radial velocity of the star from the input file.

# Statement of need


# Acknowledgements

We acknowledge contributions from


# References


