---
title: 'SteParKin: a Python code to determine Galactic-space velocity components, evaluate membership to stellar kinematic groups and assign stellar populations'
tags:
  - Python
  - astronomy
authors: 
name: David Montes
    orcid: 0000-0002-7779-238X
    affiliation: 1 
name: Hugo M. Tabernero
    orcid: 0000-0002-8087-4298
    affiliation: 1 
name: Carlos Cifuentes
    orcid: 0000-0003-1715-5087
    affiliation: 2
name: Miriam Cortés-Contreras
    orcid: 0000-0003-3734-9866 
    affiliation: 1 
name: Jose Antonio Caballero 
    orcid: 
    affiliation: 2
affiliations:
name: Departamento de Física de la Tierra y Astrofísica and IPARCOS-UCM (Instituto de Física de Partículas y del Cosmos de la UCM), Facultad de Ciencias Físicas, Universidad Complutense de Madrid, E-28040, Madrid, Spain
index: 1
name: Centro de Astrobiología (CSIC-INTA), ESAC, Camino Bajo del Castillo s/n, 28692 Villanueva de la Cañada, Madrid, Spain
index: 2
date: 23 November 2023
bibliography: paper.bib
---

# Summary

Steparkin is a code that uses the Galactic-space velocity components (U, V, W) of the stars to evaluate their membership to young (< 1 Gyr) kinematic moving groups and associations (5 by default) and to assign their stellar populations as proposed by Bensby et al. 2003, 2005. The code takes as input a file that contains the equatorial coordinates (in J2000 and epoch 2000), proper motions, parallax/distance, and radial velocity (along with their corresponding errors) of a star and retrieves a CSV file and a series of plots (UV and UW planes and the Toomre diagram) as outputs. The CSV file contains the Galactic-space velocity components (U, V, W), their corresponding errors (computed following the approach described in JohnsonSoderblom1987), the associations and/or groups to which the star is candidate, and its most probable stellar population along with the name, the value, and error in radial velocity of the star from the input file.

# Scheme

The code is constituted by three main functions:

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

The standard output file structure is as follows: 

| Column | Description |
| :------:|:-----------|
| NAME |  That provided in the input | 
| RV |  The value provided in the input | 
| ERV |  The value provided in the input | 
| U |  U Galactic-space velocity component in km s-1 (positive towards the Galactic center)  | 
| eU |  error in the U Galactic-space velocity component in km s-1 | 
| V |  V Galactic-space velocity component in km s-1 (positive towards the North Galactic Pole)  | 
| eV |  error in the V Galactic-space velocity component in km s-1 | 
| W |  W Galactic-space velocity component in km s-1 | 
| eW |  error in the W Galactic-space velocity component in km s-1 | 
| lab_YD | indicates if the stars belongs (YD) or not (NYD) to the young disk region defined by Eggen |
| lab_skg | kinematic moving group, association or range of them to which the star is candidate. Default values are Local Association (LA), Hyades Supercluster (HS), IC 2391 (IC), Ursa Major (UMA) and Castor (CAS), their overlaps (LA/IC, LA/CAS, LA/IC/CAS, and IC/CAS) or neither moving group nor association (NNN). |
| lab_kin | stellar population: thin disk (D), thick disk (TD), in between them (D-TD), or halo (H). |

<div><li style="text-align: justify"><b>spk_graphs</b>: This function plots Boettlinger diagrams (the UV and UW planes), with and without zoom in on the Eggen's young disk, and the Toomre diagram. This function also makes use of the parameter file association_parameters.csv. The parameters of the function are: </li></div>

```Python
"""
    Parameters
    ----------

    df_out :  DataFrame
        The DataFrame resulting from spk_groups or, at least, one containing
        the following columns: U, V, and W.
    file_name : str, optional
        A string to be included in the names of the output figures. The default
        value is empty.
    dict_colors_groups : dict or dict-like str, optional
    dict_colors_intergroup_stars : dict or dict-like str, optional
        dict_colors_groups and dict_colors_intergroup_stars are python
        dictionaries or dictionary-like strings containing the name (key) and
        color (value) given by the user to the stars that fall in one
        association/group or between several associations/groups apart from
        those provided by default. To "remove" a default
        association/group/intergroup from the graphics, the corresponding
        dictionary must contain "default name": None (for example, to avoid
        showing LA stars in red, dict_colors_groups must contain "LA": None).
    autocomplete_colors : boolean, optional
        In case the input DataFrame contains stars belonging to any new
        association, group, or intergroup and it has not been defined in the
        corresponding dictionary, the "autocomplete_colors" parameter sets
        whether showing these stars in a randomly generated color (True) or
        displaying them according to their stellar population if given (False).
        Caution: autocomplete_colors only avoid repetition if the colors are
        given by their hexadecimal name.
    independent : boolean, optional
        Different figures for each UV and WV plane? Default is False.
    gf_uw2 : str, optional
        The gf_uw2 parameter sets the aspect ratio of the Toomre diagram. The
        values, circumferences (default) or ellipses, indicate how lines of
        constant total velocity look in the resulting diagram.
"""
```

<div><li style="text-align: justify"><b>spk_wrapper</b>: This function serves as a wrapper for the other two. It takes as input the name of a .csv or .txt file (str) and return a .csv file and the plots (saved as pdf files). </li></div>

<div style="text-align: justify">spk_groups and spk_graphics can be used independently of each other and spk_wrapper in the Python interpretator, allowing the users a more flexible use of the individual functions. Thus, the users can skip the restriction about the input file extension, do only the assignments without doing the plots, or vice versa. </div>


## Running SteParKin

Using the command line:

```
$ python spk_wrapper.py [-h] -i I [-e E] [-fn FN] [-dcg DCG] [-dci DCI][-ac AC] [-ind IND] [-g G]

optional arguments:
  -h, --help  show this help message and exit
  -i I        the input file containing the data (.txt or .csv)
  -e E        set how the code deals with zeros in the errors: drop (default)
              or zeros
  -fn FN      set a string that all the output files must contain
  -dcg DCG    set colors for new groups or change those given by default (str-
              like python dictionary)
  -dci DCI    set colors for new intergroup stars or change those given by
              default (str-like python dictionary)
  -ac AC      should the stars candidates to any association, group, or
              intergroup not defined in the dictionaries be shown in a random
              color? Default is False
  -ind IND    Independent graphs for each plane? Default is False
  -g G        set the aspect ratio of the Toomre diagram: circumferences
              (default) or ellipses
```

Using the Python interpreter:

* Using the wrapper:

```python
>>> import spk_wrapper as spkw
>>> spkw.spk_wrapper(input_file, errors="drop", file_name="", dict_colors_groups=None,
				     dict_colors_intergroup_stars=None, autocomplete_colors=False,
				     independent=False, gf_uw2="circumferences")
```

* Running spk_groups and spk_graphs independently (note that the input is now a DataFrame in both cases):

```python
>>> import spk_groups as spkg
>>> spkg.spk_groups(df, file_name="")
```

```python
>>> import spk_graphs as spkgra
>>> spkgra.spk_groups(input_file, file_name="", dict_colors_groups=None,
				      dict_colors_intergroup_stars=None, autocomplete_colors=False,
				      independent=False, gf_uw2="circumferences")
```


# Acknowledgements

We acknowledge contributions from ...
We acknowledge financial support from the Agencia Estatal de Investigación of the Spanish Ministerio de Ciencia e Innovación through projects PID2019-109522GB-C5[1:4]/AEI/10.13039/501100011033



# References

* [Bensby et al. 2003] (http://ui.adsabs.harvard.edu/abs/2003A%26A...410..527B/abstract)
* [Bensby et al. 2005] (http://ui.adsabs.harvard.edu/abs/2005A%26A...433..185B/abstract)
* [Cortés-Contreras et al. 2020] (https://ui.adsabs.harvard.edu/abs/2020sea..confE.131C/abstract)
* [Cortés-Contreras et al. 2023], A&A, in prep.
* [Crundall et al. 2019] (https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.3625C/abstract)
* [Eggen (1984e)] (https://ui.adsabs.harvard.edu/abs/1984AJ.....89.1350E/abstract)
* [Eggen (1989c)] (https://ui.adsabs.harvard.edu/abs/1989PASP..101..366E/abstract)
* [Gagné et al. 2018] (https://ui.adsabs.harvard.edu/abs/2018ApJ...856...23G/abstract)
* [Johnson & Soderblom (1987)] (https://ui.adsabs.harvard.edu/abs/1987AJ.....93..864J/abstract)
* [Montes et al. 2001] (http://ui.adsabs.harvard.edu/abs/2001MNRAS.328...45M/abstract)
* [Ridel 2016] (https://ui.adsabs.harvard.edu/abs/2016ascl.soft01011R/abstract)


