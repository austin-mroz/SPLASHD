
<h1 align="center">
    <br>
    <img src="./figures/splashd_logo.png" alt="SPLASHD" width="200">
    <br>
    SPLASHD
    <br>
</h1>

<h4 align="center">Porous liquid size-exclusivity prediction workflow</h4>

<p align="center">
    <a href="#key-features">key features</a> •
    <a href="#installation">installation</a> •
    <a href="#examples">examples</a> •
    <a href="#citation">citation</a> •
    <a href="#acknowledgements">acknowledgements</a> •
    <a href="#license">license</a>
</p>

<br>

## key features
<table>
<tr>
<td>

`SPLASHD` is a porous liquid (PL) size-exclusivity prediction workflow to 
accelerate the discovery of PLs.

</td>
<td>

<img src="./figures/vis.png" alt="vis" width="1000">

</td>
</tr>

<tr>
<td colspan="2">

`SPLASHD` is comprised of 3 main steps

<img src="./figures/workflow.png" alt="workflow" width="2000">

* **system setup** -- a library of candidate solvents and a library of 
candidate porous motifs are generated and optimized in separate workflows 
* **size analysis** -- conformers of each candidate compound are extracted, 
and their size subsequently analyzed. Porous motif cavity window diameters
are calculated using <a href="https://github.com/marcinmiklitz/pywindow">PyWindow</a>.
Solvent sizes are calculated usng <a href="https://github.com/austin-mroz/SMORES">SMORES</a>.
* **size-exclusivity prediction** -- size exclusivity is predicted from the
size analyses.
</td>
</tr>
</table>

## installation

### from environment.yml

```bash

foo@bar$ conda env create -n splashd python=3.12
foo@bar$ conda activate splashd
foo@bar$ pip install .
foo@bar$ pip install git+https://github.com/austin-mroz/SMORES.git@implement-atomlite
foo@bar$ conda install -c conda-forge openbabel
foo@bar$ conda install -c conda-forge xtb
foo@bar$ pip install git+https://github.com/austin-mroz/stko.git

```

### from scratch

## examples

## citation

## acknowledgements

## license

We release this software under the conditions of the MIT license.
