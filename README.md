## About This Repo 

Hi! This repository contains the code for my personal website, built using [Quarto](https://quarto.org). 

The website is developed with [Netlify](https://app.netlify.com/), hosted at [ashlynnlu.netlify.app](https://ashlynnlu.netlify.app).


The files are structured as follows:

## _quarto.yml

This .yml file contains the overall structure for my webpage. 

## index.qmd

The [index.qmd](index.qmd) Quarto document contains the structure for the homepage, with a short introduction of myself.

## About

The about page [about.qmd](about.qmd) contains an introduction of myself. 

## posts/

This directory contains the code for all my statistical writings/article shown in the 'Writing' page. Format can be found in [writing.qmd](writing.qmd).

Inside each files, there are:

### index.qmd 
The Quarto document to produce the article. Each post can be previewed by running 'quarto preview "./posts/[POST NAME]/index.qmd"' in terminal.

### data/

This directory contains all raw and derived datasets. Further information for each dataset can be found in metadata.txt inside the file.

The derived data are all processed from the raw data by running the R files in the src/data-cleaning file. 

### src/

This directory contains the data cleaning process and helper functions. 

- The file `data-cleaning\` involves two R files which cleans the raw data and produce derived data.
- The file `helper-function\` include functions aid for analysis/plotting.

### images/

This directory contains the figures used in each article, and the code used to produce the plots. 


## photo/

This directory contains my photography shown in 'Photo' page. Format can be found in [photo.qmd](photo.qmd). The [lightbox extension](https://github.com/quarto-ext/lightbox) is used for the layout.


## To reproduce

Clone or download this repo.

The webpage can be preview locally by running 'quarto preview "./index.qmd"' in terminal.

