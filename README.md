README
================
Richard Sprague
June 06, 2018

Getting Started
---------------

This package helps you read and manipulate microbiome files in the open source package Phyloseq.

Prerequisites
=============

If this is your first time to use Phyloseq, you will need to load a few prerequisites:

    > install.packages(c(""ggplot2","dplyr", "tidyr", "readxl"))
    > source("https://bioconductor.org/biocLite.R")
    > biocLite("phyloseq")

Once these are installed, you’ll need the `actino` package, which you can download from this repo using the devtools package:

    > library(devtools)
    > install_github("richardsprague/actino")

Quick Start
===========

Actino includes a few built-in phyloseq objects, including `kombucha.genus`. The following commands will plot the abundances:

    > library("actino")
    > kombucha.csv #  a pre-loaded CSV 
    > plot_bar(prune_taxa(taxa_sums(kombucha.genus)>100000,kombucha.genus), fill = "Genus")

``` r
library("actino")
library(phyloseq)
plot_bar(prune_taxa(taxa_sums(kombucha.genus)>100000,kombucha.genus), fill = "Genus")
```

![](README_files/figure-markdown_github/plotKombucha-1.png)

To load a new phyloseq object from a JSON file you will need a mapfile, which is an Excel file containing a row for each sample, and several columns including "SSR", "Date", "Label" and more. (Figure @ref(mapfile))

| Username        |    SSR| Label   | Date       |
|:----------------|------:|:--------|:-----------|
| Richard Sprague |   7423| Label 1 | 2018-01-01 |
| Richard Sprague |  55309| Label 2 | 2018-01-02 |
| Richard Sprague |   9742| Label 3 | 2018-01-03 |
| Richard Sprague |  55669| Label 4 | 2018-01-04 |
| Richard Sprague |  78762| Label 5 | 2018-01-05 |

If your data is in a CSV file, you'll need it in a form with the following columns

    > experiment_to_phyloseq("myCSVfile.csv",mapfile)

This will generate a valid phyloseq object.

Load some data
--------------

The package includes some sample data, which you can load like this:

``` r
library(actino)
data("kombucha.csv")
head(kombucha.csv)
```

    ##             tax_name     tax_rank count count_norm taxon parent ssr
    ## 1               root         root 77768    1035401     1      0   0
    ## 2           Bacteria superkingdom 75109    1000000     2 131567   0
    ## 3      Campylobacter        genus     2         26   194  72294   0
    ## 4     Flavobacterium        genus    60        798   237  49546   0
    ## 5     Alcaligenaceae       family    38        505   506  80840   0
    ## 6 Enterobacteriaceae       family    87       1158   543  91347   0

The data, loaded into the variable `kombucha.csv` is a dataframe created from a CSV file with three columns: `tax_name`, `ssr`, and `count`. It represents the SSRs associated with a week-long experiment of a user (me!) who drank several liters of kombucha each day to see the effect on his microbiome.

This is the raw data. To be useful for Phyloseq, it needs an accompanying file, called a mapfile, that maps attributes to each sample. The mapfile is a dataframe where the columns correspond to attributes of the data. You can have as many columns as you like, but one of the columns (usually the first one) must stand for the name of the sample. In most uBiome situations, the samples will be referred to by their SSR, and the `actino` package will expect a column named `ssr` in its mapfile. (Note: lowercase)

Let's create a very simple mapfile for `kombucha.csv`:

``` r
kombucha <- kombucha.csv  # simplify the name to make it easier to type

ssrs <- unique(kombucha$ssr)

ssrs
```

    ## [1] 0

As you can see, we have a total of 15 unique SSRs, corresponding to 15 unique samples. Let's construct a very simple map file.

The first column represents the sample name (aka SSR). We'll add a second column that represents the date of the sample. Finally, let's add a third column `color` that we’ll pretend is associated with some feature of the samples. Color of the stool? Color of the kombucha drink that day? Okay, if you don’t like my example features, so go ahead and add some other columns if you like.

``` r
mapfile<-data.frame(ssr=ssrs,date = Sys.Date()+seq(1,length(ssrs)), color = colors()[1:length(ssrs)])

mapfile
```

    ##   ssr       date color
    ## 1   0 2018-06-07 white

Although this example is a bit contrived, the result for Phyloseq is an entirely valid mapfile. Alternatively you could have made the map file manually too; Remember that in R it’s easy to create any dataframe from a CSV or Excel file. If you create one by hand, simply read it into R using the `read.csv` command. Although in this demo, we are loading some pre-existing data, you could just as easily have loaded straight from the CSV file of your choice like this:

    kombucha <- read.csv("myKombuchaDataFile.csv")
    mapfile <- read.csv("myKombuchaMapFile.csv")

References and useful ways to get started
=========================================

The paper that announced it to the world: <http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217>

Here is an excellent beginners guide:

<http://joey711.github.io/phyloseq-demo/phyloseq-demo.html>

How to import from other data formats: <http://joey711.github.io/phyloseq/import-data.html>

Don’t forget to check the [Phyloseq online manual](https://rdrr.io/bioc/phyloseq/): <https://rdrr.io/bioc/phyloseq/>

FAQ: <https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#other-issues-related-the-biom-format>
