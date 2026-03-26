#!/usr/bin/env Rscript
options(width=180)

# load libs
library("here")
library("tidyverse")
library("glue")


# make cols
blast.cols <- list(
    sampleid="c",
    dbid="c",
    blastEvalue="d",
    blastLength="d",
    blastPident="d",
    blastNident="d",
    blastBitscore="d"
    )

# create col names and types
col.names <- names(blast.cols)
col.types <- do.call(cols, blast.cols)

# load blast results
blast.results <- readr::read_tsv(here::here("temp/cox1-blast-local.out"), col_names=col.names, col_types=col.types, show_col_types=FALSE) |>
    tidyr::separate_wider_delim(dbid,delim="|",names=c("processid","binid"))


# taxonomy
bold.tax <- readr::read_tsv(here::here("temp/source/BOLDistilled_COI_Jan2026_TAXONOMY.tsv"),show_col_types=FALSE) |> 
    dplyr::select(bin,family,species) |> 
    dplyr::rename(binid=bin)

blast.results |> dplyr::left_join(bold.tax) |> 
    dplyr::group_by(sampleid) |> 
    dplyr::slice_max(blastBitscore,with_ties=TRUE) |> 
    dplyr::ungroup() |>
    dplyr::mutate(dateSearched=lubridate::today()) |>
    print(n=Inf)