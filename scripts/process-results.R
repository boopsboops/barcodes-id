#!/usr/bin/env Rscript
options(width=180)

# load libs
library("here")
library("tidyverse")
library("glue")

# load sample codes
sample.codes <- readr::read_tsv(here::here("temp/sample-codes.tsv"), col_types=list("c","c"), show_col_types=FALSE)

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

# load bold taxonomy
bold.tax <- readr::read_tsv(here::here("temp/source/BOLDistilled_COI_Jan2026_TAXONOMY.tsv"),show_col_types=FALSE) |> 
    dplyr::select(bin,family,species) |> 
    dplyr::rename(blastid=species) |>
    dplyr::rename(binid=bin)

# load sintax results
sintax.cols <- c("sampleid","idsProbs","strand","ids")
sintax.table <- readr::read_tsv(here::here("temp/cox1-sintax-out.tsv"),col_names=sintax.cols,show_col_types=FALSE) |>
    tidyr::separate_wider_delim(cols=idsProbs,delim=",s:",names=c("higherid","speciesid")) |>
    tidyr::separate_wider_delim(cols=speciesid,delim="(",names=c("sintaxid","sintaxBootstrap")) |> 
    dplyr::mutate(sintaxBootstrap=stringr::str_replace_all(sintaxBootstrap,"\\)","")) |> 
    dplyr::select(sampleid,sintaxid,sintaxBootstrap)

# clean and filter blast resuts
blast.results.clean <- blast.results |> dplyr::left_join(bold.tax) |> 
    dplyr::group_by(sampleid) |> 
    dplyr::slice_max(blastBitscore,with_ties=TRUE) |> 
    dplyr::ungroup() |>
    dplyr::mutate(dateSearched=lubridate::today())

# annotate blast and join with sintax and sample ids
blast.results.clean.ann <- blast.results.clean |> 
    dplyr::left_join(sintax.table) |> 
    dplyr::left_join(sample.codes) |>
    dplyr::mutate(sampleid=glue::glue("{sampleid}-{sampleCode}")) |>
    dplyr::mutate(identification=if_else(blastid==sintaxid,blastid,"no identification")) |>
    dplyr::select(sampleid,processid,binid,family,identification,blastEvalue,blastLength,blastPident,sintaxBootstrap,dateSearched)

# check and write out
blast.results.clean.ann |> print(n=Inf)
blast.results.clean.ann |> readr::write_csv(here::here("temp/cox1-id-results.csv"))

# print table for report
blast.results.clean.ann |> dplyr::mutate(output=glue::glue("{sampleid}_FishF1/R1: {identification} (blast length {blastLength}, identity {round(blastPident,digits=2)}, sintax bootstrap {sintaxBootstrap})")) |> pull(output)
