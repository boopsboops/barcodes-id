#!/usr/bin/env Rscript
options(width=180)

# load libs
library("here")
library("tidyverse")
library("glue")
library("traits")

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
blast.results <- readr::read_tsv(here::here("temp/16s-blast-local.out"), col_names=col.names, col_types=col.types, show_col_types=FALSE) |>
    dplyr::mutate(dbid=stringr::str_replace_all(dbid,"\\..+",""))

# clean and filter blast resuts
blast.results.clean <- blast.results |> 
    dplyr::group_by(sampleid) |> 
    dplyr::slice_max(blastBitscore,with_ties=TRUE) |> 
    dplyr::ungroup() 

# load ncbi
ncbi.taxonomy <- traits::ncbi_byid(ids=dplyr::pull(dplyr::distinct(blast.results.clean,dbid),dbid)) |> 
    tibble::as_tibble() |> 
    dplyr::select(taxon,acc_no) |>
    dplyr::rename(blastid=taxon,dbid=acc_no) |>
    dplyr::mutate(dbid=stringr::str_replace_all(dbid,"\\..+",""))

# annotate with taxonomy and dereplicated
blast.results.clean.ann <- blast.results.clean |> 
    dplyr::left_join(ncbi.taxonomy) |>
    dplyr::mutate(blastid=stringr::str_replace_all(blastid,"Gadus finnmarchicus","Gadus chalcogrammus")) |>
    dplyr::mutate(dateSearched=lubridate::today()) |> 
    dplyr::group_by(sampleid,blastEvalue,blastLength,blastPident,blastNident,blastBitscore,blastid,dateSearched) |>
    dplyr::summarise(blastuniq=dplyr::n_distinct(dbid)) |>
    dplyr::ungroup()


# load sintax results
sintax.cols <- c("sampleid","idsProbs","strand","ids")
sintax.table <- readr::read_tsv(here::here("temp/16s-sintax-out.tsv"),col_names=sintax.cols,show_col_types=FALSE) |>
    tidyr::separate_wider_delim(cols=idsProbs,delim=",s:",names=c("higherid","speciesid")) |>
    tidyr::separate_wider_delim(cols=speciesid,delim="(",names=c("sintaxid","sintaxBootstrap")) |> 
    dplyr::mutate(sintaxBootstrap=stringr::str_replace_all(sintaxBootstrap,"\\)","")) |> 
    dplyr::mutate(sintaxid=stringr::str_replace_all(sintaxid,"_.+","")) |> 
    dplyr::select(sampleid,sintaxid,sintaxBootstrap)

# annotate blast and join with sintax and sample ids
blast.results.clean.merged <- blast.results.clean.ann |> 
    dplyr::left_join(sintax.table) |> 
    dplyr::left_join(sample.codes) |>
    dplyr::mutate(sampleid=glue::glue("{sampleid}-{sampleCode}")) |>
    #dplyr::mutate(identification=if_else(blastid==sintaxid,blastid,"no identification")) |>
    dplyr::select(sampleid,blastid,blastEvalue,blastLength,blastPident,blastuniq,sintaxid,sintaxBootstrap,dateSearched)

# check and write out
blast.results.clean.merged |> print(n=Inf)
blast.results.clean.merged |> readr::write_csv(here::here("temp/16s-id-results.csv"))

# print table for report
blast.results.clean.merged |> dplyr::mutate(output=glue::glue("{sampleid}_16sarL/16sbrH: blastid {blastid} (blast length {blastLength}, identity {round(blastPident,digits=2)}, blast unique seqs {blastuniq}), sintaxid {sintaxid} (sintax bootstrap {sintaxBootstrap})")) |> pull(output)
