:- use_module(library(pio)).
:- use_module(library(dcg/basics)).
:- use_module(csv_load_mod).
:- portray_text(true).
:- set_prolog_flag(double_quotes, codes).
:- set_prolog_flag(back_quotes,string).
:- set_prolog_stack(global, limit(15*10**9)).

:-[csv_driver].

:-load_clinical.
:-load_microbiome.
:-load_gene_expression.
:-load_otu_bacteria.
:-load_img.
:-load_cog.
:-load_ncbi_to_cog.

:-edit('debug.pl').
:-edit('csv_load_mod.pl').
:-edit('csv_driver.pl').
:-edit('basic_queries.pl').


