:- use_module(library(pio)).
:- use_module(library(dcg/basics)).
:- portray_text(true).
:- set_prolog_flag(double_quotes, codes).
:- set_prolog_flag(back_quotes,string).
:- set_prolog_stack(global, limit(4*10**9)).


:- [library(csv)].

:- dynamic samples/3.
:- dynamic column_keys/1.

prepare_db(File) :-
    retractall(column_keys(_)),
    retractall(samples(_,_,_)),
    ( nonvar(File) ; File = '/tmp/test.csv' ),
    forall(read_row(File, Row), store_row(Row)).

store_row(Row) :-
    Row =.. [row|Cols],
    (   column_keys(ColKeys)
    ->  Cols = [RowKey|Samples],
        maplist(store_sample(RowKey), ColKeys, Samples)
    ;   (Cols = [_H|T],assertz(column_keys(T)))
    ).

store_sample(RowKey, ColKey, Sample) :-
    assertz(samples(RowKey, ColKey, Sample)).


read_row(File, Row) :-
    csv_read_file_row(File, Row, []),
    writeln(read_row(Row)).
