:-module(csv_load_mod,[prepare_db/3]).
:-use_module(library(csv)).

prepare_db(File, Column_Key,Relation) :-
    Column_Key_Term =.. [Column_Key,_],
    Relation_Term =.. [Relation,_,_,_],
    retractall(Column_Key_Term),
    retractall(Relation_Term),
    forall(read_row(File, Row), store_row(Row,Column_Key,Relation)).

store_row(Row,Column_Key,Relation) :-
    Column_Key_Test =.. [Column_Key,ColKeys],
    Row =.. [row|Cols],
    (   call(Column_Key_Test)
    ->  Cols = [RowKey|Values],
        maplist(store_relation(Relation,RowKey), ColKeys, Values)
    ;   ( Cols = [_H|T],
	Column_Key_Term =.. [Column_Key,T],
       assertz(Column_Key_Term))
    ).

store_relation(Relation,RowKey, ColKey, Values) :-
    Relation_Term =.. [Relation,RowKey,ColKey,Values],
    assertz(Relation_Term).

read_row(File, Row) :-
    csv_read_file_row(File, Row, []).
    %writeln(read_row(Row)).
