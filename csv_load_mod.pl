:-module(csv_load_mod,[prepare_db/4]).
:-use_module(library(csv)).

prepare_db(CM,File, Column_Key,Relation) :-
    Column_Key_Term =.. [Column_Key,_],
    Relation_Term =.. [Relation,_,_,_],
    retractall(CM:Column_Key_Term),
    retractall(CM:Relation_Term),
    forall(read_row(File, Row), store_row(CM,Row,Column_Key,Relation)).

store_row(CM,Row,Column_Key,Relation) :-
    Column_Key_Test =.. [Column_Key,ColKeys],
    Row =.. [row|Cols],

    (   call(CM:Column_Key_Test)
    ->  Cols = [RowKey|Values],
        maplist(store_relation(CM,Relation,RowKey), ColKeys, Values)
    ;   ( Cols = [_H|T],
	Column_Key_Term =.. [Column_Key,T],
       CM:assertz(Column_Key_Term))
    ).

store_relation(CM,Relation,RowKey, ColKey, Values) :-
    Relation_Term =.. [Relation,RowKey,ColKey,Values],
    CM:assertz(Relation_Term).

read_row(File, Row) :-
    csv_read_file_row(File, Row, []).
    %writeln(read_row(Row)).
