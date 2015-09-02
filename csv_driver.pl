:- use_module(csv_load_mod).


load_clinical :-
	context_module(CM),
	prepare_db(CM,'Clinical/Psoriasis.csv',clinical_factors,clinical).


load_microbiome :-
	context_module(CM),
	prepare_db(CM,'Microbiome/main_table.csv', mirobiome_ids,microbiome).
