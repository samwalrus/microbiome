:- use_module(csv_load_mod).

%:-load_all.

load_all:-
	load_clinical,
	%load_microbiome,
	%load_gene_expression,
	load_otu_bacteria,
	load_img,
	load_cog,
	load_ncbi_to_cog.

load_clinical :-
	context_module(CM),
	prepare_db(CM,'Clinical/Psoriasis.csv',clinical_factors,clinical).

load_microbiome :-
	context_module(CM),
	prepare_db(CM,'Microbiome/main_table.csv', mirobiome_ids,microbiome).

load_gene_expression :-
	context_module(CM),
	prepare_db(CM,'Affy/affy_table.csv', affy_ids, gene_expression).

load_otu_bacteria:-
	context_module(CM),
	prepare_db(CM,'Microbiome/otu_bacteria.csv',taxonmy,otu_tax).

load_img:-
	context_module(CM),
	prepare_db(CM,'IMG/IMG_1_and_2.tsv',img_fields,img_data).

load_cog:-
	context_module(CM),
	prepare_db(CM,'COG/cog2003-2014.csv',cog_fields,cog_data).

load_ncbi_to_cog:-
	context_module(CM),
	prepare_db(CM,'COG/code_ncbitax_name.csv',ncbi_cog,code_ncbi_name).




