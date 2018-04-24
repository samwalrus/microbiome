:-dynamic(class/2).
%:-consult(otu_classes).

otu_tax_list(OTU,Taxlist):-
	%otu_tax(OTU,'Kingdom',Kingdom), %not needed
	otu_tax(OTU,'Phylum',Phylum),
	otu_tax(OTU,'Class',Class),
	otu_tax(OTU,'Order',Order),
	otu_tax(OTU,'Family',Family),
	otu_tax(OTU,'Genus',Genus),
	otu_tax(OTU,'Species',Species),
	Taxlist = [Phylum,Class,Order,Family,Genus,Species].

strip_usless_tax(Full, Without):-
	once(append(Without,[''|_],Full)).

ncbi_tax_list(NCBI,Taxlist):-
	img_data(Img,'Phylum',Phylum),
	img_data(Img,'Class',Class),
	img_data(Img,'Order',Order),
	img_data(Img,'Family',Family),
	img_data(Img,'Genus',Genus),
	img_data(Img,'Species',Spec),
	img_data(Img,'NCBI Taxon ID',NCBI),
	Taxlist =[Phylum,Class,Order,Family,Genus,Spec].

otu_ncbi(Otu,Ncbi,Tax):-
	otu_tax_list(Otu,TaxFull),
	strip_usless_tax(TaxFull,Tax),
	%length(Tax,TaxLength),
	%length(Matcher,TaxLength),
	ncbi_tax_list(Ncbi,NCBITax),
	once(append(Tax,_Tail,NCBITax)).


otu_cogs(Otu,Ncbi,Cogs):-
	otu_ncbi(Otu,Ncbi,_Tax),
	ncbi_cogsets(Ncbi,_NumberOfCogs,Cogs).

otu_cogs_truth(Otu,Ncbi,CogsTruth):-
	otu_cogs(Otu,Ncbi,Cogs),
	all_cogs(AllCogs),
	maplist(logical_purity:list_memberd_truth(Cogs),AllCogs,CogsTruth).

write_header(Stream):-
	format(Stream,"Otu, Ncbi,",[]),
	all_cogs(AllCogs),
	forall(member(C,AllCogs),format(Stream,"~w,",[C])),
	format(Stream,"Class ~n",[]).

write_data_line(Stream,Otu,Ncbi,CogsTruth,Class):-
	format(Stream,"~w, ~w, ",[Otu,Ncbi]),
	forall(member(C,CogsTruth), format(Stream,"~w,",[C])),
	%class(Otu,Class),
	format(Stream,"~w ~n",[Class]).

write_data_line(Stream,Otu,Ncbi,CogsTruth,Class):-
	format(Stream,"~nFailed with ~w ~w ~n~n~n", [Otu,Ncbi]).

%write_data:-
	%forall(class(Otu,Class),(otu_cogs_truth(Otu,Ncbi,CogsTruth), write_data_line(Otu,Ncbi,CogsTruth,Class))).

write_data2(Stream):-
	forall((class(Otu,Class),otu_cogs_truth(Otu,Ncbi,CogsTruth)),write_data_line(Stream,Otu,Ncbi,CogsTruth,Class)).

write_both(File):-
	setup_call_cleanup(open(File,write,Out,[encoding(utf8)]),
	write_header(Out),
	write_data2(Out)),
	close(Out).


class_write(File):-
	otus(Otus,_),
	setup_call_cleanup(open(File, write, Out, [encoding(utf8)]),
			   forall(member(Otu,Otus),
				  (   otu_class(Otu,Class),
				      format("class('~w',~w).~n ",[Otu,Class]),

				  format(Out,"class('~w',~w).~n ",[Otu,Class])
				  )
				 ),
			   close(Out)).


assert_class:-
	otus(Otus,_),
	forall(member(Otu,Otus), (otu_class(Otu,Class),format("~w ~w ~n",[Otu,Class]))).
%class(X,unknown).

otu_class(Otu,les):-
	otus_in_les(OTUS_Les,_),
	otus_in_nonles(OTUS_nonLes,_),
	ord_intersection(OTUS_Les,OTUS_nonLes,_Inter,JustLes),
	member(Otu,JustLes).

otu_class(Otu,both):-
	otus_in_les(OTUS_Les,_),
	otus_in_nonles(OTUS_nonLes,_),
	ord_intersection(OTUS_Les,OTUS_nonLes,Inter,_JustLes),
	member(Otu,Inter).
otu_class(Otu,nonles):-
	otus_in_les(OTUS_Les,_),
	otus_in_nonles(OTUS_nonLes,_),
	ord_intersection(OTUS_nonLes,OTUS_Les,_Inter,JustNonLes),
	member(Otu,JustNonLes).




%maplist(logical_purity:list_memberd_truth([c1,c3,c4]),[c1,c2,c3,c4],Q).
%Q = [true, false, true, true].
%maplist(logical_purity:list_memberd_truth(CogsInOneOtu),AllCogs,Q).



otu_test(Prop,Val):-
	otu_tax(69664,Prop,Val).


spec_test(X,Spec):-
	img_data(X,'Phylum','Firmicutes'),
	img_data(X,'Class','Clostridia'),
	img_data(X,'Order','Clostridiales'),
	img_data(X,'Family','Clostridiaceae'),
	img_data(X,'Genus','Clostridium'),
	img_data(X,'Species',Spec).

imgId_tax(X,[P,C,O,F,G,S]):-
	img_data(X,'Phylum',P),
	img_data(X,'Class',C),
	img_data(X,'Order',O),
	img_data(X,'Family',F),
	img_data(X,'Genus',G),
	img_data(X,'Species',S).


otu_otutax_img(Otu,_OtuTax,_ImgId,_Max):-
        otu_mostspecifictax(Otu,noClassification),
	fail. % We dont want to general cases. maybe change from fail or do some other way.
otu_otutax_img(Otu,OtuTax,ImgId,_Max):-
	otu_mostspecifictax(Otu,kingdom),
	fail.
otu_otutax_img(Otu,OtuTax,ImgId,phylum):-
	otu_mostspecifictax(Otu,phylum),
	otu_tax(Otu,'Phylum',P),
	imgId_tax(ImgId,[P,C,O,F,G,S]),
	OtuTax=[P,C,O,F,G,S].

otu_otutax_img(Otu,OtuTax,ImgId,class):-
	otu_mostspecifictax(Otu,class),
	otu_tax(Otu,'Phylum',P),
	otu_tax(Otu,'Class',C),
	imgId_tax(ImgId,[P,C,O,F,G,S]),
	OtuTax=[P,C,O,F,G,S].

otu_otutax_img(Otu,OtuTax,ImgId,order):-
	otu_mostspecifictax(Otu,order),
	otu_tax(Otu,'Phylum',P),
	otu_tax(Otu,'Class',C),
	otu_tax(Otu,'Order',O),
	imgId_tax(ImgId,[P,C,O,F,G,S]),
	OtuTax=[P,C,O,F,G,S].

otu_otutax_img(Otu,OtuTax,ImgId,family):-
	otu_mostspecifictax(Otu,family),
	otu_tax(Otu,'Phylum',P),
	otu_tax(Otu,'Class',C),
	otu_tax(Otu,'Order',O),
	otu_tax(Otu,'Family',F),
	imgId_tax(ImgId,[P,C,O,F,G,S]),
	OtuTax=[P,C,O,F,G,S].

otu_otutax_img(Otu,OtuTax,ImgId,genus):-
	otu_mostspecifictax(Otu,genus),
	otu_tax(Otu,'Phylum',P),
	otu_tax(Otu,'Class',C),
	otu_tax(Otu,'Order',O),
	otu_tax(Otu,'Family',F),
	otu_tax(Otu,'Genus',G),
	imgId_tax(ImgId,[P,C,O,F,G,S]),
	OtuTax=[P,C,O,F,G,S].


otu_otutax_img(Otu,OtuTax,ImgId,species):-
	otu_mostspecifictax(Otu,species),
	otu_tax(Otu,'Phylum',P),
	otu_tax(Otu,'Class',C),
	otu_tax(Otu,'Order',O),
	otu_tax(Otu,'Family',F),
	otu_tax(Otu,'Genus',G),
	otu_tax(Otu,'Species',S),
	imgId_tax(ImgId,[P,C,O,F,G,S]),
	OtuTax=[P,C,O,F,G,S].






q1(Xs,Y,Z,L):-setof(X,Z^(microbiome(X,Y,Z), Z>0),Xs), length(Xs,L).

q2:-setof(X,Y^Z^(microbiome(X,Y,Z), Z>0),Xs), length(Xs,L).



otus(OTUS,Count):-
	setof(OTU,
	      State^Count^Mars^(clinical(Mars,lesional,State),
			  microbiome(OTU,Mars,Count),Count>1),
	      OTUS),length(OTUS,Count).


otus_in_les(OTUS,Count):-
	setof(OTU,
	      Count^Mars^(clinical(Mars,lesional,'LES'),
			  microbiome(OTU,Mars,Count),Count>1),
	      OTUS),length(OTUS,Count).
otus_in_nonles(OTUS,Count):-
	setof(OTU,
	      Count^Mars^(clinical(Mars,lesional,'NON_LES'),
			  microbiome(OTU,Mars,Count),Count>1),
	      OTUS),length(OTUS,Count).


ncbi_cog(Ncbi,Cog):-
	code_ncbi_name(Code,taxid,Ncbi),
	code_ncbi_name(Code,'FTP-name',Name),
	cog_data(Domain_Id,'genome-name',Name),
	cog_data(Domain_Id,'cog-id',Cog).

ncbi_cogsets(N,L,Cs):-
	setof(C,ncbi_cog(N,C),Cs),
	length(Cs,L).

all_cogs(Cogs):-
	setof(Cog,D^cog_data(D,'cog-id',Cog),Cogs).


otu_mostspecifictax(Otu,noClassification):-
	otu_tax(Otu,'Kingdom',blank).

otu_mostspecifictax(Otu,kingdom):-
	otu_tax(Otu,'Kingdom',K),
	dif(K,blank),
	otu_tax(Otu,'Phylum',blank).

otu_mostspecifictax(Otu,phylum):-
	otu_tax(Otu,'Kingdom',K),
	dif(K,blank),
	otu_tax(Otu,'Phylum',P),
	dif(P,blank),
	otu_tax(Otu,'Class',blank).


otu_mostspecifictax(Otu,class):-
	otu_tax(Otu,'Kingdom',K),
	dif(K,blank),
	otu_tax(Otu,'Phylum',P),
	dif(P,blank),
	otu_tax(Otu,'Class',C),
	dif(C,blank),
	otu_tax(Otu,'Order',blank).

otu_mostspecifictax(Otu,order):-
	otu_tax(Otu,'Kingdom',K),
	dif(K,blank),
	otu_tax(Otu,'Phylum',P),
	dif(P,blank),
	otu_tax(Otu,'Class',C),
	dif(C,blank),
	otu_tax(Otu,'Order',O),
	dif(O,blank),
	otu_tax(Otu,'Family',blank).

otu_mostspecifictax(Otu,family):-
	otu_tax(Otu,'Kingdom',K),
	dif(K,blank),
	otu_tax(Otu,'Phylum',P),
	dif(P,blank),
	otu_tax(Otu,'Class',C),
	dif(C,blank),
	otu_tax(Otu,'Order',O),
	dif(O,blank),
	otu_tax(Otu,'Family',F),
	dif(F,blank),
	otu_tax(Otu,'Genus',blank).

otu_mostspecifictax(Otu,genus):-
	otu_tax(Otu,'Kingdom',K),
	dif(K,blank),
	otu_tax(Otu,'Phylum',P),
	dif(P,blank),
	otu_tax(Otu,'Class',C),
	dif(C,blank),
	otu_tax(Otu,'Order',O),
	dif(O,blank),
	otu_tax(Otu,'Family',F),
	dif(F,blank),
	otu_tax(Otu,'Genus',G),
	dif(G,blank),
	otu_tax(Otu,'Species',blank).

otu_mostspecifictax(Otu,species):-
	otu_tax(Otu,'Kingdom',K),
	dif(K,blank),
	otu_tax(Otu,'Phylum',P),
	dif(P,blank),
	otu_tax(Otu,'Class',C),
	dif(C,blank),
	otu_tax(Otu,'Order',O),
	dif(O,blank),
	otu_tax(Otu,'Family',F),
	dif(F,blank),
	otu_tax(Otu,'Genus',G),
	dif(G,blank),
	otu_tax(Otu,'Species',S),
	dif(S,blank).



family_bag(OTU,NCBI,Family):-
	otu_tax(OTU,'Family',Family),
	img_data(IMG_ID,'Family',Family),
	img_data(IMG_ID,'NCBI Taxon ID',NCBI).

q1:- family_bag(Otu,Ncbi,Family), ncbi_cogsets(Ncbi,SetSize,Cs).

gaz(IMG_ID,NCBI_ID,Gram):-
	img_data(IMG_ID, 'NCBI Taxon ID', NCBI_ID),
	img_data(IMG_ID, 'Gram Staining', Gram).


gaz_write(File):-
	setup_call_cleanup(open(File, write, Out, [encoding(utf8)]),
			   forall((otu_ncbi(OTU,Ncbi,Tax), gaz(_Img,Ncbi,Gram)),
				  format(Out,"OTU=~w, Ncbi=~w, tax=~w, Gram=~w, ~n",[OTU,Ncbi,Tax,Gram])),
			   close(Out)).






%genus bags
write_bagsizes(File,Class):-
	setup_call_cleanup(open(File, write, Out, [encoding(utf8)]),
			   forall((otu_class(Otu,Class),otu_mostspecifictax(Otu,T),setof(Img,Otutax^otu_otutax_img(Otu,Otutax,Img,T),Imgs),length(Imgs,L)),
				  format(Out,"OTU=~w, BagSize,~w, Class,~w,MostspecificTax,~w~n",[Otu,L,Class,T])),
			   close(Out)).
