:- use_module(csv_load_mod).
:- use_module(library(clpfd)).
:- use_module(library(clpb)).
:- use_module('reif.pl').
:- set_prolog_stack(global, limit(10*10**9)).
:- dynamic memo_/1.

memo(Goal) :-
    (    memo_(Goal) -> true
    ;    Goal,
         assertz(memo_(Goal))
    ).


:-[prolog_terms_for_bags]. %instance/3.
:-[otu_classes].
:-[microbiome_features].

%These are also in csv_driver file but here we dont need the rest of it.
load_otu_bacteria:-
	context_module(CM),
	prepare_db(CM,'Microbiome/otu_bacteria.csv',taxonmy,otu_tax).

load_img:-
	context_module(CM),
	prepare_db(CM,'IMG/IMG_1_and_2.tsv',img_fields,img_data).

%You need to load these at the begining but slow make with here.
%:-load_otu_bacteria.
%:-load_img.
/* Dcg predicates */
list([])     --> [].
list([L|Ls]) --> [L], list(Ls).

concatenation([]) --> [].
concatenation([List|Lists]) -->
        list(List),
        concatenation(Lists).

bool01_t(1,true).
bool01_t(0,false).

#>=(X,Y,Truth) :- X #>= Y #<==> B, bool01_t(B,Truth).

numlist(0,[]).
numlist(Top,[Top|Rest]):-
    Top #>0,
    Next #= Top -1,
    numlist(Next,Rest).

%!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




number_of_crossoverpoints(3).
mutation_prob(0.05).



otu_tax_list(OTU,Taxlist):-
	otu_tax(OTU,'Phylum',Phylum),
	otu_tax(OTU,'Class',Class),
	otu_tax(OTU,'Order',Order),
	otu_tax(OTU,'Family',Family),
	otu_tax(OTU,'Genus',Genus),
	otu_tax(OTU,'Species',Species),
	Taxlist = [Phylum,Class,Order,Family,Genus,Species].


posOtus(P):-
	setof(Otu,T^instance(Otu,T,les),P).

negOtus(N):-
	setof(Otu,T^instance(Otu,T,nonles),NonLes),
        setof(Otu,T^instance(Otu,T,both),Both),
	once(append(NonLes,Both,N)).

f(Id,f(At,Op,V)):-
	f(Id,Op,At,V).

numberbiggerthan1_t(Number,Value):-
    if_(#>=(Number,1),Value=1,Value=0).

heuristic_chromosome_fitness(rate_diff,C,Fitness):-
        memo(posOtus(PosOtus)),
        memo(negOtus(NegOtus)),
        maplist(f,C,Features),
        length(PosOtus,PositiveLength),
        length(NegOtus,NegativeLength),
        otus_features_filtered(PosOtus,Features,TruePosOtus),
        otus_features_filtered(NegOtus,Features,FalsePosOtus),
        length(TruePosOtus,TruePositives),
        length(FalsePosOtus,FalsePositives),
        Fitness  is (TruePositives/PositiveLength)-(FalsePositives/NegativeLength).

% This predicate filters the list of otus on each recursion for each
% feature
otus_features_filtered([],_,[]).
otus_features_filtered(Otus,[],Otus).
otus_features_filtered(Otus,Features,Filterd):-
    Features =[H|Rest],
    otupartition(H,Otus,Filterd0,_),
    length(Otus,L),
    writeln("NumberofOtus:"),
    writeln(L),
    length(Filterd0,L1),
    %instance(H,_,Class),
    writeln("Class:"),
    writeln(Class),
    writeln("NumberOfTrues:"),
    writeln(L1),
    otus_features_filtered(Filterd0,Rest,Filterd).


otupartition(Feature,List,Ts,Fs):-
    otupartition_pos_ts_fs_feature(List,Ts,Fs,Feature).

otupartition_pos_ts_fs_feature([],[],[],_).
otupartition_pos_ts_fs_feature([X|Xs0],Ts,Fs,Feature):-
    Feature = f(At,Op,FValue),
    instance(X,XBag,_),
    length(XBag,BagSize),
    %size_max_bag(BagSize,300,SmallBag),
    %writeln("SmallBagSize"),
    %length(SmallBag,SBL),
    %writeln(SBL),
    if_(bag_feature_t(BagSize,XBag,f(At,Op,FValue)),
        (   Ts=[X|Ts0],Fs=Fs0),
       (    Fs=[X|Fs0],Ts=Ts0)),
    otupartition_pos_ts_fs_feature(Xs0,Ts0,Fs0,f(At,Op,FValue)).


individual_atvalues(I,AtValues):-
    findall(A-V,img_data(I,A,V),AtValues).

individualcovered_feature_t(I,f(At,_Op,V),T):-
    memo(individual_atvalues(I,AtValueList)),
    if_(memberd_t(At-V,AtValueList),T=true,T=false).


bag_feature_t(_Size,[],_,false).
bag_feature_t(Size,[OneFromBag|RestOfBag],F,Truth):-
    length(RestOfBag,L),
    %writeln(L),
    if_(individualcovered_feature_t(OneFromBag,F),(writeln(Size-L),Truth=true),bag_feature_t(Size,RestOfBag,F,Truth)).


size_max_bag(S,M,B):-
    M #=<S,
    length(B,M).
size_max_bag(S,M,B):-
    S #< M,
    length(B,S).

%!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% splitlocs is the location of crossover points, Length is length of
% chromeosome, recurssion and cut to make sure we get an answer.
splitlocs_length(Sorted,L):-
    number_of_crossoverpoints(N),
    length(C,N),
    maplist(random_between(1,L),C),
    all_different(C),
    sort(C,Sorted),!.

splitlocs_length(S,L):-
    splitlocs_length(S,L).


x_y_xypair(X,Y,X-Y).

my_get_assoc(A,K,V):-
    get_assoc(K,A,V).
assoc_keys_values(A,K,V):-
    maplist(my_get_assoc(A),K,V).

%cross over
mum_dad_c1_c2(M,D,C1,C2):-
    length(M,LM),
    numlist(1,LM,Index),
    maplist(x_y_xypair,Index,M,MPair),
    maplist(x_y_xypair,Index,D,DPair),
    splitlocs_length(SplitLocs,LM),
    list_splitlocs_splits(Index,SplitLocs,SplitIndexes),
    list_to_assoc(MPair,MAssoc),
    list_to_assoc(DPair,DAssoc),
    maplist(assoc_keys_values(MAssoc),SplitIndexes,MumSplit),
    maplist(assoc_keys_values(DAssoc),SplitIndexes,DadSplit),
    s_L1_L2_LA_LB(MumSplit,DadSplit,C1Split,C2Split),
    phrase(concatenation(C1Split), C1),
    phrase(concatenation(C2Split),C2),!.


%split locs should be ordered.
%List is the list of indexes to the list
list_splitlocs_splits(L,SplitLocs,S):-
    list_splitlocs_splits(L,SplitLocs,S,1),!.

list_splitlocs_splits(L,[],[L],_).
list_splitlocs_splits(L,SplitLocs,[[H|Split1]|RestSplits],Count):-
        L=[H|T],
        SplitLocs =[First|Rest],
        H #< First,
        Count2 #= Count+1,
        list_splitlocs_splits(T,[First|Rest],[Split1|RestSplits],Count2).
list_splitlocs_splits(L,SplitLocs,[[H]|RestSplits],Count):-
        L=[H|T],
        SplitLocs =[First|Rest],
        H = First,
        Count2 #= Count+1,
        list_splitlocs_splits(T,Rest,RestSplits,Count2).

numberoffeatures(4072).


list_listindexed(List,ListIndexed):-
    length(List,Len),
    numlist(Len,Indexes),
    maplist(x_y_xypair,List,Indexes,ListIndexed).


newvalue_index_old1_new1(New,Index,_-Index,New-Index).
newvalue_index_old1_new1(_New,Index2,Old-Index,Old-Index):-
    dif(Index2,Index).


old_index_value_new(Old,Index,Value,New2):-
    list_listindexed(Old,OldIndexed),
    maplist(newvalue_index_old1_new1(Value,Index),OldIndexed,New),
    maplist(x_y_xypair,New2,_,New).


chromosome_hillclimbmutated(Range,C1,C2a,C_Fitnesses):-
    %Range =4,
    %Chose random allele
    list_listindexed(C1,C1Index),
    random_member(Feature-Index,C1Index),
    Value #= Feature,
    numberoffeatures(NumberOfFeatures),
    findall(NewValue,(
                Max #=NumberOfFeatures,
                Min #=1,
                Top #=Value+Range,
                Bottom #= Value-Range,

                NewValue in Min..Max,
                NewValue #=<Top,
                NewValue #\=Value,
                NewValue #>=Bottom,
                label([NewValue])
            ),
                NewValues),

    maplist(old_index_value_new(C1,Index),NewValues,NewCs),
    maplist(heuristic_chromosome_fitness(rate_diff),NewCs,ChildFitness),
    maplist(x_y_xypair,ChildFitness,NewCs,C_Fitnesses),

    pop_sortedpop(C_Fitnesses,[C2|Rest]),
    writeln('Mutated Sorted:'),
    maplist(portray_clause,[C2|Rest]),
    C2 = _ValueOfC2-C2a.
    %numlist +-25 of value %if values goout of range ignore.
    %maplist(heuristic_chromosome_fitness blah). %should this look at improvements in context or indvidually. I think in context.


% these have to many choice points. Work out largest and smallest and
% insure these are the last value?
value_max_count_maxfactor_tops(Value,Max,Count,MaxFactor,[NewValue|Tops]):-
    Count #=<MaxFactor,
    NewValue #= Value + (10^Count),
    NewValue #< Max,
    Count1 #= Count +1,
    value_max_count_maxfactor_tops(Value,Max,Count1,MaxFactor,Tops).

value_max_count_maxfactor_tops(_,_,_,_,[]).


value_min_count_maxfactor_bottoms(Value,Min,Count,MaxFactor,[NewValue|Bot]):-
    Count #=<MaxFactor,
    NewValue #= Value - (10^Count),
    NewValue #> Min,
    Count1 #= Count +1,
    value_min_count_maxfactor_bottoms(Value,Min,Count1,MaxFactor,Bot).

value_min_count_maxfactor_bottoms(_,_,_,_,[]).

mylength(S,L):-
    length(L,S).
%Generate a population
%Have constraints and use random labelling.
population_size_len_type(P,Size,L,binary):-
    length(P,Size),
    maplist(mylength(L),P),
    l_randomseeds(Size,Seeds),
    maplist(random_labeling,Seeds,P).


population_size_len_type(P,S,L,number(Min-Max)):-
    length(P,S),
    maplist(random_list(L,Min,Max),P).


random_list(Len,Min,Max,List):-
    length(List,Len),
    maplist(random_between(Min,Max),List).


l_randomseeds(L,S):-
    length(S,L),
    maplist(random_between(1,1000),S).

%recalculating fitnesses.
pop_sortedpop(P,SPR):-
    keysort(P,SP),
    reverse(SP,SPR).


%tournament selection
my_random_member(List,M):-
    random_member(M,List).

pop_tournementselected(P,S):-
    length(TournementMembers,6),
    maplist(my_random_member(P),TournementMembers),
    %should I make this a set so that we dont have duplicates If we do we need to check not singular?
    writeln('Tournement members:'),
    maplist(portray_clause,TournementMembers),
    %pop_sortedpop(TournementMembers,[Best-_X,SecondBest-_Y|_Rest]),
    pop_sortedpop(TournementMembers,[_X-Best,_Y-SecondBest|_Rest]),
    S=[Best,SecondBest].

s_L1_L2_LA_LB(L1,L2,LA,LB):-
    interweave_L1_L2_LA_LB(odd,L1,L2,LA,LB),!.

interweave_L1_L2_LA_LB(_,[],[],[],[]).
interweave_L1_L2_LA_LB(odd,[L1_H|L1_T],[L2_H|L2_T],[L2_H|LA_T],[L1_H|LB_T]):-
    interweave_L1_L2_LA_LB(even,L1_T,L2_T,LA_T,LB_T).
interweave_L1_L2_LA_LB(even,[L1_H|L1_T],[L2_H|L2_T],[L1_H|LA_T],[L2_H|LB_T]):-
    interweave_L1_L2_LA_LB(odd,L1_T,L2_T,LA_T,LB_T).

%Use this to generate our base line performance.
random_population_evaluation(PopSize,GenomeLength,Type,P_Fitnesses):-
    population_size_len_type(P,PopSize,GenomeLength,Type),
    writeln("Init pop:"),
    maplist(heuristic_chromosome_fitness(rate_diff),P,Fitnesses),
    maplist(x_y_xypair,Fitnesses,P,P_Fitnesses),
    maplist(portray_clause,P_Fitnesses).





%type is either `binary` or `number(Min-Max)`
start(Generations,PopSize, GenomeLength,NewBest,Type):-
    %poscpgs(PosCpgs), %remember this is wrong -- why!!
    %negcpgs(NegCpgs),
    population_size_len_type(P,PopSize,GenomeLength,Type),
    writeln("Init pop:"),
    maplist(heuristic_chromosome_fitness(rate_diff),P,Fitnesses,Sums),
    maplist(x_y_xypair,Fitnesses,P,P_Fitnesses),
    maplist(x_y_xypair,Sums,P,P_Sums),
    maplist(portray_clause,P_Fitnesses),
    writeln("Sums:"),
    maplist(portray_clause,P_Sums),
    pop_tournementselected(P_Fitnesses,[Best,SecondBest]),
    mum_dad_c1_c2(Best,SecondBest,C1,C2),
    maplist(heuristic_chromosome_fitness(rate_diff),[C1,C2],ChildFitness,ChildSums),
    maplist(x_y_xypair,ChildFitness,[C1,C2],C_Fitnesses),
    append(P_Fitnesses,C_Fitnesses,NewPop),
    pop_sortedpop(NewPop,Sorted),
    writeln('Sorted pop:'),
    maplist(portray_clause,Sorted),
    writeln('Sorted pop:'),
    maplist(portray_clause,Sorted),

    phrase(concatenation([Middle,[_Last1,_Last2]]),Sorted), %Kill Worst two
    Generation2 #= Generations -1,
    run(Generation2,Middle,NewBest),!.

run(0,[Best|_Rest],Best).
run(Generations,Population,NewBest):-
    pop_tournementselected(Population,[Best,SecondBest]),
    mum_dad_c1_c2(Best,SecondBest,C1a,C2),
    chromosome_hillclimbmutated(4,C1a,C1b,_C_FitnessesX), % Mutate Child1?
    maplist(heuristic_chromosome_fitness(rate_diff),[C1b,C2],ChildFitness),
    maplist(x_y_xypair,ChildFitness,[C1b,C2],C_Fitnesses),
    append(Population,C_Fitnesses,NewPop),
    pop_sortedpop(NewPop,Sorted),
    phrase(concatenation([Middle,[_Last1,_Last2]]),Sorted), %Kill last two
    Generation2 #= Generations -1,
    format("Generation ~w~n",[Generation2]),
    maplist(portray_clause,Sorted),
    run(Generation2,Middle,NewBest).

