%% Population und Generationen
g_max = 20;        % Maximale Anzahl der Generationen
p = 50;             % Populationsgröße
% Nach der Umweltselektion bleiben immer p Individuen: Neue Eltern
% Aus diesen p Eltern werden immer p Kinder erzeugt



%% Parameter der Fitnessfunktion
anz_eig = 1;            % hier Anzahl der Eignungsgrade angeben


% Optimierung bei gegebener Taktzeit
t_takt = 0;             % Wenn hier eine Zeit t_takt ~= 0 eingetragen ist, entfällt die Taktzeit in der Fitnessfunktion
                        

% Gewichtung: Gesamtfitness, es können beliebige Zahlen eingegeben werden
w(1,1) = 0.5;       % Gewichtung Takzeit
w(2,1) = 0;       % Gewichtung Auslastung

w(5,1) = 0.5;       % Gewichtung 1. Eignungsgrad
w(6,1) = 0;       % Gewichtung 2. Eignungsgrad
w(7,1) = 0;       % Gewichtung 3. Eignungsgrad
w(8,1) = 0;       % Gewichtung 4. Eignungsgrad
w(9,1) = 0;       % Gewichtung 5. Eignungsgrad

% Werkzeuge und Pausen ein- und ausschalten
w(3,1) = 0;       % Gewichtung Pausen: 1(ein) oder 0(aus)
w(4,1) = 0;       % Gewichtung Werkzeuge: 1(ein) oder 0(aus)

% Auslastung: Verteilung
a(1,1) = 0.5;       % Mensch
a(2,1) = 0.5;       % Roboter

% Pausen: Verteilung und maximale Belohung
max_belohn_pausen = 0.2;
b(1,1) = 0.5;       % Mensch, wenige Pausen gut
b(2,1) = 0.5;       % Roboter, möglichst verteilt

% Pausen: Verteilung und maximale Belohung
max_belohn_werkzeuge = 0.2;
t(1,1) = 0.5;       % Mensch, wenige Pausen gut
t(2,1) = 0.5;       % Roboter, möglichst verteilt

% Eignung 1: Verteilung
e(1,1) = 0.5;
e(2,1) = 0.5;

% Eignung 2: Verteilung
e2(1,1) = 0.5;
e2(2,1) = 0.5;

% Eignung 3: Verteilung
e3(1,1) = 0.5;
e3(2,1) = 0.5;

% Eignung 4: Verteilung
e4(1,1) = 0.5;
e4(2,1) = 0.5;

% Eignung 5: Verteilung
e5(1,1) = 0.5;
e5(2,1) = 0.5;

% FÜR ERGONOMIE: Verteilung des Eignungsgrades für den Roboter auf 0 und für
% den Menschen auf 1 setzen


%% Abbruchbedingungen

% Abbbruchbedingnung 1: Keine Verbesserung des Champion
abbruch_typ1 = 0;           % Einschalten (1) und Ausschalten (0) 
abbruch_gen = 30;           % Differenz zwischen letztem Champion und Vergleichschampion
abbruch_wert = 0.01;        % Minimal zulässige Fitnessdifferenz der Vergleichschampion

% Abbruchbedingung 2: Zeitüberschreitung
abbruch_typ2 = 0;           % Einschalten (1) und Ausschalten (0)            
abbruch_time = 2;           % Maximale Zeit in Minuten 


%% Parameter der Decodierung
% Art der Decodierung
dec_typ = 2;
    % 1: Semi-aktive Zeitpläne
    % 2: Aktive und Non-Delay Zeitpläne

% für aktive und Non-Delay Zeitpläne: sigma festlegen
dc_sigma = 0.5; % Wertebereiche [0 1], 0 entspricht Non-Delay, 1 enspricht aktiv


%% Parameter des genetischen Algorithmus
% ERSETZUNGSSTRATEGIE
%   Hier kann festgelegt werden, welche und wie viele Individuen in die
%   Umweltselektion gelangen

% Plus- und Kommaselektion: Bestandteile von PK festlegen
mu_el = 1.0;          % Anteil der Eltern für die Umweltselektion
mu_ki = 1.0;          % Anteil der Kinder für die Umweltselektion
% Summe aus p_el und p_ki muss p_el + p_ki >= p sein, sonst zu wenig Individuen für
% nächste Gen. Bei einer Summe von p_el + p_ki = p wird die Selektion Sinnlos


% ELTERNMUTATION
el_mut_whk = 1.0;   % Alle, da nur eine Kopie der Population mutiert wird (sollte immer so bleiben)

% Auswahl des Verfahrens
%Sequenz:
el_mutation_seq = 1;
    % 1 = Vertauschende Mutation (swap_mutation)
        el_mut_whk_swap = 0.05;
    % 2 = Mischende Mutation (mixed_mutation)(el_mut_whk)
    % 3 = Invertierende Mutation (invert_mutation)(el_mut_whk)
       


%Zuordnung:
el_mutation_allo = 1;
    % 1 = 1-Bit-Mutation (ein_bit_mutation)
        el_mut_whk_einbit = 0.05;
    % 2 = Invertierende Mutation (invert_mutation)(el_mut_whk)
    % 3 = Mischende Mutation (mixed_mutation)(el_mut_whk)
   


% Ersetzung der Eltern durch mutierte Eltern
anteil_mut_eltern = 0.5;        % bei 1.0 werden alle Eltern mutiert (auch der Champion!!!)
                                % bei 0.0 werden die Eltern nicht mutiert    



% ELTERNSELEKTION
% Auswahl der Fitness
fitness_typ = 4;
    % 2: tat. Fitness
    % 4: rangbasierte Fitness 
    % 6: Shared Fitness (muss über Crowding Distance zusätzlich eingestellt
    % werden)


% Auswahl des Verfahrens
eltern_selektion = 1;       % bisher keine Alternative
    % 1: Stochastik Universal Sampling



% CROSSOVER
% Auswahl des Verfahrens

%Sequenz:
crossover_seq = 1;
PermuCrossWhk = 0.4;    % Whk für Rekombination, wenn nicht wird das erstgewählte (zufällig) Elter gewählt
    % 1 = Ordnungscrossover
    % 2 = Partially Mapped Crossover


% Zuordnung:
crossover_allo = 1;
AlloCrossWhk = 0.4;    % % Whk für Rekombination, wenn nicht wird das erstgewählte (zufällig) Elter gewählt
    % 1 = 2-Punkt-Crossover
    % 2 = Uniformer-Crossover
        co_unif_whk = 0.5; % Zusatz nur für Uniformer Crossover     
        % Sollte immer auf 0.5 stehen (keines der Eltern wird bevorzugt)
 
   
        
% MUTATION

% Auswahl des Verfahrens
%Sequenz:
mutation_seq = 1;
    % 1 = Vertauschende Mutation (swap_mutation)
        mut_whk_swap = 0.1;
    % 2 = Mischende Mutation (mixed_mutation)
    % 3 = Invertierende Mutation (invert_mutation)
        PermuMutationWhk = 1;    % für mixed und invert   


%Zuordnung:
mutation_allo = 1;
    % 1 = 1-Bit-Mutation (ein_bit_mutation)
        mut_whk_einbit = 0.1;
    % 2 = Invertierende Mutation (invert_mutation)
    % 3 = Mischende Mutation (mixed_mutation)
        AlloMutationWhk = 1;     % für invert und mixed


% UMWELTSELEKTION
% Auswahl der Fitness
umwelt_fitness = 3;
    % 1: Gesamtfitness
    % 2: Shared Fitness
    % 3: Rang mit Crowding Distance


% Auswahl des Verfahrens
umwelt_selektion = 1;
    % 1: Q-stufige-Tunierselektion, bisher keine Alternative
 
% Ein- und Ausschalten des Champion
champion_confi = 1;
    % 1: Champion an
    % 0: Champion aus

% Anzahl der Tuniere je Generation
q = 3;


% Ränge und Nischentechniken
anz_raenge = 3;     % Anzahl der Ränge (Einfluss auf Eltern- und Umweltselektion)    
R_max = 1.6;        % Empfehlungswert aus Quelle:

% nur relevant wenn Fitness Sharing in einer der Funktionen verwendet wird
hamming_abstand = 5;    % gemittelter Hamming Abstand  
delta_startzeit = 5;    % Zeitabstand einer OperationX(t_start) des Individuum1 zu der OperationX des Individuums2 (Mittelwert)

