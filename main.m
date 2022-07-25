clear
clc
close all


%% INITIALISIERUNG
run('config')   % Laden der Konfiguration
run('trans')    % Umrechungen einige Konfigurationen auf Programmvariablen
run('data')     % Auslesen der Daten

inputname = inputdlg('Name der Prozess Datei', 'Prozessname');
prozessname = char(inputname);

% Deklaration einiger Laufvariablen

% z1: Laufvariable der Populationsgröße p
% z2: Laufvariable der Prozessanzahl n
% s: Spalten einer Matrix (meist bis n)


% Matrizen für Plots
fit_plot = zeros(p,g_max);          % Matrix zum Plotten der Fitness (Jedes Individuum)
R_plot = zeros(anz_raenge,g_max);   % Matrix zum Plotten der Ränge 
HA = zeros(p,g_max);                % Matrix zum Plotten des Hamming-Abstands (Jedes Individuum)
CD = zeros(p,g_max);                % Matrix zum Plotten des Zeitabstandes (Jedes Individuum)
champ_plot = zeros(1,g_max);        % Matrix zum Plotten der Champion(Nur der Champion)
eig_plot = zeros(2,g_max);          % Matrix zum Plotten der Eignung(Mensch/Roboter)
ausl_plot = zeros(2,g_max);         % Matrix zum Plotten der Auslastung (Mensch/Roboter)
ges_eig_plot = zeros(p,g_max);      % Matrxi zum Plotten der gesamt Eignung (Jedes Individuum)
ges_ausl_plot = zeros(p,g_max);     % Matrix zum Plotten der gesamt Auslastung (Jedes Individuum)


%% INI: Erstellung der ersten Generation, zufallsgeneriert

g = 1;                  % Initialisierungsgeneration, erhält die Nr.1
P_ini = cell (p,1);     % Population, jede Zelle stellt ein Individuum dar

for z1 = 1:p
    P_ini{z1}(1,:) = randperm(n);
end

for z1 = 1:p
    P_ini{z1}(2,:) = round(rand(1,n));
end

for z1 = 1:p
    P_ini{z1}(3,:) = ones(1,n) - P_ini{z1}(2,:);
end
% 

% Ausgabe der Individuen im struct zur Einsicht / Kontrolle

for position = 1:p
    Population.(['Individuum' num2str(position)]) = struct('Sequenz',P_ini{position}(1,:),'Mensch',P_ini{position}(2,:),'Roboter',P_ini{position}(2,:));
end

P_ausgabe = zeros(3*p,n);

for z1=1:p
    P_ausgabe(z1*3-2,:) = P_ini{z1,1}(1,:);
    P_ausgabe(z1*3-1,:) = P_ini{z1,1}(2,:);
    P_ausgabe(z1*3,:) = P_ini{z1,1}(3,:);
end
   
xlswrite('Population.xls',P_ausgabe)


%% INI: Eingabe von Individuen

x = inputdlg('Soll ein manuell erstelltes Individuum hinzufügt werden? (j/n)', 'Eingabe von Individuen');

position = 1;
while x{1,1}(1,1) ~= 'n'

    [xlsfile,path2xls] = uigetfile('*.xls','Bitte Individuum auswählen');
    fid = fopen (fullfile(path2xls,xlsfile),'r'); 
    [num] = xlsread(xlsfile);

    P_ini{position,1}(1,:) = num(1,:);
    P_ini{position,1}(2,:) = num(2,:);
    P_ini{position,1}(3,:) = num(3,:);
    
    position  = position + 1;   % neue Position für weiteres Individuum
    
    x = inputdlg('Soll ein weiteres manuell erstelltes Individuum hinzufügt werden? (j/n)', 'Eingabe von Individuen');

end


%% INI: Eingabe manuelle Zuordnungen
correct = inputdlg('Sollen Operationen oder Permutationen fest zuordnet werden? (j/n)', 'Feste Zuordnungen');
if correct{1,1} == 'j'
    
    % Eingabe für die Korrektur
    [P_ini,geneM,geneR,permutation] = manuelle_allocation_eingabe(P_ini,n);

    % Bestimmung von festen Werten für Zuordnung und Permutation für ausgewählte Gene
    [P_ini] = manuelle_allocation(P_ini,p,geneM,geneR,permutation);
    
end


%% INI: Eingabe von kollaborativen Operationen
kolla = inputdlg('Sollen Operationen kollaborativ ausgeführt werden? (j/n)', 'Zuordnungsgen');
if kolla{1,1}(1,1) == 'j'
    prompt = {'Operationensnummern (Leerzeichen als Trennung)'};
    title = 'Kollaborative Operationen';
    dims = [1 35];
    answer = inputdlg(prompt,title,dims);
    
    op_kolla(1,:) = str2num(regexprep(answer{1,1}, '\D+', ' '));
else
    op_kolla = zeros(1,1);     % Wenn nicht, dann Nullvektor erstellen
end


%% INI: Plausibilitätstest
% eventulle falsche Eingabe von Individuen

[wrong_seq,wrong_allo] = plausibilitaet (P_ini,p,n);

if ~isempty(wrong_seq)
    msgbox(['Das Sequenzchromosom der Individuen ',num2str(wrong_seq),' der Population P wurde falsch erzeugt'],'Fehler')
end

if ~isempty(wrong_allo)
    msgbox(['Das Zuordnungchromosom der Individuen ',num2str(wrong_allo),' der population P wurde falsch erzeugt'],'Fehler')
end

%% INI: Kollaborative Operationen
if kolla{1,1} == 'j'
    for z1 = 1:p
        for s = 1:length(op_kolla)
            P_ini{z1,1}(2,op_kolla(1,s)) = 1;
            P_ini{z1,1}(3,op_kolla(1,s)) = 1;
        end
    end
end


%% INI: DECODIERUNG P_ini
% Wählbar: Semi-aktive und hybrides scheduling

% Semi-aktive:
if dec_typ == 1
    P_ini = decodierung_semi_aktiv_opt (P_ini,n,T,A_ges,Prae,op_kolla);
end

% hybrid:
if dec_typ == 2    
    P_ini = decodierung_opt (P_ini,n,T,A_ges,Prae,op_kolla);
end


%% INI: FITNESSFUNKION P_ini
if t_takt == 0
    P_ini = fitness(P_ini,n,T,E,E2,E3,E4,E5,WZ);
end

if t_takt > 0
    P_ini = fitness_takt(P_ini,n,T,E,E2,E3,E4,E5,WZ);
end


%% INI: Plotten der Fitness von P_ini
P_plot = sortrows(P_ini,-2);        % P soll nicht sortiert werden

for z1 = 1:p
    fit_plot(z1,g) = P_plot{z1,2}(1,1);     % Spalte = Generation, Zeile = Individuum
end


%% INI: Ranking P_ini
P_ini = ranking(P_ini,p,anz_raenge,R_max);


%% INI: Plotten der Ränge von P_ini
% Für das Plotten der Ränge und der Anz der Indi pro Rang
R = zeros(z1,1);
for z1=1:p
     R(z1,1) = P_ini{z1,5}(1,1);
end

for rang = 1:anz_raenge
    R_plot(rang,g) = sum(R==rang);
end


%% INI: Crowding - Distance P_ini
P_ini = crowding_distance(P_ini,p);


%% INI: Plotten der Crowding - Distance
for z1 = 1:p
    HA(z1,g) = mean(P_ini{z1,7}(:,1));  % Hamming Abstand
    CD(z1,g) = mean(P_ini{z1,7}(:,2));  % Delta Zeit
end


%% INI: Fitness Sharing von P_ini
% Ermittlung der Nachbarn jedes Individuums und teilen der Fitness 
P_ini = fitness_sharing(P_ini,p,hamming_abstand,delta_startzeit);


%% INI: Plotten der Eignung
% Bisher wird nur Eignung 1 geplottet (2 bis 5 sind ein Bonus und müssten
% manuelle geändert oder hinzugefügt werden

eig_plot(1,g) = P_plot{1,3}(5,1);   % Eignung Mensch
eig_plot(2,g) = P_plot{1,3}(5,2);   % Eignung Roboter

ausl_plot(1,g) = P_plot{1,3}(2,1);  % Auslastung Mensch
ausl_plot(2,g) = P_plot{1,3}(2,2);  % Auslastung Roboter

for z1 = 1:p
    ges_eig_plot(z1,g) = P_plot{z1,3}(5,1) * e(1,1) + P_plot{z1,3}(5,2) * e(2,1);     % Spalte = Generation, Zeile = Individuum
    ges_ausl_plot(z1,g) = P_plot{z1,3}(2,1) * a(1,1) + P_plot{z1,3}(2,2) * a(2,1);
end


%% INI: Gantt bestes Indi (nur erste Generation)
% ermittlung des besten Individuums
fitnessChampIni = 0.01;
for z1 = 1:p
    fit = P_ini{z1,2}(1,1);
    if fit > fitnessChampIni
        fitnessChampIni = fit;             % Fitnesswert des besten
        bestesIndiIni = z1;            % Zeile des besten
    end
end

champ_plot(1,g) = fitnessChampIni;      % Fitness des Champ in PLot eintragen


% Erstellung des Gantt-Diagramms
[Duration,Positions,mensch_sort,roboter_sort,anz_mensch,anz_roboter] = gantt(P_ini,bestesIndiIni);

if length(mensch_sort)> length(roboter_sort)
    anz_gaps = length(mensch_sort);   
elseif length(roboter_sort)>= length(mensch_sort)
    anz_gaps = length(roboter_sort);  
end

width = 0.2;
GAPS1 = zeros(1,anz_gaps); 
fig1 = figure('Name','Bestes Individuum der Initialisierungspopulation');
G1 = barh(Positions,Duration,width,'stacked','w');
hold on

set(gca,'yticklabel',{'Mensch','Roboter'});
for ind = 1:anz_gaps
    GAPS1(1,ind) = ind*2-1;    % Spalten bestimmen, die Gabs sind (also keine Prozesse)
end
set(G1(GAPS1),'Visible','off')    % Gabs nicht sichtbar machen

% Prozessnummer in bars eintragen
for ind = 1:anz_mensch 
        Pos1m = mensch_sort(5,ind)+(mensch_sort(6,ind)/2);   % Startzeit + Dauer/2
        text(Pos1m-0.3,1,num2str(mensch_sort(1,ind)))          
end

for ind = 1:anz_roboter
        Pos1r = roboter_sort(5,ind)+(roboter_sort(6,ind)/2);
        text(Pos1r-0.3,2,num2str(roboter_sort(1,ind)))
end
set(fig1,'Color',[1 1 1])
eig_mensch = ['Eignung 1: ',num2str(P_ini{bestesIndiIni,3}(5,1))];
ausl_mensch = ['Fitness Auslastung: ',num2str(P_ini{bestesIndiIni,3}(2,1))];

eig_roboter = ['Eignung 1: ',num2str(P_ini{bestesIndiIni,3}(5,2))];
ausl_roboter = ['Fitness Auslastung: ',num2str(P_ini{bestesIndiIni,3}(2,2))];

dim = [0.15 0.18 0.25 0.08];
str = {eig_mensch,ausl_mensch};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

dim = [0.15 0.73 0.25 0.08];
str = {eig_roboter,ausl_roboter};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
hold off


%% BEGINN DER GENERATIONEN SCHLEIFE
P_neu = P_ini;

for z1 = 1:p
    P_neu{z1,10}(1,1) = 1;
end

tic

%%
while g < g_max
    
    if g > 1
        clearvars Eltern Paare K PK
    end
    
    %% Übertragen der Individuen von P_neu zu P
    
    P(:,1) = P_neu(:,1);        % Chromosome
    P(:,2) = P_neu(:,2);        % Gesamtfitness
    P(:,3) = P_neu(:,3);        % Teilfitness
    P(:,4) = P_neu(:,4);        % Rangbasierte Fitness
    P(:,5) = P_neu(:,5);        % Rang übergeben
    P(:,6) = P_neu(:,6);        % Shared Fitness
    P(:,7) = P_neu(:,7);        % Crowding Distance
    P(:,8) = P_neu(:,8);        % Anzahl der Nachbarn
    % in Spalte 9 werden die Siege eingetragen
    P(:,10) = P_neu(:,10);
    
    clearvars P_neu             % P_neu löschen
    
    
    %% ELTERNMUTATION
    
    P_mut(:,1) = P(:,1);
    
    for z1 = 1:p
        P_mut{z1,1}(4,:)= zeros(1,n);
        P_mut{z1,1}(5,:)= zeros(1,n);
        P_mut{z1,1}(6,:)= zeros(1,n);
    end
    
    if anteil_mut_eltern ~= 0   
    % folgender Abschnitt wird nur benötigt, wenn Eltern zum Teil mutiert
    % werden sollen, dann ist anteil_mut_eltern größer als 0
    
    
        %% ELTERNMUTATION: Sequenzchromosom

        % evtl. als Funktion
        % erster gewählter Ansatz: Vertauschende Mutation
        if el_mutation_seq == 1
            [P_mut,MonitoringPermuMutation] = swap_mutation(P_mut,el_mut_whk_swap,p,n);
        end


        % zweiter gewählter Ansatz: Mischende Mutation
        if el_mutation_seq == 2
            genom = 1;                  % Sequenzgenom
            el_mut_whk = 1.0;           % Jedes Elternindividuum soll mutiert werden
            P_mut = mixed_mutation(P_mut,p,n,genom,el_mut_whk);
        end
        
        
        % dritter Ansatz: Invertierende Mutation
        if el_mutation_seq == 3
            genom = 1;                  % Sequenzgenom
            el_mut_whk = 1.0;           % Jedes Elternindividuum soll mutiert werden
            P_mut = invert_mutation(P_mut,p,n,genom,el_mut_whk);
        end


        %% ELTERNMUTATION: Zuordnungschromosom

        % erster gewählter Ansatz: 1-Bit-Mutation
        if el_mutation_allo == 1
            [P_mut,MonitoringAlloMutation] = ein_bit_mutation (P_mut,el_mut_whk_einbit,p,n);
        end


        % zweiter gewählter Ansatz: Invertierende Mutation
        if el_mutation_allo == 2
            genom = 2;                  % Sequenzgenom
            el_mut_whk = 1.0;           % Jedes Elternindividuum soll mutiert werden
            P_mut = invert_mutation(P_mut,p,n,genom,el_mut_whk);
        end


        % möglicher dritter Ansatz: Mischende Mutation
        if el_mutation_allo == 3
            genom = 2;                  % Sequenzgenom
            el_mut_whk = 1.0;           % Jedes Elternindividuum soll mutiert werden
            P_mut = mixed_mutation(P_mut,p,n,genom,el_mut_whk);
        end

        for z1 = 1:p            
            P_mut{z1,1}(3,:) = ones(1,n) - P_mut{z1,1}(2,:);    % Roboter aus Mensch Zuordnung berechnen
        end


        %% ELTERNMUTATION: Korrekturzuordnung
        % Sequenz- und Zuordnungschromosom korrigieren
        if correct{1,1} == 'j' 
            P_mut = manuelle_allocation(P_mut,p,geneM,geneR,permutation);
        end


        %% ELTERNMUTATION: Plausiblitätstest

        [wrong_seq,wrong_allo] = plausibilitaet (P_mut,p,n);

        if ~isempty(wrong_seq)
            msgbox(['Das Sequenzchromosom der Individuen ',num2str(wrong_seq),' der Population K wurde falsch erzeugt'],'Fehler')
        end

        if ~isempty(wrong_allo)
            msgbox(['Das Zuordnungchromosom der Individuen ',num2str(wrong_allo),' der population K wurde falsch erzeugt'],'Fehler')
        end
        
        
        %% ELTERNMUTATION: Kollaborative Operationen
        if kolla{1,1} == 'j'
            for z1 = 1:p
                for s = 1:length(op_kolla)
                    P_mut{z1,1}(2,op_kolla(1,s)) = 1;
                    P_mut{z1,1}(3,op_kolla(1,s)) = 1;
                end
            end
        end


        %% ELTERNMUTATION: Decodierung
        % Wählbar: Semi-aktive und hybrides scheduling

        % Semi-aktive:
        if dec_typ == 1
            P_mut = decodierung_semi_aktiv_opt (P_mut,n,T,A_ges,Prae,op_kolla);
        end

        % hybrid:
        if dec_typ == 2    
            P_mut = decodierung_opt (P_mut,n,T,A_ges,Prae,op_kolla);
        end
        
        
        %% ELTERNMUTATION: Fitness  
        if t_takt == 0
            [P_mut] = fitness(P_mut,n,T,E,E2,E3,E4,E5,WZ);
        end

        if t_takt > 0
            [P_mut] = fitness_takt(P_mut,n,T,E,E2,E3,E4,E5,WZ);
        end


        %% ELTERNMUTATION: ENDE
        P_mut = sortrows(P_mut,-2);     % Sortieren nach Fitness
        
        for z1 = 1:p
            P_mut{z1,5}(1,1) = 1;       % Generation der mutierten Eltern auf 1 setzen
        end
    end
    
    %% ELTERNSELEKTION
    % Wählbar: Fitness, anhand der selektiert wird
    % Wählbar: Verfahren
    
  
    %% ELT-SELEK: Selektionsverfahren
    
    % Erster gewählter Ansatz: Stochastic Universal Sampling
    
    % ERKLÄRUNG:
    % zunächst muss die Fitness gewählt werden (in der config) auf dessen Basis die
    % "Tortensegmente" berechnet werden:
    
    % Tatsächliche Fitness des Individuums: fitness_typ = 2
    % Rangbasierte Fitness des Individuums: fitness_typ = 4
    % Sharing Fitness des Individuums: fitness_typ = 6
    
    if eltern_selektion == 1
        [Eltern,P,Pointer,SUS_Fitness] = SUS_elternselektion (P,p,fitness_typ);
    end
    
    
    % Zweiter gewählter Ansatz: 
    
    % if eltern_selektion == 2
    %    
    % end
    
    
    %% CROSSOVER
    
    %% CROSS: Initialisierung

    % Bildung von p Paaren
    Paare = cell(p,2);
    ZuordnungEltern(:,1) = randperm(2*p);

    for z1 = 1:p
        partner1 = find(ZuordnungEltern(:,1)==2*z1-1);
        partner2 = find(ZuordnungEltern(:,1)==2*z1);
        Paare{z1,1} = Eltern{partner1,1}(1:3,:); 
        Paare{z1,2} = Eltern{partner2,1}(1:3,:);
    end

    K = cell(p,1);


    %% CROSS: Sequenzchromosom

    % gewählter Ansatz: Ordered Crossover (Ordnungsrekombination)
    if crossover_seq == 1
        K = ordered_crossover(K,Paare,p,n,PermuCrossWhk);
    end

    % % gewählter Ansatz: Partially Mapped Crossover
    if crossover_seq == 2
        K = pmx_crossover(K,Paare,p,n,PermuCrossWhk);
    end
    
        
    %% CROSS:Zuordnungschromosom
    
    % gewählter Ansatz: 2-Punkt-Crossover
    if crossover_allo == 1
        K = two_point_crossover(K,Paare,p,n,AlloCrossWhk);
    end
    
    % zweiter Ansatz: Uniformer-Crossover
    if crossover_allo == 2
        K = uniform_crossover(K,Paare,p,n,co_unif_whk,AlloCrossWhk);
    end
    
    for z1 = 1:p            
        K{z1,1}(3,:) = ones(1,n) - K{z1,1}(2,:);    % Roboter aus Mensch Zuordnung berechnen
    end
    
    
    %% CROSS: Plausiblitätstest
    
    [wrong_seq,wrong_allo] = plausibilitaet (K,p,n);
    
    if ~isempty(wrong_seq)
        msgbox(['Das Sequenzchromosom der Individuen ',num2str(wrong_seq),' der Population K wurde falsch erzeugt'],'Fehler')
    end
    
    if ~isempty(wrong_allo)
        msgbox(['Das Zuordnungchromosom der Individuen ',num2str(wrong_allo),' der population K wurde falsch erzeugt'],'Fehler')
    end
    
        
    %% MUTATION

    %% MUTA: Sequenzchromosom

    % erster gewählter Ansatz: Vertauschende Mutation
    if mutation_seq == 1
        [K,MonitoringPermuMutation] = swap_mutation(K,mut_whk_swap,p,n);
    end
    
    % zweiter gewählter Ansatz: Mischende Mutation
    if mutation_seq == 2
        genom = 1;                  % Sequenzgenom
        K = mixed_mutation(K,p,n,genom,PermuMutationWhk);
    end
    
    % dritter Ansatz: Invertierende Mutation
    if mutation_seq == 3
        genom = 1;                  % Sequenzgenom
        K = invert_mutation(K,p,n,genom,PermuMutationWhk);
    end
    
    
    %% MUTA: Zuordnungschromosom
    
    % erster gewählter Ansatz: 1-Bit-Mutation
    if mutation_allo == 1
        [K,MonitoringAlloMutation] = ein_bit_mutation (K,mut_whk_einbit,p,n);
    end
    
    
    % zweiter gewählter Ansatz: Invertierende Mutation
    if mutation_allo == 2
        genom = 2;                  % Sequenzgenom
        K = invert_mutation(K,p,n,genom,AlloMutationWhk);
    end
    
    
    % möglicher dritter Ansatz: Mischende Mutation
    if mutation_allo == 3
        genom = 2;                  % Sequenzgenom
        K = mixed_mutation(K,p,n,genom,AlloMutationWhk);
    end
    
    for z1 = 1:p            
        K{z1,1}(3,:) = ones(1,n) - K{z1,1}(2,:);    % Roboter aus Mensch Zuordnung berechnen
    end
    
    
    %% KINDER: Korrekturzuordnung
    % Sequenz- und Zuordnungschromosom korrigieren
    if correct{1,1} == 'j' 
        K = manuelle_allocation(K,p,geneM,geneR,permutation);
    end
    
    %% KINDER: Plausiblitätstest

    [wrong_seq,wrong_allo] = plausibilitaet (K,p,n);

    if ~isempty(wrong_seq)
        msgbox(['Das Sequenzchromosom der Individuen ',num2str(wrong_seq),' der Population K wurde falsch erzeugt'],'Fehler')
    end

    if ~isempty(wrong_allo)
        msgbox(['Das Zuordnungchromosom der Individuen ',num2str(wrong_allo),' der population K wurde falsch erzeugt'],'Fehler')
    end
    
    
    %% KINDER: Kollaborative Operationen
    if kolla{1,1} == 'j'
        for z1 = 1:p
            for s = 1:length(op_kolla)
                K{z1,1}(2,op_kolla(1,s)) = 1;
                K{z1,1}(3,op_kolla(1,s)) = 1;
            end
        end
    end
    
    
    %% KINDER: Decodierung
    % Wählbar: Semi-aktive und hybrides scheduling
    
    % Semi-aktive:
    if dec_typ == 1
        K = decodierung_semi_aktiv_opt (K,n,T,A_ges,Prae,op_kolla);
    end
    
    % hybrid:
    if dec_typ == 2    
        K = decodierung_opt (K,n,T,A_ges,Prae,op_kolla);
    end
    
    
    %% KINDER: Fitness  
    if t_takt == 0
        K = fitness(K,n,T,E,E2,E3,E4,E5,WZ);
    end
    
    if t_takt > 0
        K = fitness_takt(K,n,T,E,E2,E3,E4,E5,WZ);
    end
    
    
    %% UMWELTSELEKTION
    
    for z1 = 1:p
        K{z1,10}(1,1) = 0;     % Alle Kinder erhalten Generation 0
    end
    
    % Schlechte Eltern durch beste mutierte Eltern ersetzen
    P = sortrows(P,-2);                 % P ist immer noch sortiert, aber nach Anzahl der Siege, daher neu sortieren nach Fitness
                                        % und mutierte Eltern müssen ebenfalss berücksichtigt werden
    % Werte von P löschen, die nicht mehr gebraucht werden
    P(:,4) = [];            % Rangbasierte Fitness 
    P(:,4) = [];            % Rang übergeben
    P(:,4) = [];            % Shared Fitness
    P(:,4) = [];            % Crowding Distance
    P(:,4) = [];            % Anzahl der Nachbarn
    
    P((p-anz_mut_eltern)+1:end,:) = P_mut(1:anz_mut_eltern,:);    % beste neu mutierte Eltern hinzufügen
    
    % Plus - Selektion: Eltern und Kinder (alt)
%     PK(1:p,1) = P(1:p,1);           % p Eltern in PK eintragen 
%     PK(1:p,2) = P(1:p,2);           % Gesamtfitness hinzufügen
%     PK(1:p,3) = P(1:p,3);           % Fitnesswerte hinzufügen

%     PK(p+1:2*p,1) = K(1:p,1);       % p Kinder in PK eintragen
%     PK(p+1:2*p,2) = K(1:p,2);       % Gesamtfitness hinzufügen
%     PK(p+1:2*p,3) = K(1:p,3);       % Fitnesswerte hinzufügen

    % Neu
    if p_el ~= 0        % wenn p_el = 0, dann PK nur aus Kindern
        PK(1:p_el,1) = P(1:p_el,1);             % p_el Eltern in PK eintragen 
        PK(1:p_el,2) = P(1:p_el,2);             % Gesamtfitness hinzufügen
        PK(1:p_el,3) = P(1:p_el,3);             % Fitnesswerte hinzufügen
        PK(1:p_el,10) = P(1:p_el,5);
    end
    
    PK(p_el+1:(p_el+p_ki),1) = K(1:p_ki,1);     % p_ki Kinder in PK eintragen
    PK(p_el+1:(p_el+p_ki),2) = K(1:p_ki,2);     % Gesamtfitness hinzufügen
    PK(p_el+1:(p_el+p_ki),3) = K(1:p_ki,3);     % Fitnesswerte hinzufügen
    PK(p_el+1:(p_el+p_ki),10) = K(1:p_ki,10); 
    
    p_sel = length(PK);                         % Populationsgröße für die Umweltselektion
    
    
    %% UMW-SELEK: Ranking PK
    PK = ranking(PK,p_sel,anz_raenge,R_max);
    
    
    %% UMW-SELEK: Crowding Distance PK
    PK = crowding_distance(PK,p_sel);
    
    
    %% UMW-SELEK: Fitness Sharing PK
    PK = fitness_sharing(PK,p_sel,hamming_abstand,delta_startzeit);
    
    
    %% UMW-SELEK: Champion bestimmen
    % bei einem Champion besteht die Gefahr eines Superindividuums
    fitnessChamp = 0;
    for z1 = 1:p_sel
        fit = PK{z1,2}(1,1);
        if fit > fitnessChamp
            fitnessChamp = fit;         % Fitnesswert des besten
            champ = z1;                 % Zeile des Besten
        end
    end
    
    
    %% UMW-SELEK: Selektionsverfahren
    
    % Erster gewählter Ansatz: Q-Stufige-Tunierselektion
    if umwelt_selektion == 1
        PK = q_tunier_selektion (PK,p_sel,q,champ,umwelt_fitness,champion_confi);
    end
    
    P_neu = sortrows(PK,[-9 -2]);       % Sortieren nach Anzahl der Siege und Fitness und speichern in P_neu
    
    
    
    % Zweiter gewählter Ansatz: 
    % if umwelt_selektion == 2
    %     PK = q_tunier_selektion (PK,p,q);
    % end
    
    
    %% UMW-SELEK: Neue Generation (g+1) erstellen
  
    P_neu(p+1:end,:) = [];      % löschen der schlechten Individuen
                                % GENERATION G+1 ERSTELLT
    
    for z1 = 1:p
        P_neu{z1,10}(1,1) = P_neu{z1,10}(1,1) + 1; 
    end                          
                             
    clearvars MonitoringPermuMutation MonitoringAlloMutation 
    clearvars partner1 partner2 gegner1 gegner2 ZuordnungEltern
    clearvars wrong_allo wrong_seq
    clearvars ind index fit
       
    
    %% NÄCHSTE GENERATION
    g = g + 1               % Ausgabe offen, um Fortschritt anzeigen zu lassen
    
    
    %% NEXY GEN: Plotten des Champion
    champ_plot(1,g) = fitnessChamp;
    
    
    %% NEXT GEN: Plotten der Fitness von P_neu
    P_plot = sortrows(P_neu,-2);
    
    for z1 = 1:p
        fit_plot(z1,g) = P_plot{z1,2}(1,1);     % Spalte = Generation, Zeile = Individuum
    end
    
    
    %% NEXT GEN: Ranking P
    P_neu = ranking(P_neu,p,anz_raenge,R_max);
    
    
    %% NEXT GEN: Plotten der Ränge 
    % Für das Plotten der Ränge und der Anz der Indi pro Rang
    R = zeros(z1,1);
    for z1=1:p
         R(z1,1) = P_neu{z1,5}(1,1);
    end
    
    for rang = 1:anz_raenge
        R_plot(rang,g) = sum(R==rang);
    end
    
    
    %% NEXT GEN: Crowding - Distance P
    P_neu = crowding_distance(P_neu,p);
    
    
    %% NEXT GEN: Plotten der Crowding - Distance
    for z1 = 1:p
        HA(z1,g) = mean(P_neu{z1,7}(:,1));
        CD(z1,g) = mean(P_neu{z1,7}(:,2));
    end
    
    
    %% NEXT GEN: Fitness Sharing von P
    % Ermittlung der Nachbarn jedes Individuums und teilen der Fitness 
    P_neu = fitness_sharing(P_neu,p,hamming_abstand,delta_startzeit);
    
    
    %% NEXT GEN: Plotten der Eignung und Auslastung
    % P-plot ist nach Fintess soritert, das erste Indi ist das beste
    eig_plot(1,g) = P_plot{1,3}(5,1);   % Eignung Mensch
    eig_plot(2,g) = P_plot{1,3}(5,2);   % Eignung Roboter
    
    ausl_plot(1,g) = P_plot{1,3}(2,1);  % Auslastung Mensch
    ausl_plot(2,g) = P_plot{1,3}(2,2);  % Auslastung Roboter
    
    for z1 = 1:p
        ges_eig_plot(z1,g) = P_plot{z1,3}(5,1) * e(1,1) + P_plot{z1,3}(5,2) * e(2,1);     % Spalte = Generation, Zeile = Individuum
        ges_ausl_plot(z1,g) = P_plot{z1,3}(2,1) * a(1,1) + P_plot{z1,3}(2,2) * a(2,1);
    end
    
    
    %% ABBRUCHKRITERIUM
    
    % Abbruchkriterium: Kein Verbesserung des Champion
    if abbruch_typ1 == 1
        abbruch_champ = champ_plot; % Umspeichern
        abbruch_champ(abbruch_champ == 0) = []; % Nullen löschen
        sz_abbruch_champ = length(abbruch_champ);

        if sz_abbruch_champ > abbruch_gen      % Erst anfangen wenn i = abbruch_gen Champion vorhanden
            abbruch_diff = abbruch_champ(end) - abbruch_champ(end-abbruch_gen); % Fitnessdifferenz berechnen
            if abbruch_diff < abbruch_wert  % Wenn Diff kleiner ist als min. zul., dann
                g = g_max;                  % g auf max. Generationen setzen
            end
        end
    end
    
    
    % Abbruchkriterium: Zeit
    if abbruch_typ2 == 1
        time = toc;
        if time > time_max
            g = g_max;
        end
    end
    
    %% ENDE DER GENERATIONENSCHLEIFE
    
    
end

%% AUSGABE: finale Population
P_final = P_neu;


%% AUSGABE: GANTT bestes Individuum

bestIndi = 0.01;
for z1 = 1:p
    fit = P_final{z1,2}(1,1);
    if fit > bestIndi
        bestIndi = fit;             % Fitnesswert des besten
        bestesIndi = z1;            % Zeile des besten
    end
end

champ_plot(1,g) = bestIndi;

[Duration,Positions,mensch_sort,roboter_sort,anz_mensch,anz_roboter] = gantt(P_final,bestesIndi);

if length(mensch_sort)> length(roboter_sort)
    anz_gaps = length(mensch_sort);
elseif length(roboter_sort)>= length(mensch_sort)
    anz_gaps = length(roboter_sort);
end

width = 0.2;
GAPS2 = zeros(1,anz_gaps);

fig2 = figure('Name','Bestes Individuum der letzten Generation');
G2 = barh(Positions,Duration,width,'stacked','w');
hold on

set(gca,'yticklabel',{'Mensch','Roboter'});
for ind = 1:anz_gaps
    GAPS2(1,ind) = ind*2-1;
end
set(G2(GAPS2),'Visible','off')


for ind = 1:anz_mensch
        Pos2m = mensch_sort(5,ind)+ (mensch_sort(6,ind)/2);
        text(Pos2m,1,num2str(mensch_sort(1,ind)))
end


for ind = 1:anz_roboter
        Pos2r = roboter_sort(5,ind)+(roboter_sort(6,ind)/2);
        text(Pos2r,2,num2str(roboter_sort(1,ind)))
end
set(fig2,'Color',[1 1 1])
eig_mensch = ['Eignung 1: ',num2str(P_final{bestesIndi,3}(5,1))];
ausl_mensch = ['Fitness Auslastung: ',num2str(P_final{bestesIndi,3}(2,1))];

eig_roboter = ['Eignung 1: ',num2str(P_final{bestesIndi,3}(5,2))];
ausl_roboter = ['Fitness Auslastung: ',num2str(P_final{bestesIndi,3}(2,2))];

dim = [0.15 0.18 0.25 0.08];
str = {eig_mensch,ausl_mensch};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

dim = [0.15 0.73 0.25 0.08];
str = {eig_roboter,ausl_roboter};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

hold off


%% AUSGABE: Entwicklung des GA

figure('Name','Fitness')
bar3(fit_plot)

figure('Name','Ränge')
bar3(R_plot)

figure('Name','Hamming-Abstand')
bar3(HA)

figure('Name','Delta-Zeit')
bar3(CD)

figure('Name','Fitness des Champion jeder Population')
plot(champ_plot)

figure('Name','Gesamteignung jedes Individuums')
bar3(ges_eig_plot)

figure('Name','Gesamtauslastung jedes Individuums')
bar3(ges_ausl_plot)

eig_plot = eig_plot.';
figure('Name','Eignung des Champion')
bar(eig_plot)
lgd2 = legend({'Mensch','Roboter'},'Location','southoutside');

ausl_plot = ausl_plot.';
figure('Name','Auslastung des Champion')
bar(ausl_plot)
lgd3 = legend({'Mensch','Roboter'},'Location','southoutside');

%% SPEICHERN WORKSPACE

filename = [prozessname '_workspace.mat'];
save(filename);
