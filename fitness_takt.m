function [Population] = fitness_takt(Population,n,T,E,E2,E3,E4,E5,WZ) 

%% Initialisierung
run('config')
run('trans')

% Erstellung der Fitnessmatrix für jedes Individuum  
for z1 = 1:p
    Population{z1,2} = zeros(1,1);
    Population{z1,3} = zeros(9,2);
end


%% Pausen und Auslastung: Vorbereitung
pausen_mensch = cell(p,1);
pausen_roboter = cell(p,1);

% Prozesse bestimmen notwendig um Pausen zu ermitteln: duration_
for z1=1:p
    prozesse_m = find(Population{z1,1}(2,:)==1);                % Suche Prozesse auf dem Mensch
    mensch = [prozesse_m ; Population{z1,1}(1:6,prozesse_m)];   % Prozesse des Menschen abspeichern
    [~, index] = sort(mensch',1);                               % 
    mensch_sort = mensch(:,index(:,5));                         % Prozesse sortieren nach der Zeit
 
    prozesse_r = find(Population{z1,1}(3,:)==1);                % Suche Prozesse auf dem Roboter
    roboter = [prozesse_r ; Population{z1,1}(1:6,prozesse_r)];  % Prozesse des Roboters abspeichern
    [~, index] = sort(roboter',1);                              %
    roboter_sort = roboter(:,index(:,5));                       % Prozesse sortieren nach der Zeit
    
    [~,anz_mensch] = size(mensch_sort);                         % Anzahl der Prozesse Mensch
    [~,anz_roboter] = size(roboter_sort);                       % Anzahl der Prozesse Roboter
    
    
    duration_mensch = zeros(1,anz_mensch*2);
    duration_roboter = zeros(1,anz_roboter*2);
    
    
    for s = 1:anz_mensch
        if s == 1
            duration_mensch(1,1) = mensch_sort(5,1);        % Wartezeit beim Start berechnen und in Spalte 1
            duration_mensch(1,2) = mensch_sort(7,1) - mensch_sort(5,1); % Prozessdauer von Prozess 1 in Spalte 2
        else
            duration_mensch(1,(2*s-1)) = mensch_sort(5,s) - mensch_sort(7,(s-1));   % Pause zwischen s (ab 2 gezählt) und s-1 berechnen
            duration_mensch(1,(2*s)) = mensch_sort(7,s) - mensch_sort(5,s);         % Prozesszeit von s
        end
    end
    
    for s = 1:anz_roboter
        if s == 1
            duration_roboter(1,1) = roboter_sort(5,1);      % Wartezeit beim Start berechnen
            duration_roboter(1,2) = roboter_sort(7,1) - roboter_sort(5,1);
        else
            duration_roboter(1,(2*s-1)) = roboter_sort(5,s) - roboter_sort(7,(s-1));
            duration_roboter(1,(2*s)) = roboter_sort(7,s) - roboter_sort(5,s);
        end
    end
    
    % Matrizen für Pausen erstellen
    pausen_mensch{z1,1} = zeros(1,anz_mensch);      % Die Wartezeit am Ende wird später hinzugefügt
    pausen_roboter{z1,1} = zeros(1,anz_roboter);
    
    
    for s = 1:anz_mensch
        pausen_mensch{z1,1}(1,s) = duration_mensch(1,s*2-1);  % Ungerade Stellen sind Pausen
    end
    
    
    for s = 1:anz_roboter
        pausen_roboter{z1,1}(1,s) = duration_roboter(1,s*2-1);    % Ungerade Stellen sind Pausen
    end
end


%% Anzahl der Pausen
% Mensch: möglist eine Pausenzeit
% Roboter: möglichst kurze kleine Pausenzeiten
% Belohnungssystem

if w(3,1) == 1
    anz_pausen = zeros(p,2);

    for z1=1:p
        anz_pausen(z1,1) = length(find(pausen_mensch{z1,1}(1,:)>0));   % Anz Pausen beim Mensch
        anz_pausen(z1,2) = length(find(pausen_roboter{z1,1}(1,:)>0));  % Anz Pausen beim Roboter
    end

    % Wartezeit zum Schluss zählt auch als Pause, muss dazugezählt werden
    for z1 = 1:p
        zeiten_mensch = Population{z1,1}(2,:) .* Population{z1,1}(6,:);     % Alle Prozesse des Menschen 
        zeiten_roboter = Population{z1,1}(3,:) .* Population{z1,1}(6,:);    % Alle Prozesse des Roboters
        belegung_mensch = max(zeiten_mensch);                               % Maximum Mensch
        belegung_roboter = max(zeiten_roboter);                             % Maximum Roboter


        if belegung_mensch > belegung_roboter
            anz_pausen(z1,2) = anz_pausen(z1,2) + 1;  % Roboter wartet
        elseif belegung_mensch < belegung_roboter
            anz_pausen(z1,1) = anz_pausen(z1,1) + 1;  % Mensch wartet
        end
    end


    % Mensch
    for z1 = 1:p
        pausen_mensch_max = max(anz_pausen(:,1));   % Maximum der Pausen aller Indi in dieser Generation
        pausen_mensch_min = min(anz_pausen(:,1));   % Minimum der Pausen allen Indi in dieser Generation
        pausen_mensch_range = max_belohn_pausen / (pausen_mensch_max - pausen_mensch_min);  % Range der Belohnung bestimmen

        % Belohnung = Max_Belohnung - (Ist_Pausen - Min_Pausen)*Range
        Population{z1,3}(3,1) = max_belohn_pausen - ( (anz_pausen(z1,1)-pausen_mensch_min) * pausen_mensch_range);  % Belohnung ermitteln
    end


    % Roboter
    for z1 = 1:p
        pausen_roboter_max = max(anz_pausen(:,2));
        pausen_roboter_min = 1;
        pausen_roboter_range = max_belohn_pausen / (pausen_roboter_max - pausen_roboter_min);

        if anz_pausen(z1,2) == 0
            Population{z1,3}(3,2) = max_belohn_pausen;
        else
            Population{z1,3}(3,2) = (anz_pausen(z1,2) - 1)*pausen_roboter_range;
        end
    end
end


% Alter ansatz
% for z1=1:p
%     % Anpassung Pausen Mensch
%     Population{z1,3}(3,1) = 1-0.02*(anz_pausen{z1,1}(1,1)-1);
%     
%     % Anpassung Pausen Roboter
%     anz_roboter = length(find(Population{z1,1}(3,:)==1));   % Ermittlung der Prozesse auf R
%     anz_gabs = anz_roboter - 1;                             % Eine Pause weniger als Prozesse
%     if anz_pausen{z1,2}(1,1) > 0                            % bei min. einer Pause
%       Population{z1,3}(3,2) = 1-0.02*(anz_gabs - anz_pausen{z1,2}(1,1));
%     elseif anz_pausen{z1,2}(1,1) == 0                       % keine Pause, Optimum
%       Population{z1,3}(3,2) = 1;
%     end
% end



%% Auslastung 
% Fitnesswert zwischen 0 und 1
fitness_auslast = cell(p,2);

for z1=1:p
    
    warten_roboter = 0;
    warten_mensch = 0;
    
    taktzeit = max(Population{z1,1}(6,:));
    
    zeiten_mensch = Population{z1,1}(2,:) .* Population{z1,1}(6,:);     % Alle Prozesse des Menschen 
    zeiten_roboter = Population{z1,1}(3,:) .* Population{z1,1}(6,:);    % Alle Prozesse des Roboters
    belegung_mensch = max(zeiten_mensch);                               % Maximum Mensch
    belegung_roboter = max(zeiten_roboter);                             % Maximum Roboter
    
    % Berechnung der Wartezeiten, da diese nicht in pausen enthalten sind
    if belegung_mensch > belegung_roboter                       % Roboter zuerst fertig
        warten_roboter = belegung_mensch - belegung_roboter;    % Wartezeit des Roboters
    elseif belegung_mensch < belegung_roboter                   % Mensch zuerst fertig
        warten_mensch = belegung_roboter - belegung_mensch;     % Wartezeit des Menschen
    end
    
    % 1-(pausen + wartezeit / taktzeit)
    fitness_auslast{z1,1}(1,1) = 1-(sum(pausen_mensch{z1,1}(1,:))+warten_mensch)/taktzeit;
    fitness_auslast{z1,2}(1,1) = 1-(sum(pausen_roboter{z1,1}(1,:))+warten_roboter)/taktzeit;
end



for z1=1:p
    % Auslastung Mensch (Maximum bei 0.97)
    if fitness_auslast{z1,1}(1,1) <= 0.97
        Population{z1,3}(2,1) = 1/0.97*fitness_auslast{z1,1}(1,1);
    elseif fitness_auslast{z1,1}(1,1) > 0.97
        Population{z1,3}(2,1) = -(1/0.97)*fitness_auslast{z1,1}(1,1)+2;
    end
    
    % Auslastung Roboter (Maximum bei 1.00)
    Population{z1,3}(2,2) = fitness_auslast{z1,2}(1,1);
end


%% Werkzeuge
% Fitness als Belohnung
same_werkz = zeros(p,2);
WZ(isnan(WZ)) = 0;

if max(WZ)>0.5 && w(4,1) == 1  % Wenn WZgruppen eingetragen sind und die Gewichtung auf 1 steht

    for z1 = 1:p

        prozesse_m = find(Population{z1,1}(2,:)==1);                % Suche Prozesse auf dem Mensch
        mensch = [prozesse_m ; Population{z1,1}(1:6,prozesse_m)];   % Prozesse des Menschen abspeichern
        [~, index] = sort(mensch',1);                               % 
        mensch_sort = mensch(:,index(:,5));                         % Prozesse sortieren nach der Zeit

        prozesse_r = find(Population{z1,1}(3,:)==1);                % Suche Prozesse auf dem Roboter
        roboter = [prozesse_r ; Population{z1,1}(1:6,prozesse_r)];  % Prozesse des Roboters abspeichern
        [~, index] = sort(roboter',1);                              %
        roboter_sort = roboter(:,index(:,5));                       % Prozesse sortieren nach der Zeit


        % Mensch
        [~,sN]=size(mensch_sort);   % Sortierte Operationen auf der Ressource
        for s = 1:(sN-1)
            op1 = mensch_sort(1,s); 
            op2 = mensch_sort(1,s+1);
            if WZ(op1,1) == WZ(op2,1) && WZ(op1,1) ~= 0
                same_werkz(z1,1) = same_werkz(z1,1) + 1;    % Anzahl der gleichen WZ speichern
            end
        end


        % Roboter 
        [~,sN]=size(roboter_sort);  % Sortierte Operationen auf der Ressource
        for s = 1:(sN-1)
            op1 = roboter_sort(1,s);
            op2 = roboter_sort(1,s+1);
            if WZ(op1,1) == WZ(op2,1) && WZ(op1,1) ~= 0
                same_werkz(z1,2) = same_werkz(z1,2) + 1;    % Anzahl der gleichen WZ speichern
            end
        end

    end

    max_werk_m = max(same_werkz(:,1));  % Maximum der gleichen Werkzeuge bestimmen
    min_werk_m = min(same_werkz(:,1));  % Minimum der gleichen Werkzeuge bestimmen


    max_werk_r = max(same_werkz(:,2));  % Maximum der gleichen Werkzeuge bestimmen
    min_werk_r = min(same_werkz(:,2));  % Minimum der gleichen Werkzeuge bestimmen

    range_werk_m = max_belohn_werkzeuge / (max_werk_m - min_werk_m);    % Bereich zwischen Maximum und Minimum bestimmen
    range_werk_r = max_belohn_werkzeuge / (max_werk_r - min_werk_r);    % Bereich zwischen Maximum und Minimum bestimmen


    for z1 = 1:p
        Population{z1,3}(4,1) = same_werkz(z1,1) * range_werk_m;    % Belohnung eintragen
        Population{z1,3}(4,2) = same_werkz(z1,2) * range_werk_r;    % Belohnung eintragen
    end
    
else
    
    for z1 = 1:p
        Population{z1,3}(4,1) = 0;
        Population{z1,3}(4,2) = 0;
    end
    
end
    
  
%% Taktzeit
% Berechnung der normierten Taktzeit
    % Berechnung der maximale Zeit des Prozesses
    % Summe alle Prozesszeiten 1-n (wobei immer die Prozesszeit von
    % Mensch und Roboter gewählt wird, die länger dauert)
f_zeit_ob = 0;
for z2 = 1:n
    p_max = max(T(z2,:));
    f_zeit_ob = f_zeit_ob + p_max;
end

f_zeit_ub = 0;
for z2 = 1:n
    p_min = min(T(z2,:));
    f_zeit_ub = f_zeit_ub + p_min;
end

f_zeit_ub = f_zeit_ub/3;


% Fitnesswert jedes Individuums bzgl. der Taktzeit
for z1 = 1:p
    Population{z1,3}(1,1) = 1 - ((max(Population{z1,1}(6,:))-f_zeit_ub) / (f_zeit_ob-f_zeit_ub));
end



%% Eignungsgrad
% Berechung des Eignungsgrades
if anz_eig >= 1
    for z1 = 1:p
        E_Mensch = Population{z1,1}(2,:)*E(:,1) / sum(Population{z1,1}(2,:)==1);
        E_Roboter = Population{z1,1}(3,:)*E(:,2) / sum(Population{z1,1}(3,:)==1);
        Population{z1,3}(5,1) = E_Mensch;
        Population{z1,3}(5,2) = E_Roboter;
    end
end

if anz_eig >= 2
    for z1 = 1:p
        E_Mensch = Population{z1,1}(2,:)*E2(:,1) / sum(Population{z1,1}(2,:)==1);
        E_Roboter = Population{z1,1}(3,:)*E2(:,2) / sum(Population{z1,1}(3,:)==1);
        Population{z1,3}(6,1) = E_Mensch;
        Population{z1,3}(6,2) = E_Roboter;
    end
end

if anz_eig >= 3
    for z1 = 1:p
        E_Mensch = Population{z1,1}(2,:)*E3(:,1) / sum(Population{z1,1}(2,:)==1);
        E_Roboter = Population{z1,1}(3,:)*E3(:,2) / sum(Population{z1,1}(3,:)==1);
        Population{z1,3}(7,1) = E_Mensch;
        Population{z1,3}(7,2) = E_Roboter;
    end
end

if anz_eig >= 4
    for z1 = 1:p
        E_Mensch = Population{z1,1}(2,:)*E4(:,1) / sum(Population{z1,1}(2,:)==1);
        E_Roboter = Population{z1,1}(3,:)*E4(:,2) / sum(Population{z1,1}(3,:)==1);
        Population{z1,3}(8,1) = E_Mensch;
        Population{z1,3}(8,2) = E_Roboter;
    end
end

if anz_eig >= 5
    for z1 = 1:p
        E_Mensch = Population{z1,1}(2,:)*E5(:,1) / sum(Population{z1,1}(2,:)==1);
        E_Roboter = Population{z1,1}(3,:)*E5(:,2) / sum(Population{z1,1}(3,:)==1);
        Population{z1,3}(9,1) = E_Mensch;
        Population{z1,3}(9,2) = E_Roboter;
    end
end


%% Fitnessfunktion

for z1= 1:p
    
    % Verteilung der Einzelwerte berechnen
    Taktzeit = Population{z1,3}(1,1);
    Auslastung = Population{z1,3}(2,1) * a(1,1) + Population{z1,3}(2,2) * a(2,1);
    Pausen = Population{z1,3}(3,1) * b(1,1) + Population{z1,3}(3,2) * b(2,1);       % b: breaks
    Werkzeuge = Population{z1,3}(4,1) * t(1,1) + Population{z1,3}(4,2) * t(2,1);    % t: tools
    
    Eignung = Population{z1,3}(5,1) * e(1,1) + Population{z1,3}(5,2) * e(2,1);
    Eignung2 = Population{z1,3}(6,1) * e2(1,1) + Population{z1,3}(6,2) * e2(2,1);
    Eignung3 = Population{z1,3}(7,1) * e3(1,1) + Population{z1,3}(7,2) * e3(2,1);
    Eignung4 = Population{z1,3}(8,1) * e4(1,1) + Population{z1,3}(8,2) * e4(2,1);
    Eignung5 = Population{z1,3}(9,1) * e5(1,1) + Population{z1,3}(9,2) * e5(2,1);
    
    
    % Gesamtfitnessfunktion
    Population{z1,2}(1,1) =  Auslastung * w(2,1) + Eignung * w(5,1) + Eignung2 * w(6,1) + Eignung3 * w(7,1) + Eignung4 * w(8,1) + Eignung5 * w(9,1);
    
    
    % Belohnung
    Belohnung_Pausen = Pausen * Population{z1,2}(1,1) * w(3,1);
    Belohnung_Werkzeuge = Werkzeuge * Population{z1,2}(1,1) * w(4,1);
    Population{z1,2}(1,1) = Population{z1,2}(1,1) + Belohnung_Pausen + Belohnung_Werkzeuge; % Belohnung der Fitness
end


% Abfrage, ob Taktzeit überschritten wird
for z1 = 1:p
    if max(Population{z1,1}(6,:)) > t_takt % wird Taktzeit überschritten?
        Population{z1,2}(1,1) = Population{z1,2}(1,1) / 2;  % wenn ja, Fitness halbieren
    end
end

