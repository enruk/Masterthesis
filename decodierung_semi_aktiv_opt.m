function [Population] = decodierung_semi_aktiv_opt (Population,n,T,A_ges,Prae)

run('config')
run('trans')

%% Initialisierung der Decodierung

D = cell (p,1);             % Duration / (Prozess-)Dauer
S_p = cell(p,1);            % Enth�lt f�r jede Operation k die Prozesszeiten der Vorg�nger (nur 1 Schritt zur�ck)
S_min = cell(p,1);          % Aufsummieren der Prozesszeit aller Vorg�nger, um t_min zu erhalten
for z1 = 1:p                 % Auff�llen
    S_p{z1,1} = zeros(n,n);
    S_min{z1,1} = zeros(n,n);
end


% Matrizen / Vektoren f�r GT-Algorithmus
A = zeros(n,1);             % Einheitsvektor A der Anfangsoperationen (Prozessabh�ngig)
A_x = cell(p,1);            % Vektor der Anfangsoperationen 


F = cell(p,1);              % Matrix der Fertigungszeiten jeder Operation
for z1 = 1:p
    F{z1,1} = zeros(n,1);
end
   
FA_x = cell(p,1);           % Matrix FA_x der Fertigungszeiten der Operationen in A

for z1 = 1:p
    FA_x{z1,1} = zeros(n,1);
end   


B_x = cell(p,1);            % Matrix B_x der Operationen, die Element von A sind und auf der
                            % gleichen Ressource wie O' sind
for z1 = 1:p
    B_x{z1,1} = zeros(n,1);
end


V_x = cell(p,1);            % Matrix V_x der Vorg�nger, zur �berpr�fung ob ein Prozess ausgef�hrt
                            % werden kann (0), schon ausgef�hrt wurde (NaN) oder noch warten muss
                            % (>0)      
for z1 = 1:p
    V_x{z1,1} = A_ges;
end

for z1 = 1:p
    V_x{z1}(V_x{z1}==0) = zeros;
end

operation = cell(p,1);      % Vektor der Operationen f�r jedes Individuum
                            % 1,1 ist Operation'
                            % 2,1 ist Operation''
                            % 3,1 ist Operation*
for z1 = 1:p
    operation{z1,1} = zeros(3,1);
end



% Startzeiten Berechnung
AO = cell(p,1);             % Vektor f�r Startoperationen zur Berechnung der Startzeiten
VO = cell(p,1);             % Vektor der Vorg�nger zur Ermittlung von AO



%% Giffler-Thompson Algorithmus
% Decodierung: Gew�hlter Ansatz: Hybrides scheduling


%% Schritt 0: Vorbereitung und Initialisierung
    

% Erstellung der Prozesszeitenmatrix D f�r jedes Individuum
    % T sind die Prozessdaten, ausgelesen in data.m aus der Excel
    % Matrizen D mit T bef�llen
    for z1 = 1:p
        D{z1}(:,1) = T(1:n,1);
        D{z1}(:,2) = T(1:n,2);
    end

    % f�r den Menschen (Prozesse auf dem Roboter null setzen)
    for z1 = 1:p
         dm = D{z1}(:,1);
         m = Population{z1}(2,:)';
         dm(m==0) = 0;
         D{z1}(:,1) = dm;
    end

    % f�r den Roboter (Prozesse auf dem Menschen null setzen)
    for z1 = 1:p
         dr = D{z1}(:,2);
         r = Population{z1}(3,:)';
         dr(r==0) = 0;
         D{z1}(:,2) = dr;
    end
    

%% Erstellung fr�hstm�gliche Startzeiten (S_min) und Fertigungszeiten F 

% Erstellung der fr�hstm�glichen Startzeiten der einzelnen Prozesse

    % f�r Prozess k werden die Prozesszeiten der Vorg�nger Prozesse k-1
    % (nicht k-2) eingetragen in S_p

    for k = 1:p
        for z = 1:n
            for s = 1:n
                if Prae(z,s) == 1 
                    if Population{k}(2,s) == 1 
                        S_p{k}(z,s) = D{k}(s,1);
                    end
                    if Population{k}(3,s) == 1
                        S_p{k}(z,s) = D{k}(s,2);
                    end
                end
            end
        end
    end


    % fr�hstm�gliche Startzeit ist das Maximum von 
    % t = (Startzeit Vorprozess + Dauer Vorprozess, 
    % Startzeit + Prozesszeit des letzten Prozesses auf der Ressource)
    % 
    % Startzeit + Prozesszeit des letzten Prozesses auf der Ressource ist
    % nicht bekannt, da noch keine Reihenfolge fest zugeordnet wurde
    %
    % Startzeit Vorprozess + Dauer Vorprozess kann vorl�ufig bestimmt
    % werden
    
    % WICHTIG: Sobald erste Prozese im Giffler-Thompson Algorithmus
    % festgelegt werden, muss S_min aktualisiert werden
    % Genauere Erl�uterung in Masterarbeit
    
    % Prozessdauer des Vorprozesses + sp�teste Startzeit des Vorprozesses
    % ist gleich der fr�hstm�glichen Startzeit Prozessabh�ngig

%     for z1 = 1:p
%         for z2 = 1:n
%             for s = 1:n
%                 if S_p{z1,1}(z2,s) ~= 0
%                     S_min{z1,1}(z2,s) = S_p{z1,1}(z2,s) + max(S_min{z1,1}(s,:)); 
%                 end
%             end
%         end
%     end
    
    
    %% 
    for z1 = 1:p
        
        % AO und VO nullen
        AO{z1,1} = zeros(n,1);
        VO{z1,1} = A_ges;
        VO{z1}(VO{z1}==0) = zeros;


        % Allgmeines Vorgehen: Der Vorranggraph wird von vorne nach hinten abgearbeitet
        % Anfangsoperationen ermitteln
        for z2 = 1:n
            if max(S_p{z1,1}(z2,1:n)) < 0.5   % Bedingung f�r Anfangsoperation 
               AO{z1,1}(z2,1) = 1;
            end
        end


        anz_op = 0;   
        while anz_op < n

            % Startzeit der Prozesse in AO: t = t_v + p_v
            for z2 = 1:n
                if AO{z1,1}(z2,1) == 1  % Anfangsoperation gefunden
                    anf_op = z2;        % Anfangsoperation neu abspeichern
                    for z3 = 1:n        % Nachfolger suchen
                        if S_p{z1,1}(z3,anf_op) ~= 0        % Bedingung f�r Nachfolger gefunden
                            S_min{z1,1}(z3,anf_op) = S_p{z1,1}(z3,anf_op) + max(S_min{z1,1}(anf_op,:)); % Maximum des Vorg�ngers draufrechnen
                        end
                    end
                VO{z1,1}(anf_op,:) = NaN;
                VO{z1,1}(VO{z1,1}==anf_op) = 0;
                AO{z1,1}(anf_op,1) = NaN;
                end
            end

            % AO neu bestimmen
            for z2 = 1:n
                if max(VO{z1}(z2,:)) == 0      % Fall 1: alle Vorg�nger wurden ausgef�hrt
                    AO{z1}(z2,1) = 1;          % Dann: Operation zu AO hinzuf�gen
                end
                if isnan(max(VO{z1}(z2,:)))      % Fall 2: Operation wurde schon ausgef�hrt
                    AO{z1}(z2,1) = NaN;          % Dann: Als NaN kennzeichnen
                end
                if max(VO{z1}(z2,:)) > 0      % Fall 3: alle Vorg�nger wurden noch nicht ausgef�hrt
                    AO{z1}(z2,1) = 0;          % Dann: Operation nicht zu AO hinzuf�gen
                end
            end

            % Wenn AO leer, dann fertig
            anzahl_op = find(isnan(AO{z1}(1:n,1)));
            [anz_op,~] = size(anzahl_op);
        end
    end
    
    
    %%
    % Fertigungszeiten: Maximum der vorl�ufigen Startzeiten + Prozessdauer
    for z1 = 1:p
        for z2 = 1:n
            F{z1,1}(z2,1) = max(S_min{z1,1}(z2,:)) + D{z1,1}(z2,1) + D{z1,1}(z2,2);
        end
    end

%% GT - Schritt 1: Vektor A bestimmen mit Anfangsoperationen
    % Erstellung der Matrix der Anfangsoperationen A (f�r alle Individuen
    % hier noch gleich gleich, da nur vom Prozess abh�ngig, nicht vom Individuum)

    for z2 = 1:n
        if max(A_ges(z2,:)) < 0.1
            A(z2,1) = 1;
        end
    end

    % Matrix A_x der Anfangsoperationen f�r jedes Individuum
    % theoretisch ist keine Cell notwendig, eine Matrix reicht
    % hier gew�hlt: jedes Indivdiuum bekommt eine eigene Matrix f�r A_x,
    % B_x, FA_x, V_x
    for z1 = 1:p
        A_x{z1,1} = A;
    end


%% GT - Schritt 2: Beginn der Schleife

 
% Erst gesamter Ablauflauf Individuum 1, dann Individuum 2 usw.)

for z1=1:p
    
    % Vorbereitung f�r neues Individuum
    t_beleg_m = 0;     % Aktuelle Belegungszeit Mensch zur�cksetzen
    t_beleg_r = 0;     % Aktuelle Belegungszeit Roboter zur�cksetzen
    anz_nan = 0;

    %%
    
    while anz_nan < n
        
               
        operation{z1} = zeros(3,1);
        % Beim semi-aktiven werden die ersten beiden Zeilen nicht
        % benutzt
         
        
        %% GT - Schritt 2: Operation O* ermitteln

        permu = 1000;
        for z2 = 1:n
            if A_x{z1}(z2,1) == 1
                
                permu_neu = Population{z1}(1,z2);
                
                if permu_neu < permu
                    permu = permu_neu;
                    operation{z1}(3,1) = z2;    % Position von O* speichern
                end
            end
        end
        
        %% Ressource und Belegungszeit bestimmen
        
        % Ermitteln der Ressource von O*
        if Population{z1}(2,operation{z1}(3,1)) == 1
            ressource = 2;      % Ressource = Mensch
        else 
            ressource = 3;      % Ressource = Roboter
        end
        
        
        % Belegungszeit ermitteln
        % Belegungszeit wird nur zur Sicherheit bestimmt
        % Die Startzeiten aller Operationen auf einer Ressource werden
        % aktualisiert und k�nnen so im Normalfall nicht vor der
        % Beleugngszeit beginnnen
        
        if ressource == 2
            t_beleg_res = t_beleg_m;
        else
            t_beleg_res = t_beleg_r;
        end
        
            % Operation O* aus A l�schen
        A_x{z1}(operation{z1}(3,1),1) = NaN;      % Aufgaben die durchgef�hrt sind mit NaN kennzeichnen
        V_x{z1}(V_x{z1}==operation{z1}(3,1)) = 0;
        V_x{z1}(operation{z1}(3,1),:) = NaN;
        
        %% GT - Schritt 3: O* dem Ablaufplan hinzuf�gen


              % Startzeit in P eintragen
        if max(S_min{z1}(operation{z1}(3,1),:)) < t_beleg_res
            Population{z1}(4,operation{z1}(3,1)) = t_beleg_res;
        end
            
        if max(S_min{z1}(operation{z1}(3,1),:)) >= t_beleg_res
            Population{z1}(4,operation{z1}(3,1)) = max(S_min{z1}(operation{z1}(3,1),:)); 
        end
        
            % Prozesszeit in P eintragen
        Population{z1}(5,operation{z1}(3,1)) = D{z1,1}(operation{z1}(3,1),1) + D{z1,1}(operation{z1}(3,1),2);
        
            % Fertigungszeit in P eintragen
        Population{z1}(6,operation{z1}(3,1)) = Population{z1}(4,operation{z1}(3,1)) + Population{z1}(5,operation{z1}(3,1));
            
            % Belegungsplan aktualisieren
        if ressource == 2
            t_beleg_m = Population{z1}(6,operation{z1}(3,1)); % t_beleg_m + D{z1}(operation{z1}(3,1),1) + D{z1}(operation{z1}(3,1),2);
        else
            t_beleg_r = Population{z1}(6,operation{z1}(3,1)); %t_beleg_r + D{z1}(operation{z1}(3,1),1) + D{z1}(operation{z1}(3,1),2);
        end
        
        
        
        %% GT - Schritt 4: Nachfolger von O* zu A hinzuf�gen
        
            % Nachfolger von O* zu A hinzuf�gen
            
            
        for z2 = 1:n
            if max(V_x{z1}(z2,:)) == 0      % Fall: alle Vorg�nger wurden ausgef�hrt
                A_x{z1}(z2,1) = 1;          % Dann: Operation zu A hinzuf�gen
            end
            if isnan(max(V_x{z1}(z2,:)))      % Fall: alle Vorg�nger wurden ausgef�hrt
                A_x{z1}(z2,1) = NaN;          % Dann: Operation zu A hinzuf�gen
            end
            if max(V_x{z1}(z2,:)) > 0      % Fall: alle Vorg�nger wurden noch nicht ausgef�hrt
                A_x{z1}(z2,1) = 0;          % Dann: Operation nicht zu A hinzuf�gen
            end
        end
               
        %% GT - Schritt 5: Belegungszeiten aktualisieren
        
        % Anpassen der Matrix S_min
        for z2 = 1:n
            if ~isnan(A_x{z1}(z2,1)) 
                if Population{z1}(ressource,z2) == 1
                    S_min{z1}(z2,n+1) = Population{z1}(6,operation{z1}(3,1));    % Neue Startzeit f�r Operationen derselben Ressource
                end
            end
        end
        
%         for z2 = 1:n
%             if S_min{z1}(z2,operation{z1}(3,1)) > 0
%                 S_min{z1}(z2,n+1) = Population{z1}(6,operation{z1}(3,1));        % Neue Startzeit f�r Operationen die Nachfolger von O* sind
%             end
%         end
        
        
        %% GT - Schritt 6: Neuberechnung der Matrix S_min
            % 
%         for z2 = 1:n
%             for s = 1:n
%                 if S_p{z1}(z2,s) ~= 0
%                     S_min{z1}(z2,s) = S_p{z1}(z2,s) + max(S_min{z1}(s,:));
%                 end
%             end
%         end
        
        % VO und AO nullen
        AO{z1,1} = zeros(n,1);
        VO{z1,1} = A_ges;
        VO{z1}(VO{z1}==0) = zeros;
        

        % Allgmeines Vorgehen: Der Vorranggraph wird von vorne nach hinten abgearbeitet
        % Anfangsoperationen ermitteln
        for z2 = 1:n
            if max(S_p{z1,1}(z2,1:n)) < 0.5   % Bedingung f�r Anfangsoperation 
               AO{z1,1}(z2,1) = 1;
            end
        end


        anz_op = 0;   
        while anz_op < n

            
            % Startzeit der Prozesse in AO: t = t_v + p_v
            for z2 = 1:n
                if AO{z1,1}(z2,1) == 1  % Anfangsoperation gefunden
                    anf_op = z2;        % Anfangsoperation neu abspeichern
                    for z3 = 1:n        % Nachfolger suchen
                        if S_p{z1,1}(z3,anf_op) ~= 0        % Bedingung f�r Nachfolger gefunden
                            S_min{z1,1}(z3,anf_op) = S_p{z1,1}(z3,anf_op) + max(S_min{z1,1}(anf_op,:)); % Maximum des Vorg�ngers draufrechnen
                        end
                    end
                VO{z1,1}(anf_op,:) = NaN;
                VO{z1,1}(VO{z1,1}==anf_op) = 0;
                AO{z1,1}(anf_op,1) = NaN;
                end
            end


            % AO neu bestimmen
            for z2 = 1:n
                if max(VO{z1}(z2,:)) == 0      % Fall 1: alle Vorg�nger wurden ausgef�hrt
                    AO{z1}(z2,1) = 1;          % Dann: Operation zu AO hinzuf�gen
                end
                if isnan(max(VO{z1}(z2,:)))      % Fall 2: Operation wurde schon ausgef�hrt
                    AO{z1}(z2,1) = NaN;          % Dann: Als NaN kennzeichnen
                end
                if max(VO{z1}(z2,:)) > 0      % Fall 3: alle Vorg�nger wurden noch nicht ausgef�hrt
                    AO{z1}(z2,1) = 0;          % Dann: Operation nicht zu AO hinzuf�gen
                end
            end

            
            % Wenn AO leer, dann fertig
            anzahl_op = find(isnan(AO{z1}(1:n,1)));
            [anz_op,~] = size(anzahl_op);
        end
        
        %% GT - Schritt 7: Neuberechnung der Matrix F
            % Neuberechnung der Matrix F
        
        for z2 = 1:n
            F{z1}(z2,1) = max(S_min{z1}(z2,:)) + D{z1}(z2,1) + D{z1}(z2,2);
        end
        

        %% GT - Ende: Abbruchkriterium neu berechnen
            % Abbruchkriterium anpassen
        anzahl_nan = find(isnan(A_x{z1}(1:n,1)));
        [anz_nan,~] = size(anzahl_nan);
        
        % n�chster Durchlauf des GT Algorithmus
    end
    
    
    % hier m�sste dann das Zur�cksetzen der Matrizen A_x, B_x, FA_x und V_x
    % stattfinden, wenn diese nicht als Zelle definiert werden
    
    % n�chstes individuum    
end

