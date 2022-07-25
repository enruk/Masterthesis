function [Population] = decodierung_opt (Population,n,T,A_ges,Prae,op_kolla)

run('config')
run('trans')

%% Initialisierung der Decodierung

D = cell (p,1);             % Duration / (Prozess-)Dauer
S_p = cell(p,1);            % Enth�lt f�r jede Operation k die Prozesszeiten der Vorg�nger (nur 1 Schritt zur�ck)
S_min = cell(p,1);          % Aufsummieren der Prozesszeit aller Vorg�nger, um t_min zu erhalten
for z1 = 1:p                % Auff�llen
    S_p{z1,1} = zeros(n,n);
    S_min{z1,1} = zeros(n,n+1);
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
% Decodierung: W�hlbar zwischen semi-aktiv und hybrid in config


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
         dm = D{z1}(:,1);           % Dauer Mensch
         m = Population{z1}(2,:)';  % Verteilung Mensch 
         dm(m==0) = 0;              % Wenn Verteilungschromosom =0, dann Dauer in D auf 0 setzen
         D{z1}(:,1) = dm;           % neu abspeichern
    end

    % f�r den Roboter (Prozesse auf dem Menschen null setzen)
    for z1 = 1:p
         dr = D{z1}(:,2);           % Dauer Roboter
         r = Population{z1}(3,:)';  % Verteilung Roboter
         dr(r==0) = 0;              % Wenn Verteilungschromosom =0, dann Dauer in D auf 0 setzen
         D{z1}(:,2) = dr;           % neu abspeichern
    end
    
    
    % kollaborative Operationen besitzen nun doppelte Zeit
    for z1 = 1:p
        for s = 1:length(op_kolla)
            D{z1}(op_kolla(1,s),1) = D{z1}(op_kolla(1,s),1) / 2;    % daher jeweils teilen
            D{z1}(op_kolla(1,s),2) = D{z1}(op_kolla(1,s),2) / 2;
        end
    end
    
    
%% Erstellung der fr�hstm�glichen Startzeiten (S_min) und Fertigungszeiten F 


% Erstellung der fr�hstm�glichen Startzeiten der einzelnen Prozesse

    % f�r Prozess k werden die Prozesszeiten der Vorg�nger Prozesse k-1
    % (nicht k-2) eingetragen in S_p

    for z1 = 1:p
        for z = 1:n         % Zeilen hochz�hlen
            for s = 1:n         % Spalten hochz�hlen
                if Prae(z,s) == 1   % Wenn Pr�zedenzbeziehung
                    if Population{z1}(2,s) == 1     % Dauer vom Mensch 
                        S_p{z1}(z,s) = D{z1}(s,1);
                    end
                    if Population{z1}(3,s) == 1     % oder vom Roboter   
                        S_p{z1}(z,s) = D{z1}(s,2);
                    end
                end
            end
        end
    end

    %%

    % fr�hstm�gliche Startzeit ist das Maximum von 
    % t = (Startzeit Vorprozess + Dauer Vorprozess, 
    % Startzeit + Prozesszeit des letzten Prozesses auf der Ressource)
    % 
    % Startzeit + Prozesszeit des letzten Prozesses auf der Ressource ist
    % nicht bekannt, da noch keine Reihenfolge fest zugeordnet wurde
    %
    % Startzeit Vorprozess + Dauer Vorprozess kann vorl�ufig bestimmt
    % werden
    %
    % WICHTIG: Sobald erste Prozese im Giffler-Thompson Algorithmus
    % festgelegt werden, muss S_min aktualisiert werden
    % Genauere Erl�uterung in Masterarbeit
    %
    % Prozessdauer des Vorprozesses + sp�teste Startzeit des Vorprozesses
    % ist gleich der fr�hstm�glichen Startzeit Prozessabh�ngig    
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
                if max(VO{z1}(z2,:)) == 0       % Fall 1: alle Vorg�nger wurden ausgef�hrt
                    AO{z1}(z2,1) = 1;           % Dann: Operation zu AO hinzuf�gen
                end
                if isnan(max(VO{z1}(z2,:)))     % Fall 2: Operation wurde schon ausgef�hrt
                    AO{z1}(z2,1) = NaN;         % Dann: Als NaN kennzeichnen
                end
                if max(VO{z1}(z2,:)) > 0        % Fall 3: alle Vorg�nger wurden noch nicht ausgef�hrt
                    AO{z1}(z2,1) = 0;           % Dann: Operation nicht zu AO hinzuf�gen
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
    t_beleg_m = 0;      % Aktuelle Belegungszeit Mensch zur�cksetzen
    t_beleg_r = 0;      % Aktuelle Belegungszeit Roboter zur�cksetzen
    anz_nan = 0;        % Abbruckbedingungen zur�cksetzen
    
    %%
    
    while anz_nan < n
        
        % Operation O',O'' und O* nullen 
        operation{z1} = zeros(3,1);
        
        
        %% GT - Schritt 2.H1: Ermittlung von O'
        
            % Ermittlung der Operation mit der k�rzesten Fertigungszeit 
        for z2 = 1:n
            if A_x{z1}(z2,1) == 1   % Wenn Anfangsoperation
                FA_x{z1}(z2,1) = F{z1}(z2,1);   % Dann in FA_x eintragen
            else
                FA_x{z1}(z2,1) = NaN;           % Andernfalls nicht
            end
        end

        FA_x{z1}(FA_x{z1}==0)=NaN;      % Nullen durch NaN ersetzen in FA_x
        [~,y] = min(FA_x{z1}(:,1));     % Minimale Fertigungszeit suchen
        operation{z1}(1,1) = y;         % Nummer von O' speichern
        
        %% GT - Schritt 2.H2: Vektor B bestimmen
        
            % Abfrage, ob O' eine kollaborative Operation ist
        answer = any(op_kolla == operation{z1}(1,1));
        if answer == 1          % Wenn ja, 
            if t_beleg_m > t_beleg_r 
                ressource = 3;  % k�rzer belegte Ressource wird bevorzugt
            else 
                ressource = 2;
            end
        end
        
            % Ermitteln der Ressource von Operation_Strich
        if answer == 0              % Wenn nein,
            if Population{z1}(2,operation{z1}(1,1)) == 1    % Ressource bestimmen
                ressource = 2;      % Ressource = Mensch
            else 
                ressource = 3;      % Ressource = Roboter
            end
        end
        
            % Erstellung des Vektors B
            % Bedingungen: Operations ist Anfangsoperation und auf der gleichen
            % Ressource wie Operation_Strich
        B_x{z1} = zeros(n,1);
        
        for z2 = 1:n
            if A_x{z1}(z2,1) == 1 && Population{z1}(ressource,z2) == 1
                B_x{z1}(z2,1) = 1;
            end
        end 
        
        %% GT - Schritt 2.H3: Ermittlung von O''
        
        
            % Ermittlung der Operation von B mit der fr�hstm�glichen
            % Startzeit
        t_min = 1000;       % Vergleichsparameter hochsetzen
        for z2 = 1:n
            if B_x{z1}(z2,1) == 1
                t_neu = max(S_min{z1}(z2,:));
                if t_neu < t_min
                    t_min = t_neu;
                    operation{z1}(2,1) = z2;    % Nummer von O'' speichern
                end
            end
        end

        %% GT - Schritt 2.H4: Operationen au�erhalb des Zeitfensters entfernen
        
             % Ermittlung der Operationen, die die Bedingung 
             % t < t''+sigma((t'+p')-t'') nicht erf�llen und l�schen aus B
             
        % Belegungszeit der Ressource separat speichern, damit nicht immer
        % die Ressource abgefragt werden muss
        if ressource == 2       % Die Ressource wurde vorher bestimmt
            t_beleg_res = t_beleg_m;    % Mensch
        else
            t_beleg_res = t_beleg_r;    % Roboter
        end

        
        
        % M�glich, dass dieser Schritt unn�tigt ist
        % Die Belegungszeit sind in den Startzeiten ja ber�cksichtigt
        if  max(S_min{z1}(operation{z1}(2,1),:)) < t_beleg_res  
            t_2str = t_beleg_res;
        else
            t_2str =  max(S_min{z1}(operation{z1}(2,1),:));
        end
        
        
        % Entfernen der Operationen au�erhalb des Zeitfensters
        for z2 = 1:n
            if B_x{z1}(z2,1) == 1
                if max(S_min{z1}(z2,:)) > t_2str + dc_sigma * ( max(F{z1}(operation{z1}(1,1),1)) -  t_2str )
                    B_x{z1}(z2,1) = 0;
                end               
            end
        end
        
        
        %% GT - Schritt 2.H5: Operation O* ermitteln
        
        % Giffler-Thompson: Schritt 2.H5
        permu = 1000;
        for z2 = 1:n
            if B_x{z1}(z2,1) == 1
                
                permu_neu = Population{z1}(1,z2);
                
                if permu_neu < permu
                    permu = permu_neu;
                    operation{z1}(3,1) = z2;    % Position von O* speichern
                end
            end
        end
            
            % Operation O* aus A l�schen
        A_x{z1}(operation{z1}(3,1),1) = NaN;        % Operation O* mit NaN kennzeichnen (wurde ausgef�hrt)
        V_x{z1}(V_x{z1}==operation{z1}(3,1)) = 0;   % Dort wo O* Vorg�nger war, als ausgef�hrt kennzeichnen
        V_x{z1}(operation{z1}(3,1),:) = NaN;        % Operation O* besitzt keine Vorg�nger mehr und wurde ausgef�hrt
        
        
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
            % Wenn kollaborativ, beide Ressourcen aktualisieren
        answer = any(op_kolla == operation{z1}(3,1));
        if answer == 1
            t_beleg_m = Population{z1}(6,operation{z1}(3,1));
            t_beleg_r = t_beleg_m;
        end
            
            % Wenn nicht kollaborativ, nur die betreffende
        if answer  == 0
            if ressource == 2
                t_beleg_m = Population{z1}(6,operation{z1}(3,1)); % t_beleg_m + D{z1}(operation{z1}(3,1),1) + D{z1}(operation{z1}(3,1),2);
            else
                t_beleg_r = Population{z1}(6,operation{z1}(3,1)); %t_beleg_r + D{z1}(operation{z1}(3,1),1) + D{z1}(operation{z1}(3,1),2);
            end
        end
        
        
        %% GT - Schritt 4: Nachfolger von O* zu A hinzuf�gen
        
            % Nachfolger von O* zu A hinzuf�gen
            
            
        for z2 = 1:n
            if max(V_x{z1}(z2,:)) == 0          % Fall: alle Vorg�nger wurden ausgef�hrt
                A_x{z1}(z2,1) = 1;              % Dann: Operation zu A hinzuf�gen
            end
            if isnan(max(V_x{z1}(z2,:)))        % Fall: Operation wurde schon ausgef�hrt
                A_x{z1}(z2,1) = NaN;            % Dann: Operation nicht zu A hinzuf�gen
            end
            if max(V_x{z1}(z2,:)) > 0           % Fall: alle Vorg�nger wurden noch nicht ausgef�hrt
                A_x{z1}(z2,1) = 0;              % Dann: Operation nicht zu A hinzuf�gen
            end
        end
               
        %% GT - Schritt 5: Belegungszeiten aktualisieren
        
            % Anpassen der Matrix S_min
            
            % Wenn kollaborativ, alle verbleibenden Operationen erhalten
            % neue Belegungszeit
        answer = any(op_kolla == operation{z1}(3,1));
        if answer == 1
            for z2 = 1:n
                if ~isnan(A_x{z1}(z2,1)) 
                    S_min{z1}(z2,n+1) = Population{z1}(6,operation{z1}(3,1));    % Neue Startzeit f�r Operationen derselben Ressource
                end
            end
        end
        
            % Wenn nicht kollaborativ, nur die Operationen auf derselben
            % Ressource wie O* neue Belegungszeit
        if answer == 0
            for z2 = 1:n
                if ~isnan(A_x{z1}(z2,1)) 
                    if Population{z1}(ressource,z2) == 1
                        if S_min{z1}(z2,n+1) < Population{z1}(6,operation{z1}(3,1))     % Belegungszeit nur aktualisieren, wenn neue gr��er als alte ist (f�r kollaborative Operationen wichtig, da auf zwei Ressourcen)
                            S_min{z1}(z2,n+1) = Population{z1}(6,operation{z1}(3,1));    % Neue Startzeit f�r Operationen derselben Ressource
                        end
                    end
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
%         for z2 = 1:n          % altes Vorgehen
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

