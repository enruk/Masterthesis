function Kinder = invert_mutation(Kinder,p,n,genom,mut_whk)

for z1 = 1:p
    random = rand(1);                       % Zufallszahl bestimmen
    
    if random < mut_whk                     % Wenn random kleiner als Mutationswahrscheinlichkeit ist, dann wird mutiert
        
        schnitt1 = round(n*rand(1));        % Schnitt1 zufällig wählen
        while schnitt1 == 0
            schnitt1 = round(n*rand(1));
        end
        
        schnitt2 = round(n*rand(1));        % Schnitt2 zufällig wählen
        while schnitt1 == schnitt2 || schnitt2 == 0
            schnitt2 = round(n*rand(1));    % Schnitt2 wenn er gleich schnitt1 ist
        end
        
        if schnitt2 < schnitt1              % Schnitte tauschen wenn Schnitt2 vor Schnitt1
            [schnitt1,schnitt2] = deal(schnitt2,schnitt1);
        end
        
        abschnitt = Kinder{z1,1}(genom,schnitt1:schnitt2);  % Abschnitt aus Kind auslesen
        
        abschnitt = fliplr(abschnitt);                      % Abschnitt invertieren
        
        Kinder{z1,1}(genom,schnitt1:schnitt2) = abschnitt;  % Abschnitt Kind wieder Kind zurückgeben
    end
    
end