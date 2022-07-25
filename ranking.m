function [Population] = ranking(Population,p,anz_raenge,R_max)

for z1 = 1:p
    Population{z1,4}(1,1) = 0;
    Population{z1,5}(1,1) = 0;
end

    % Bestes Individuum bestimmen
best = 0.01;
for z1 = 1:p
    fit = Population{z1,2}(1,1);
    if fit > best
        best = fit;
    end
end

    % Schlechtestes Individuum bestimmen
worst = 1;
for z1 = 1:p
    fit = Population{z1,2}(1,1);
    if fit < worst
        worst = fit;
    end
end



range = best - worst;
rang_range = range / anz_raenge;    % Berechnung der Breite eines Rangs

        % Äußere Schleife: Ränge hochzählen (1 bis anz_raenge)
        % Innere Schleife: Individuen hochzählen
for rang = 1:anz_raenge
    for z1 = 1:p
        if Population{z1,2}(1,1) >= worst + (rang-1)*rang_range 
            Population{z1,5}(1,1) = rang;
        end
    end
end



% angepasste Fitness wird in Spalte 4 geschrieben

for z1=1:p     % Nach Baker (Buttelmann, Lohmann 2004, S.155)
    Population{z1,4}(Population{z1,5}==anz_raenge) = R_max; % höchster Rang
    if anz_raenge ~= 1
        Population{z1,4}(Population{z1,5}==1) = 2 - R_max;      % niedrigster Rang
        Population{z1,4}(1,1) = (2-R_max) + (R_max - (2-R_max))*(Population{z1,5}(1,1)-1)/(anz_raenge-1);   % Zwischenränge
    end
end