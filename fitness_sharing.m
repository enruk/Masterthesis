function [Population] = fitness_sharing (Population,p,hamming_abstand,delta_startzeit)

for z1 = 1:p                    % Nachbarn zur�cksetzen / nullen
    Population{z1,8}(1,1) = 0;
end


% Z�hlen der Nachbarn

    % �u�ere Schleife: Individuen hochz�hlen
for z1 = 1:p       % Individuen 1 hochz�hlen
    anz_Nachbar = 0;
    rang = Population{z1,5}(1,1);    % Speichern des Rangs des Individuums
    for z2 = 1:p            % Individuen 2 hochz�hlen
        if Population{z2,5}(1,1) == rang     % Individuum 2 mit gleichem Rang suchen wie Individuum 1
            if Population{z1,7}(z2,1) < hamming_abstand && Population{z1,7}(z2,2) < delta_startzeit && z1 ~= z2       % wenn CD klein genug
                anz_Nachbar = anz_Nachbar + 1;
                Population{z1,8}(1,1) = anz_Nachbar;                        % dann Nachbar hinzuf�gen zu Individuum 1
            elseif Population{z1,7}(z2,1) < 1 && z1 ~= z2       % wenn HA klein genug
                anz_Nachbar = anz_Nachbar + 1;
                Population{z1,8}(1,1) = anz_Nachbar;                        % dann Nachbar hinzuf�gen zu Individuum 1
            elseif Population{z1,7}(z2,2) < 2 && z1 ~= z2       % wenn CD klein genug
                anz_Nachbar = anz_Nachbar + 1;
                Population{z1,8}(1,1) = anz_Nachbar;                        % dann Nachbar hinzuf�gen zu Individuum 1
            end
        end
    end
    % n�chstes Individuum
end


    % Angepasste Fitness in Spalte 6 eintragen
for z1 = 1:p
    Population{z1,6}(1,1) = Population{z1,2}(1,1) / (Population{z1,8}(1,1) + 1);
end  