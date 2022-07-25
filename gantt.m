function [Duration,Positions,mensch_sort,roboter_sort,anz_mensch,anz_roboter] = gantt(Population,individuum)


%Gantt-Diagramm 

prozesse_m = find(Population{individuum,1}(2,:)==1);
mensch = [prozesse_m ; Population{individuum,1}(1:6,prozesse_m)];
[~, index] = sort(mensch',1);
mensch_sort = mensch(:,index(:,5));

prozesse_r = find(Population{individuum,1}(3,:)==1);
roboter = [prozesse_r ; Population{individuum,1}(1:6,prozesse_r)];
[~, index] = sort(roboter',1);
roboter_sort = roboter(:,index(:,5));

[~,anz_mensch] = size(mensch_sort);
[~,anz_roboter] = size(roboter_sort);



%%

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

% kürzere Ressource mit nullen auffüllen
if length(duration_roboter)<length(duration_mensch)
    duration_roboter(end+1:length(duration_mensch)) = 0;
elseif length(duration_roboter)>length(duration_mensch)
    duration_mensch(end+1:length(duration_roboter)) = 0;
end

Duration = [duration_mensch;duration_roboter];

Positions = [1,2];

