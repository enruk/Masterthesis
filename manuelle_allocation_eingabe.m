function [Population,geneM,geneR,permutation] = manuelle_allocation_eingabe(Population,n)



%% Zuorndung auf Ressourcen

% Mensch
str = inputdlg('Sollen Operationen dem Mensch fest zuordnet werden? (j/n)', 'Zuordnungsgen');
if str{1,1}(1,1) == 'j'
    prompt = {'Feste Operationen (Leerzeichen als Trennung)'};
    title = 'Feste Zuordnung für den Menschen';
    dims = [1 35];
    answer = inputdlg(prompt,title,dims);
    
    geneM(1,:) = str2num(regexprep(answer{1,1}, '\D+', ' '));
else
    geneM = zeros(2,1);     % Wenn nicht, dann Nullvektor erstellen
end

% Roboter
str = inputdlg('Sollen Operationen dem Roboter fest zuordnet werden? (j/n)', 'Zuordnungsgen');
if str{1,1}(1,1) == 'j'
    prompt = {'Feste Operationen (Leerzeichen als Trennung)'};
    title = 'Feste Zuordnung für den Roboter';
    dims = [1 35];
    answer = inputdlg(prompt,title,dims);
    
    geneR(1,:) = str2num(regexprep(answer{1,1}, '\D+', ' '));
else
    geneR = zeros(1,1);     % Wenn nicht, dann Nullvektor erstellen
end


%% Zuordnung auf Sequenz
str = inputdlg('Sollen einer Operation eine feste Permutation festzugeordnet werden? (j/n)', 'Sequenzgenom');
index  = 0;
permutation = zeros(2,n);

while str{1,1}(1,1) ~= 'n'
    prompt = {'Nummer der Operation','Feste Permutation'};
    title = 'Feste Zuordnung von Operationen';
    dims = [1 35;1 35];
    answer = inputdlg(prompt,title,dims);
    index = index + 1;
    permutation(1,index) = str2double(answer{1,1}(1,1)); 
    permutation(2,index) = str2double(answer{2,1}(1,1));
    str = inputdlg('Sollen einer weiteren Operation eine feste Permutation festzugeordnet werden? (j/n)', 'Sequenzgenom');
end

[~,s] = find(permutation(1,:) == 0);
if s(1,1) ~= 1
    permutation(:,s(1,1):end) = [];
elseif s(1,1) == 1
    permutation(:,2:end) = [];
end
