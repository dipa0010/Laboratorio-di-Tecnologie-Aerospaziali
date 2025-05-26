%% Analisi spessore layer da CSV — processing file singoli
% Script completo per misurazione automatica degli spessori dei layer
% da CSV esportati da Fiji (o programmi simili che permettono 
% il Plot Profile delle intensità di colore di un profilo).
% Richiede MATLAB R2019b o successivo.
% Non si è avuta la possibilità di utilizzare "Signal Processing Toolbox",
% per questo motivo è stata implementata una funzione findpeaks ad hoc per
% il caso in questione. 

close all; clear; clc;

%%
files = { ...
  'layer1_1.csv', 'layer1_2.csv', 'layer1_3.csv',...
  'layer2_1.csv','layer2_2.csv','layer2_3.csv', ...
  'layer3_1.csv','layer3_2.csv','layer3_3.csv','layer3_4.csv','layer3_5.csv' ...
};


%%
distanze = calcola_distanze_layer(files, 'distanze_layers.csv');

 

%%
function distanze = calcola_distanze_layer(files, output_csv)
% Trova tutti i massimi locali e filtra quelli principali per distanza e soglia y>100
% Stampa a video il numero di massimi locali, massimi principali e distanze misurate per ogni file

distanza_nominale = 0.05;
tolleranza = 0.025;
min_dist = distanza_nominale - tolleranza;
max_dist = distanza_nominale + tolleranza;
soglia_y = 100; % filtro su y

if nargin < 2
    output_csv = '';
end

distanze = [];
for i = 1:length(files)
    dati = readmatrix(files{i});
    x = dati(:,1);
    y = dati(:,2);

    % Trova tutti i massimi locali
    ind_max = find([false; y(2:end-1) > y(1:end-2) & y(2:end-1) > y(3:end); false]);
    x_max = x(ind_max);
    y_max = y(ind_max);

    % Ordina i massimi per posizione x
    [x_max, sort_idx] = sort(x_max);
    y_max = y_max(sort_idx);

    % FILTRO SULLA Y: Considera solo massimi con intensità > 100
    validi = y_max > soglia_y;
    x_max_valid = x_max(validi);
    y_max_valid = y_max(validi);

    % Filtraggio: trova la "catena" di massimi principali
    selezionati_x = [];
    selezionati_y = [];
    idx = 1; % parti dal primo massimo locale valido
    while idx <= numel(x_max_valid)
        selezionati_x(end+1,1) = x_max_valid(idx);
        selezionati_y(end+1,1) = y_max_valid(idx);

        % Candidati nel range di distanza
        candidati = find(x_max_valid > x_max_valid(idx) + min_dist & x_max_valid < x_max_valid(idx) + max_dist);
        if isempty(candidati), break; end
        [~, ind_ymax] = max(y_max_valid(candidati));
        idx = candidati(ind_ymax); % prossimo massimo principale
    end

    % Calcola distanze tra i massimi principali
    d = [];
    if numel(selezionati_x) > 1
        d = diff(selezionati_x);
        distanze = [distanze; d(:)];
    end

    % Stampa informazioni a schermo
    fprintf('File: %s\n', files{i});
    fprintf('  Numero di massimi locali trovati: %d\n', numel(ind_max));
    fprintf('  Numero di massimi validi (y > %d): %d\n', soglia_y, numel(x_max_valid));
    fprintf('  Numero di massimi principali selezionati: %d\n', numel(selezionati_x));
    fprintf('  Numero di distanze misurate: %d\n', numel(d));
    fprintf('------------------------------\n');

    % Plotta: tutti i massimi locali (blu), principali (rosso)
    figure;
    plot(x, y, 'b-'); hold on;
    plot(x_max, y_max, 'bo', 'MarkerSize', 6, 'LineWidth', 1.2); % tutti i massimi locali
    plot(selezionati_x, selezionati_y, 'ro', 'MarkerSize', 12, 'LineWidth', 2); % principali
    title(['File: ', strrep(files{i}, '_', '\_')], 'Interpreter', 'tex', 'FontWeight', 'bold');
    xlabel('Distanza [mm]');
    ylabel('Intensità');
    legend('Dati', 'Massimi locali', 'Massimi principali');
    hold off;
    pause(0.2);
end

if ~isempty(output_csv)
    writematrix(distanze, output_csv);
end
end