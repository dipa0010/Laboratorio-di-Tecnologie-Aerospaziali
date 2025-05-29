close all
clear
clc
%% Importazione Dati e Conversioni 
% Nomi dei file
files = { '40_1.dat','40_2.dat','40_3.dat', ...
          'c70_1.dat','c70_2.dat','c70_3.dat' };

nProve = numel(files);

% Struttura per memorizzare tempo, spostamento e forza
dati = struct('tempo',[],'spostamento',[],'forza',[]);

for k = 1:nProve
    nomeFile = files{k};
    
    % Opzioni di importazione
    opzioni = detectImportOptions(nomeFile,'FileType','text');
    opzioni.DataLines = [9, Inf];        % i dati iniziano dalla riga 9
    opzioni.VariableNamesLine = 8;       % nomi variabili in riga 8
    opzioni.SelectedVariableNames = opzioni.VariableNames(1:3);
    opzioni = setvartype(opzioni,opzioni.SelectedVariableNames,'double');
    
    % Lettura tabella
    T = readtable(nomeFile, opzioni);
    
    % Memorizza in struttura (converti forza da kN a N)
    dati(k).tempo        = T{:,1};        
    dati(k).spostamento  = T{:,2};        
    dati(k).forza        = T{:,3} * 1e3;  
end


% Conversione in Sforzo e Deformazione
% Barrette (flessione): lunghezza, larghezza, altezza in mm
l_bar = [33.46, 33.52, 33.39]*1e-3;  % [m]
b_bar = [4.69, 4.60, 4.60]*1e-3;     % [m]
h_bar  = [2.27, 2.35, 2.24]*1e-3;     % [m]

% Cilindri (compressione): altezza e diametro medio in mm
h_cil = [ 12.81, 12.68, 12.75 ] *1e-3;    % [m]
d_cil = [  6.89,  6.96,  6.87 ] *1e-3;     % [m]
infill = 0.70; % infill nei cilindri al 70% 

% --- 3) Conversione in deformazione e sforzo per ciascun provino ---
for k = 1:3
    % flessione a tre punti
    F = dati(k).forza;
    D = dati(k).spostamento * 1e-3;  % mm→m
    L = 23.5*1e-3; % support span 
    b = b_bar(k);
    h = h_bar(k);
    
    % ε = 6·d·s / L^2
    dati(k).deformazione = (6 * h * D) / (L^2);
    % σ = 3·F·L / (2·b·h^2)
    dati(k).sforzo      = 3 .* F .* L ./ (2 .* b .* h.^2);
end

for k = 4:6
    idx = k-3;
    F = dati(k).forza;
    D = dati(k).spostamento * 1e-3;  % mm→m
    h = h_cil(idx);
    d = d_cil(idx);
    A = infill * pi*(d/2).^2;  
    
    dati(k).deformazione = D / h;
    dati(k).sforzo      = F ./ A;
end

% A questo punto:
% dati(i).tempo, dati(i).deformazione, dati(i).sforzo
% contengono le curve σ–ε per ciascuna delle 6 prove.


%% Plot Forza–Spostamento per tutte le prove (3 flessione + 3 compressione)
figure('Name','Curve Forza–Spostamento – Tutte le prove','Color','w','Units','normalized','Position',[.1 .1 .8 .6]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

for k = 1:6
    nexttile;
    D = dati(k).spostamento;  % spostamento [mm]
    F = dati(k).forza;        % forza [N]
    plot(D, F, 'LineWidth',1.5);
    xlabel('Spostamento (mm)');
    ylabel('Forza (N)');
    if k <= 3
        title(sprintf('Flessione Prova %d', k));
    else
        title(sprintf('Compressione Prova %d', k-3));
    end
    grid on;
end




%% ---------------- FLESSIONE ---------------- %% 

%% Selezione manuale della Regione Lineare, calcolo E e Statistiche
% a volte potrebbe buggarsi e selezionare da solo l'inizio della zona
% lineare (inserisce in automatico la linea verde), basta stoppare il
% codice e runnarlo di nuovo

E_flex = nan(1,3);
selFlex = nan(3,2);  % [eps_start, eps_end] per ciascuna prova

figure('Name','Fless. – Seleziona \epsilon_{start} e \epsilon_{end}','NumberTitle','off');
for k = 1:3
    eps = dati(k).deformazione;   
    sig = dati(k).sforzo;
    
    subplot(3,1,k);
    plot(eps, sig/1e6, '.-','LineWidth',1.2);  
    % qui sigma in MPa
    xlabel('\epsilon');
    ylabel('\sigma (MPa)');
    ylim([-120 10]);
    title(sprintf('Flessione %d: seleziona inizio e fine della parte lineare',k));
    grid on; hold on;
    
    % selezioni manuali
    [x1,~] = ginput(1);
    plot([x1 x1], ylim, 'g--','LineWidth',1.5);
    [x2,~] = ginput(1);
    plot([x2 x2], ylim, 'r--','LineWidth',1.5);
    
    selFlex(k,:) = [x1 x2];  % memorizzo deformazioni di inizio e fine
    
    % indici dei punti cliccati
    [~,i1] = min(abs(eps - x1));
    [~,i2] = min(abs(eps - x2));
    idx1 = min(i1,i2);
    idx2 = max(i1,i2);
    
    % estraggo tratto lineare e faccio fit
    x_lin = eps(idx1:idx2);
    y_lin = sig(idx1:idx2);
    p = polyfit(x_lin, y_lin, 1);
    E_flex(k) = p(1)/1e9;  % p(1) è la pendenza = modulo elastico (Pa)
    
    plot(x_lin, polyval(p,x_lin)/1e6,'m-','LineWidth',2);
    legend('dati','\epsilon_{start}','\epsilon_{end}','fit','Location','best');
end

% Calcolo statistiche
meanE_f = mean(E_flex);
stdE_f  = std(E_flex,1);
cvE_f   = 100*stdE_f/meanE_f;

%% Print risultati
fprintf('\n ------------ PROVE DI FLESSIONE ------------ \n\n');

fprintf('\n --- Calcolo E + Statistiche --- \n\n');
fprintf('  Prova   |   E [GPa]\n');
fprintf('---------------------\n');
for i = 1:3
    fprintf('   %2d     |   %8.2f\n', i, E_flex(i));
end
fprintf('---------------------\n');
fprintf('\nMedia       : %.2f MPa\n', meanE_f);
fprintf('Dev. std    :  %.2f MPa\n', stdE_f);
fprintf('CV(%%)       :  %.1f %%\n\n', cvE_f);
fprintf('\n');

%% Resistenza a Flessione

% (3-point flexure)
L_span = (23.5)*10^(-3);              % span di prova [m]

nF = 3;                      % numero di prove di flessione
R_f = nan(1,nF);

% 1) Estrazione del picco di forza e calcolo di σ_f
for k = 1:nF
    Pmax = max(abs(dati(k).forza));  % forza di rottura [N]
    % Legge la formula ASTM C1161-18 per 3-point flexure:
    % σ_f = 3·P_max·L_span / (2·b·h^2)
    R_f(k) = 3 * Pmax * L_span / (2 * b_bar(k) * h_bar(k)^2);
end

% 2) Statistiche
mean_Rf = mean(R_f);
std_Rf  = std(R_f,1);      % deviazione standard (normalizzazione n-1)
cv_Rf   = 100 * std_Rf/mean_Rf;

% 3) Output a video
fprintf('\n--- Resistenza a flessione (MPa) ---\n');
for k = 1:nF
    fprintf('Prova %d: σ_f = %.2f MPa\n', k, R_f(k)/1e6);
end
fprintf('\nMedia       : %.2f MPa\n', mean_Rf/1e6);
fprintf('Dev. std    : %.2f MPa\n', std_Rf/1e6);
fprintf('Coeff.var.  : %.1f %%\n\n', cv_Rf);



%% ---------------- COMPRESSIONE ---------------- %% 

%% Selezione manuale della Regione Lineare, calcolo E, Severity e Statistiche
E_comp = nan(1,3);
S_cr   = nan(1,3);
selComp= nan(3,2);    % [eps_start, eps_end]

figure('Name','Compressione: seleziona inizio e fine della parte lineare','Color','w');
for idx = 1:3
    k   = idx+3;                      
    eps = abs(dati(k).deformazione);  
    sig = abs(dati(k).sforzo);        

    subplot(3,1,idx);
    plot(eps, sig/1e6, '.-','LineWidth',1.2); hold on; grid on;
    xlabel('|ε|'); ylabel('|σ| (MPa)');
    title(sprintf('Compressione %d: seleziona inizio e fine della parte lineare',idx));

    [x1,~] = ginput(1); plot([x1 x1], ylim,'g--','LineWidth',1.5);
    [x2,~] = ginput(1); plot([x2 x2], ylim,'r--','LineWidth',1.5);
    selComp(idx,:) = [x1 x2];

    [~,i1] = min(abs(eps - x1));
    [~,i2] = min(abs(eps - x2));
    i_start = min(i1,i2);
    i_end   = max(i1,i2);

    % fit compattazione
    p1 = polyfit(eps(1:i_start), sig(1:i_start), 1);
    E1 = p1(1);
    % fit elastico
    p2 = polyfit(eps(i_start:i_end), sig(i_start:i_end), 1);
    E2 = p2(1);

    E_comp(idx) = E2;
    S_cr(idx)   = 100*(1 - E1/E2);

    plot(eps(1:i_start),    polyval(p1,eps(1:i_start))/1e6,'b-','LineWidth',1.5);
    plot(eps(i_start:i_end), polyval(p2,eps(i_start:i_end))/1e6,'m-','LineWidth',1.5);
    legend('dati','start','end','fit comp.','fit elast.','Location','best');
end
          %%
E_GPa     = (E_comp/1e9);             % da Pa a GPa
                  
% Calcola statistiche su E e su Severity
meanE = mean(E_GPa);  stdE  = std(E_GPa,1);  cvE  = 100*stdE/meanE;
meanS = mean(S_cr); stdS = std(S_cr,1); cvS = 100*stdS/meanS;

%% Print risultati
fprintf('\n ------------ PROVE DI COMPRESSIONE ------------ \n\n');

fprintf('\n --- Calcolo E + Severity + Statistiche --- \n\n');
fprintf('  Prova   |   E [GPa]   | Severity [%%]\n');
fprintf('--------------------------------------\n');
for i = 1:3
    fprintf('   %2d     |   %8.2f   |    %6.1f\n', ...
        i, E_GPa(i), S_cr(i));
end
fprintf('--------------------------------------\n');
fprintf('   Mean   |   %8.2f   |    %6.1f\n', meanE, meanS);
fprintf('   Std    |   %8.2f   |    %6.1f\n', stdE, stdS);
fprintf('   CV(%%)  |   %8.1f   |    %6.1f\n', cvE, cvS);
fprintf('\n');


%% Resistenza a compressione
nC = 3;                       % numero di prove di compressione
R_comp = nan(1,nC);          % vettore per le resistenze ultimate

for idx = 1:nC
    k = 3 + idx;              % le prove di compressione sono 4,5,6 in 'dati'
    
    % prendo il vettore di sforzi assoluti (in Pa)
    sigma = abs(dati(k).sforzo);
    
    % resistenza a compressione = valore massimo di sigma
    R_comp(idx) = max(sigma);
end

% Conversione in MPa per miglior leggibilità
R_comp_MPa = R_comp/1e6;

% Statistiche
X_comp = R_comp_MPa;   % resistenze ultime in MPa
n_comp = numel(X_comp);
mu_comp   = mean(X_comp);
s_comp    = std(X_comp, 1); % normalizzazione n-1
CV_comp   = 100 * s_comp / mu_comp;


fprintf('\n--- Resistenza a compressione (Pa) ---\n');
for i = 1:nC
    fprintf('Prova %d: R_c = %.2f MPa\n', i, R_comp_MPa(i));
end
fprintf(' \nValore medio    : %.2f MPa\n', mu_comp);
fprintf('Dev. standard   : %.2f MPa\n', s_comp);
fprintf('Coeff. variazione: %.1f %%\n', CV_comp);


