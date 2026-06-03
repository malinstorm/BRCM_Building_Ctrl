%% Parametre
area     = "NO1";                      % Oslo / Øst-Norge
baseUrl  = "https://www.hvakosterstrommen.no/api/v1/prices"; % :contentReference[oaicite:0]{index=0}
endDate  = datetime(2025,12,10);       % velg sluttdato
startDate = endDate - days(29);        % 30 dager totalt

dates = startDate:endDate;

allDatetime = datetime.empty(0,1);
allPrice    = [];

%% Hent priser dag for dag
for d = dates
    url = sprintf('%s/%04d/%02d-%02d_%s.json', ...
                  baseUrl, year(d), month(d), day(d), area);

    dayData = webread(url);            % ?????????????????????????

    % JSON har 24 elementer per dag, med feltet NOK_per_kWh :contentReference[oaicite:1]{index=1}
    dayPrices = [dayData.NOK_per_kWh]; % 1x24

    % Lag tidspunkt 00:00, 01:00, ..., 23:00 for denne dagen
    dayTimes = d + hours(0:23);

    allDatetime = [allDatetime; dayTimes(:)];
    allPrice    = [allPrice; dayPrices(:)];
end

%% Lag tabell klar for Excel
T = table(allDatetime, allPrice, ...
    'VariableNames', {'Tidspunkt', 'Pris_NOK_per_kWh'});

%% Skriv til Excel-fil
filename = 'strompriser_30d_NO1.xlsx';
writetable(T, filename);

disp(['Skrev fil: ' filename]);
