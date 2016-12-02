%
% Bottomside IRI parameters 
%

graph = 15;

trange = datenum(2002,6,[11 12],[18 5],[0 30],[0 59]);

% Initial guess :
%       pars0(1) = Nm (cm-3)
%       pars0(2) = hm (km)
%       pars0(3) = H (km)
%       pars0(4) = B0 (km)
%       pars0(5) = B1 (non-dimensional)
%
pars0 = [5.5e5 220.5 50.0 50.0 3.0];

% Factor of Nm which indicates the minimum F-region height
fact_NmF2B = 0.4;

% Minimum height of data analyzed (in km)
min_height = 120.0;

root_dir = '/home';
%myFontName = 'New Century Schoolbook';
myFontName = 'new century schoolbook';
myFontSize = 14;
myMarkerSize = 6;
%gpath = [root_dir '/rilma/work/programs/matlab/cornell/eas5840/project/'];
gpath = [root_dir '/rilma/tmp/results/IRIBottom/'];
        
setvalues = struct('GPath', gpath, 'FactNm', fact_NmF2B, ...
    'FontName', myFontName, 'FontSize', myFontSize, ...
    'MarkerSize', myMarkerSize, 'MinHeight', min_height);    

    
% Read JRO ISR Oblique data
%
root_dir = '/home';

dpath = [root_dir '/rilma/database/jro/isroblique/madrigal/'];
filename = [repmat(dpath, 1, 1), 'jro020611c.001'];        

data = read_ncar_jro(filename,[]);         

isr_time = data.time + data.timezone/24;
ind = find(isr_time >= trange(1) & isr_time <= trange(2));

jro_data = struct('time', isr_time(ind), ...
    'height', data.data.range(1, :), ...
    'un_ne', data.data.un_ne(ind, :), ...
    'un_ne_err', data.data.un_ne_err(ind, :));

    nprof = numel(jro_data.time);

    pars = NaN(5, nprof); iterations = NaN(2, nprof);

    for i = 1 : nprof

        isr_oblique = struct('time', jro_data.time(i), ...
            'height', jro_data.height, ...
            'edens', jro_data.un_ne(i, :), ...
            'err_edens', jro_data.un_ne_err(i, :));

        disp(['Profile ' num2str(i, '%03i') ...
            '.   Time: ' datestr(isr_oblique.time, 0)]);

        [yy,mm,dd,hh,mi,ss] = datevec(isr_oblique.time);
        hhFloat = hh + mi / 60.0 + ss / 3600.0;

        if hhFloat <= 8.75, hrange = [200 550]; end

        if hhFloat >= 8.75 && hhFloat < 22.0, hrange = [250 700]; end

        if hhFloat >= 22.0, hrange = [200 550]; end            

        if i > 1 && i < 3, pars0 = pars(:, i - 1)'; end;
        if i >= 3, pars0 = median(pars(:, 1 : (i - 1)), 2)'; end;
        
        result = calc_profpars(isr_oblique, graph, hrange, pars0, setvalues);

        pars(:, i) = result.m;
        iterations(:, i) = result.nloop;

    end

        
