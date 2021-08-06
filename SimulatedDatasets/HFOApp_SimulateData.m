%% Add noise to simulated data

% Directory of the simulated datasets HFOApp_SimulatedData.py
wkpath = '~/Downloads/HFOAppSimuData';

% Datasets
datasets = {'t000', 't001', 't002', 't003', 't004', 't005', 't006', 't007', 't008', 't009'};

% saving directory of noisy data
savedir = fullfile( wkpath, 'SimulatedDataWithNoise');


cd( wkpath);


% fields required by HFOApp
start_clock = datetime( '1999-01-12 00:00:00');
labels = {'a'};
file_p = wkpath;

% signal-to-noise ratio levels
snr_levels = 1 : 10;

if ~exist( savedir, 'dir')
    mkdir( savedir);
end

nb_levels = length( snr_levels);
for k = 1 : nb_levels
    data = load( datasets{ k});
    
    srate = double( data.srate);
    N = length( data.mat);
    
    fact = 0.01 : 0.0001 : 2;
    s = zeros( length( fact), 1);
    for ind = 1 : length( fact)
        fprintf( '%d/%d\n', ind,  length( fact));
        s( ind) = snr( mat, fact( ind) * randn( 1, N));
    end
            
    final_snr = zeros( nb_levels, 1);
    final_fact = zeros( nb_levels, 1);
    for ind = 1 : nb_levels
        [mval, midx] = min( abs( s - snr_levels( ind)));
        final_snr( ind, 1) = s( midx);
        final_fact( ind, 1) = fact( midx);
        
        if ind < 10
            savename = [datasets{k}, '_SNR_0', num2str( ind), '.mat'];
        else
            savename = [datasets{k}, '_SNR_', num2str( ind), '.mat'];
        end 
        
        mat = data.mat + randn( 1, N) * fact( midx);
        file = fullfile( file_p, savename);
        save( fullfile( savedir, savename), 'mat', 'srate', 'labels', 'file', 'start_clock');        
    end
    
    save( fullfile( savedir, ['SNR_factor_', datasets{ k}]), 'final_snr', 'final_fact', 'snr_levels');
    
    mat = data.mat;
    savename = [datasets{k}, '_SNR_00.mat'];
    file = fullfile( file_p, savename);
    save( fullfile( savedir, savename), 'mat', 'srate', 'labels', 'file', 'start_clock');   
end


