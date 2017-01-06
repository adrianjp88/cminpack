function [] = test_cminpack()
fit_size = 5;
functionID = 1; %GAUSS_2D
n_parameters = 5;
tolerance = 0.001;
snr = 10;
noise = 'gauss';

%% set start values
initial_parameters = ones(1,n_parameters);
initial_parameters(1) = 500;
initial_parameters(2) = (fit_size-1)/2;
initial_parameters(3) = (fit_size-1)/2;
initial_parameters(4) = 1;
initial_parameters(5) = 10;

%% test case
LogNFitsMin = 0;
LogNFitsMax = 5;
sampling_factor = 10;
ranges = logspace(LogNFitsMin,LogNFitsMax,LogNFitsMax-LogNFitsMin+1);
temp = zeros(LogNFitsMax-LogNFitsMin,10/sampling_factor);
stepslog = LogNFitsMin:LogNFitsMax;
for index = 1:length(stepslog)-1
    steps = 10^(stepslog(index)) * sampling_factor;
    temp(index,:) = steps:steps:ranges(index+1);
end
n_fits = reshape(temp', [1 10/sampling_factor*(LogNFitsMax-LogNFitsMin)]); 
n_fits = [10^LogNFitsMin n_fits];

%% noise
mean_amplitude = 2 * pi * initial_parameters(1) * initial_parameters(4) * initial_parameters(4)/fit_size/fit_size;
noise_std_dev = mean_amplitude / snr;

%% test loop
for ifit = 1:length(n_fits)
    
    initial_parameters_set = repmat(initial_parameters, 1, n_fits(ifit));

    fprintf('%d fits\n', n_fits(ifit));

    %% generate data
    xpos_mean = initial_parameters(2);
    xpos_offset = rand(n_fits(ifit), 1) - 0.5;
    input_xpos = xpos_mean + xpos_offset;

    ypos_mean = initial_parameters(3);
    ypos_offset = rand(n_fits(ifit), 1) - 0.5;
    input_ypos = ypos_mean + ypos_offset;

    parameters.a = initial_parameters(1);
    parameters.x0 = input_xpos;
    parameters.y0 = input_ypos;
    parameters.sx = initial_parameters(4);
    parameters.sy = initial_parameters(4);
    parameters.b = initial_parameters(5);
    data = generate_2Dgaussians(parameters, n_fits(ifit), fit_size);
    data = data + noise_std_dev * randn(fit_size,fit_size,n_fits(ifit));
    data = permute(data,[2,1,3]);

    %% run cminpack
    [...
        parameters_cminpack,...
        info_cminpack,...
        n_iterations_cminpack,...
        time_cminpack]...
    = cminpack(data, initial_parameters_set, functionID, tolerance);

    converged_cminpack = (info_cminpack > 0) & (info_cminpack <= 3);
    calculated.a = parameters_cminpack(1:n_parameters:end).';
    calculated.x0 = parameters_cminpack(2:n_parameters:end).';
    calculated.y0 = parameters_cminpack(3:n_parameters:end).';
    calculated.s = parameters_cminpack(4:n_parameters:end).';
    calculated.b = parameters_cminpack(5:n_parameters:end).';
    
    speed_cminpack(ifit) = n_fits(ifit)/time_cminpack;
    precision_cminpack(ifit) = calculate_precision(calculated, parameters, converged_cminpack);

    print_fit_info(precision_cminpack(ifit), time_cminpack, 'cminpack', converged_cminpack, n_fits(ifit), n_iterations_cminpack);
      
end

%% plot
figure('Name', 'cminpack', 'NumberTitle', 'off');
semilogx(...
    n_fits, speed_cminpack, 'blue.-', ...
    'LineWidth', 8)
xlabel('number of fits')
ylabel('fits per second')
legend('cminpack')

grid on;
box off;
current_figure = gca;
current_figure.FontSize = 37;
current_figure.LineWidth = 4;

end

function [data] = generate_2Dgaussians(parameters, count, size)

fprintf('generating data... ');
[xi,yj] = meshgrid(0:size-1,0:size-1);

xi = repmat(xi, [1 1 1]);
yj = repmat(yj, [1 1 1]);

data = zeros(size, size, count, 'single');

for i = 1:count
    data(:,:,i) = parameters.a...
        *exp(-1/2*((xi-parameters.x0(i))/parameters.sx).^2)...
        .*exp(-1/2*((yj-parameters.y0(i))/parameters.sy).^2)...
        + parameters.b;
end

fprintf('done!\n');

end

function [precision] = calculate_precision(calculated, parameters, converged)

[calculated, parameters] = delete_failed(calculated, parameters, converged);

precision.a = std(calculated.a-parameters.a) / parameters.a;
precision.x0 = std(calculated.x0-parameters.x0);
precision.y0 = std(calculated.y0-parameters.y0);
precision.s = std(calculated.s-parameters.sx) / parameters.sx;
precision.b = std(calculated.b-parameters.b) / parameters.b;

end

function [calculated, parameters] = delete_failed(calculated, parameters, converged)

not_converged = find(converged~=1);

calculated.a(not_converged) = [];
calculated.x0(not_converged) = [];
calculated.y0(not_converged) = [];
calculated.s(not_converged) = [];
calculated.b(not_converged) = [];

parameters.x0(not_converged) = [];
parameters.y0(not_converged) = [];

end

function [] = print_fit_info(precision, time, version, converged, n_fits, iterations)

fprintf('***%s***\t', version);
fprintf('x precision: %g\t|\t', precision.x0);
fprintf('time: %g\t|\t', time);
fprintf('converged: %f \t|\t', double(size(find(converged == 1),2)) / double(n_fits));
fprintf('iterations: %g\n', sum(iterations)/n_fits);

end
