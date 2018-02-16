function [out,info, iteration_count, time]...
    = cminpack(data, initial_parameters, functionID, tolerance, x_data)

versionID = 6;

if nargin < 5
    x_data = [];
end

%% configure parameters
if functionID == 0
    n_parameters = 4;
elseif functionID == 1
    n_parameters = 5;
elseif functionID == 2
    n_parameters = 6;
elseif functionID == 6
    n_parameters = 3;
elseif functionID == 7
    n_parameters = 4;
elseif functionID == 8
    n_parameters = 8;
elseif functionID == 9
    n_parameters = 9;
end

data_size = size(data);

fit_size = data_size(1);
if (ndims(data) == 2)
    n_fits = data_size(2);
else 
    n_fits = 1;
end

if ~isa(data, 'double')
    data = cast(data,'double');
end

if ~isa(initial_parameters, 'double')
    initial_parameters = cast(initial_parameters,'double');
end

%% run cminpack
tic;
[out, info, iteration_count]...
    = cminpackMex(versionID, data, fit_size, n_fits, n_parameters, initial_parameters, functionID, tolerance, x_data);
time = toc;

out = reshape(out, n_parameters, n_fits);

end
