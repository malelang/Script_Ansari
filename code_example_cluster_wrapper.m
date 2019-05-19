diary on

% If not inside a PBS job, use 4 processors
if isempty(getenv('PBS_NP'))
    %NP = 4;
	display('You have not setup multiple core access!');
else
    NP = str2double(getenv('PBS_NP'));
end

% Initialize the pool to use
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    myPool = parpool('current', NP-2);
else
%     delete(p);
%     p = parpool(7);
end  


tic
% Call your function here

%peak_tester_VT

%Phase_Wrapper_batch('VT');
%DL_model_builder('VT');
%Phase_Wrapper_batch('VF');
DL_model_builder('VF');

stop_time = toc;
fprintf('Processing time with %d worker: %f seconds', NP, stop_time);

% Shut down the parallel pool
delete(myPool);


diary off
