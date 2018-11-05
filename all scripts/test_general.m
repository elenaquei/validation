% test function for all the examples for Saddle and Hopf bifurcation
global rescaling_saddle
global refinement_saddle
global talkative
talkative =0;

addpath(genpath('./'));

% as of 6th July
% 
% lor_cont: 1
% va der pol: 1
% Hopf: 1
% Rychkov: 1
% Lorenz: 1
% Saddle: 1
% Rossler: 0
% Hyper: 0

%possibilities = {'lor_cont','vdp','rychkov','lorenz','hopf','saddle','rossler','hyper','all'};
possibility = 'all';%possibilities(end);
profile on

if strcmp(possibility,'lor_cont')||strcmp(possibility,'all')
    try
        template_lor_con
        try_lor_cont=1;
    catch e
        fprintf(1,'There was an error! The message was:\n%s',e.message);
        error_lor_cont = e.message;
        try_lor_cont=0;
    end
    fprintf('Lor cont:%i\n',try_lor_cont);
end

if strcmp(possibility,'vdp')||strcmp(possibility,'all')
    try
        main_VanDerPol
        try_vdp=1;
    catch e
        fprintf(1,'There was an error! The message was:\n%s',e.message);
        error_vdp = e.message;
        try_vdp=0;
    end
    fprintf('Van Der Pol:%i\n',try_vdp);
end

if strcmp(possibility,'hopf')||strcmp(possibility,'all')
    try
        example_hopf % works as of 5th July 2018
        try_example_hopf=1;
    catch e
        fprintf(1,'There was an error! The message was:\n%s',e.message);
        error_hopf = e.message;
        try_example_hopf=0;
    end
    fprintf('Hopf:%i\n',try_example_hopf);
end
if strcmp(possibility,'rychkov')||strcmp(possibility,'all')
    rescaling_saddle = 500;
    refinement_saddle =300;
    try
        main_real_rychkov % works as of 5th July 2018
        try_real_rychkov=1;
    catch e
        fprintf(1,'There was an error! The message was:\n%s',e.message);
        error_real_rychkov = e.message;
        try_real_rychkov=0;
    end
    fprintf('Rychkov:%i\n',try_real_rychkov);
end
if 1==2 && strcmp(possibility,'lorenz')||strcmp(possibility,'all')
    script_algebraic_hopf;
    rescaling_saddle = 500;
    refinement_saddle =300;
    try
        lorenz84_validated_cont; % works as of 5th July 2018
        try_lorenz84=1;
    catch e
        fprintf(1,'There was an error! The message was:\n%s',e.message);
        error_lorenz84 = e.message;
        try_lorenz84=0;
    end
    fprintf('Lorenz:%i\n',try_lorenz84);
end
if strcmp(possibility,'saddle')||strcmp(possibility,'all')
    rescaling_saddle=1;
    refinement_saddle =300;
    try
        main_saddle; % works as of 5th July 2018
        try_saddle=1;
    catch e
        fprintf(1,'There was an error! The message was:\n%s',e.message);
        error_saddle = e.message;
        try_saddle=0;
    end
    fprintf('Saddle:%i\n',try_saddle);
end
if strcmp(possibility,'rossler')||strcmp(possibility,'all')
    rescaling_saddle=72;%??
    refinement_saddle =72;%??
    try
        main_rossler; % -- DOESN'T WORK as of 5th July 2018
        try_rossler=1;
    catch e
        fprintf(1,'There was an error! The message was:\n%s',e.message);
        error_rossler = e.message;
        try_rossler=0;
    end
    fprintf('Rossler:%i\n',try_rossler);
end
if strcmp(possibility,'hyper')||strcmp(possibility,'all')
    rescaling_saddle=22;%??
    refinement_saddle =22;%??
    try
        main_hyper; % -- DOESN'T WORK as of 5th July 2018
        try_hyper=1;
    catch e
        fprintf(1,'There was an error! The message was:\n%s',e.message);
        error_hyper = e.message;
        try_hyper=0;
    end
    fprintf('Hyper:%i\n',try_hyper);
end

profile viewer
profile off


fprintf('\n\nRESULTS\n\n\n')
if strcmp(possibility,'lor_cont')||strcmp(possibility,'all')
    fprintf('Lor Cont:%i\n',try_lor_cont);
    if ~try_lor_cont
        fprintf(error_lor_cont);
    end
end
if strcmp(possibility,'vdp')||strcmp(possibility,'all')
    fprintf('Van der Pol:%i\n',try_vdp);
    if ~try_vdp
        fprintf(error_vdp);
    end
end
if strcmp(possibility,'hopf')||strcmp(possibility,'all')
    fprintf('Hopf:%i\n',try_example_hopf);
    if ~try_example_hopf
        fprintf(error_hopf);
    end
end
if strcmp(possibility,'rychkov')||strcmp(possibility,'all')
    fprintf('Rychkov:%i\n',try_real_rychkov);
    if ~try_real_rychkov
        fprintf(error_real_rychkov);
    end
end
if strcmp(possibility,'lorenz')||strcmp(possibility,'all')
    fprintf('Lorenz:%i\n',try_lorenz84);
    if ~try_lorenz84
        fprintf(error_lorenz84);
    end
end
if strcmp(possibility,'saddle')||strcmp(possibility,'all')
    fprintf('Saddle:%i\n',try_saddle);
    if ~try_saddle
        fprintf(error_saddle);
    end
end
if strcmp(possibility,'rossler')||strcmp(possibility,'all')
    if ~try_rossler
        fprintf(error_rossler);
    end
    fprintf('Rossler:%i\n',try_rossler);
end
if strcmp(possibility,'hyper')||strcmp(possibility,'all')
    if ~try_hyper
        fprintf(error_hyper);
    end
    fprintf('Hyper:%i\n',try_hyper);
end