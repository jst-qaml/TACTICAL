classdef TACTICAL < handle
    properties
        %% model info
        bm              % benchmark
        mdl             % the model with problematic NN controller
        mdl_m           % the model for fault injection by mutating weights
        mdl_b           % the buggy model (with artificial bugs), for repair
        D_run           % dataset for simulation, including ps_x and ps_y
        is_nor          % is normalized or not
        D               % dataset with test suite, for fault localization
        D_size          % the size of dataset
        T
        Ts
        %% nn parameters
        net             % neural network controller
        nn_stru         % nn structure (excludes the output layer)
        layer_num       % the number of hidden layers (does not include output layer)
        traj_sz
        weight          % weight of hidden layers and output layer stored as a cell array according to the layer index
        bias            % bias of hidden layers and output layer stored as a cell array according to the layer index
        weight_list
        max_weight      % max weight value
        min_weight      % min weight value
        %% config of input signal of model
        in_name         % input signal of the model
        in_range        % the range of input signal
        in_span         % the time span of input signal
        %% config of input signal of controller
        icc_name        % input constant of controller
        ics_name        % input signal of controller
        ic_const        % input constant value of controller
        %% config of output signal of controller
        oc_name         % output signal of controller
        oc_span         % the time span of output signal
        %% specification
        phi_str
        phi
        sig_str
        %% the parameters of differential evolution algorithm
        dim                     % individual dimension
        %% parallel computing
        core_num
    end
    methods
        function this = TACTICAL(bm, mdl, mdl_m, D_run, is_nor, D, D_size, net, nn_stru, T, Ts, in_name, in_range, in_span, icc_name, ic_const, ics_name, oc_name, oc_span, phi_str, sig_str, core_num)

            this.bm = bm;
            this.mdl = mdl;
            this.mdl_m = mdl_m;
            this.D_run = D_run;
            this.is_nor = is_nor;
            this.D = D;
            this.D_size = D_size;

            this.T = T;
            this.Ts = Ts;

            this.net = net;
            this.nn_stru = nn_stru;
            if numel(nn_stru) ~= net.numLayers - 1
                error('the layer num described in nn_stru and the layer num stored in nn file are inconsistent!');
            end
            this.layer_num = net.numLayers - 1;
            this.traj_sz = T/Ts + 1;
            % assign weight and bias
            this.weight = cell(1, this.layer_num);
            this.bias = cell(1, this.layer_num);
            for li = 1:this.layer_num 
                if li == 1
                    this.weight{1,li} = this.net.IW{1,1};
                    [li_row, ~] = size(this.net.IW{1,1});
                    if nn_stru(1,li) ~= li_row
                        error('the structure described in nn_stru and the structure stored in nn file are inconsistent!');
                    end
                else
                    this.weight{1,li} = this.net.LW{li,li-1};
                    [li_row, ~] = size(this.net.LW{li,li-1});
                    if nn_stru(1,li) ~= li_row
                        error('the structure described in nn_stru and the structure stored in nn file are inconsistent!');
                    end
                end
                this.bias{1,li} = this.net.b{li};
            end

            this.in_name = in_name;
            this.in_range = in_range;
            this.in_span = cell2mat(in_span);

            this.icc_name = icc_name;
            this.ic_const = ic_const;
            this.ics_name = ics_name;

            this.oc_name = oc_name;
            this.oc_span = cell2mat(oc_span);

            this.phi_str = phi_str;
            this.phi = STL_Formula('phi', this.phi_str);
            this.sig_str = sig_str;

            % obtain the range of NN controller' weights
            weight_cp = this.weight;
  
            this.weight_list = [];
            for li = 2: this.layer_num
                this.weight_list = [this.weight_list; weight_cp{1,li}(:)];
            end

            this.min_weight = min(this.weight_list);
            this.max_weight = max(this.weight_list);

            this.core_num = core_num;
        end
        %% functions
        function [rob, Br, tau_s, ic_sig_val, oc_sig_val] = signalDiagnosis(this, mdl, in_sig, spec_i)
            % signalDiagnosis function returns the robustness, Br, tau_s, controller's input signal and controller's output signal of a system execution
            %
            % Inputs:
            %   mdl: simulink model
            %   in_sig: the input signals of the simulink model
            %   spec_i: specification index
            % Outputs:
            %   rob: robustness
            %   Br: BreachSimulinkSystem of the Simulink model
            %   tau_s: the first timestamp at which the robustness value turns negative
            %   ic_sig_val: input signal of controller
            %   oc_sig_val: output signal of controller

            Br = BreachSimulinkSystem(mdl);
            Br.Sys.tspan = 0:this.Ts:this.T;
            % in_span is a scalar value or not
            if isscalar(unique(this.in_span))
                con_type = 'UniStep';
                input_gen.cp = this.T/this.in_span(1,1);
            else
                con_type = 'VarStep';
                input_gen.cp = this.T./this.in_span;
            end

            input_gen.type = con_type;
            Br.SetInputGen(input_gen);

            if strcmp(con_type, 'UniStep')
                for i = 1:numel(this.in_name)
                    for cpi = 0:input_gen.cp - 1
                        eval(['Br.SetParam({''',this.in_name{1,i},'_u',num2str(cpi),'''}, in_sig{1,i}(1,cpi+1));']);
                    end
                end
            elseif strcmp(con_type, 'VarStep')
                for i = 1:numel(this.in_name)
                    for cpi = 0: input_gen.cp(1,i) - 1
                        eval(['Br.SetParam({''',this.in_name{1,i},'_u',num2str(cpi),'''}, in_sig{1,i}(1,cpi+1));']);
                        if cpi ~= input_gen.cp(1,i) - 1
                            eval(['Br.SetParam({''',this.in_name{1,i},'_dt',num2str(cpi),'''},this.in_span(1,i));']);
                        end
                    end
                end
            end

            Br.Sim(0:this.Ts:this.T);
            rob = Br.CheckSpec(this.phi);

            % extract ic signal and oc signal for forward impact analysis
            % (analysis the impact of a weight to the final result)
            ic_sig_val = sigMatch(Br, this.ics_name);
            oc_sig_val = sigMatch(Br, this.oc_name);

            if this.is_nor == 1
                x_gain = this.D_run.ps_x.gain;
                x_gain = diag(x_gain);
                x_offset = this.D_run.ps_x.xoffset;
                ic_sig_val = ic_sig_val - x_offset;
                ic_sig_val = x_gain * ic_sig_val;

                y_gain = this.D_run.ps_y.gain;
                y_offset = this.D_run.ps_y.xoffset;
                oc_sig_val = oc_sig_val/y_gain;
                oc_sig_val = oc_sig_val + y_offset;
            end

            % calculate the tau_s according to the given specification
            if rob > 0
                tau_s = this.T;
                return;
            end

            if strcmp(this.bm, 'ACC')
                interval_LR = {[0,50]};
            elseif strcmp(this.bm, 'AFC') && spec_i == 1
                interval_LR = {[0,30]};
            elseif strcmp(this.bm, 'AFC') && spec_i == 2
                interval_LR = {[10,30]};
            elseif strcmp(this.bm, 'WT')
                interval_LR = {[4,4.9], [9,9.9], [14,14.9]};
            elseif strcmp(this.bm, 'SC')
                interval_LR = {[30,35]};
            end

            scan_interval = zeros(1, this.traj_sz);
            neg_interval = zeros(1, this.traj_sz);
            for scan_i = 1: numel(interval_LR)
                LR = interval_LR{1, scan_i};
                LR = LR./this.Ts;
                scan_interval(1, LR(1,1)+1:LR(1,2)+1) = 1;
            end

            if strcmp(this.bm, 'ACC') && spec_i == 1
                left_neg_interval = sigMatch(Br, "d_rel") - 1.4 * sigMatch(Br, "v_ego") < 10;
                right_neg_interval = sigMatch(Br, "v_ego") > 30.1;
                neg_idx = (left_neg_interval + right_neg_interval > 0);
                neg_interval(1, neg_idx) = 1;
            elseif strcmp(this.bm, 'ACC') && spec_i == 2
                delta = sigMatch(Br, "d_rel") - 1.4 * sigMatch(Br, "v_ego");
                % i == 500, i = 501, true
                for i = 1: 451
                    if delta(1, i) < 12
                        i_behind = i + 50;
                        delta_behind = delta(1, i:i_behind);
                        if ~any(delta_behind >= 12)
                            neg_interval(1, i_behind) = 1;
                        end
                    end
                end
            elseif strcmp(this.bm, 'AFC') && spec_i == 1
                af = sigMatch(Br, "AF");
                mu = abs(af - 14.7)/14.7;
                neg_interval(mu >= 0.2) = 1;
            elseif strcmp(this.bm, 'AFC') && spec_i == 2
                af = sigMatch(Br, "AF");
                % i == 300, i = 301, true
                for i = 101: 286
                    if af(1, i) >= 1.1*14.7 || af(1, i) <= 0.9*14.7
                        i_behind = i + 15;
                        af_behind = af(1, i:i_behind);
                        if ~any(af_behind <= 1.1*14.7 & af_behind >= 0.9 *14.7)
                            neg_interval(1, i_behind) = 1;
                        end
                    end
                end
                % localize tau_s using online monitoring. Thank you, Zhenya
                % tau = 0;
                % time_ = Br.P.traj{1}.time;
                % time = round(time_, 2);
                % sig_list = cellstr(split(this.sig_str, ","));
                % X = sigMatch(Br, sig_list);
                % trace = [time; X];
                % [rob_up, rob_low] = stl_eval_mex_pw(this.sig_str, this.phi_str, trace, tau);
                % rob_neg = find(rob_up < 0);
                % tau_s = (rob_neg(1, 1) - 1) * this.Ts;
                % % extract the input and output signals of nn controller in the interval [0, tau_s]
                % ic_sig_val = ic_sig_val(:, 1:rob_neg(1, 1));
                % oc_sig_val = oc_sig_val(:, 1:rob_neg(1, 1));
                % return;
            elseif strcmp(this.bm, 'WT')
                h_error = abs(sigMatch(Br, "h_error"));
                neg_interval(h_error > 0.86) = 1;
            elseif strcmp(this.bm, 'SC')
                pressure = sigMatch(Br, "pressure");
                neg_interval(pressure < 87 | pressure > 87.5) = 1;
            end
            int_interval = neg_interval .* scan_interval;
            if ~ismember(1,int_interval)
                sys_time = datevec(datestr(now));
                % record current bug
                buglog_filename = 'buglog';
                for i = 1:numel(sys_time)
                    buglog_filename = strcat(buglog_filename, '_', num2str(sys_time(1,i)));
                end
                buglog_filename = [buglog_filename, '.mat'];
                save(buglog_filename, 'Br', 'mdl', 'in_sig');
                disp('robustness > 0');
                rob = rob + 100000;
                tau_s = this.T;
                return;
            end
            negl_idx = find(int_interval == 1);
            %
            tau_s = (negl_idx(1,1) - 1) * this.Ts;
            % extract the input and output signals of nn controller in the interval [0, tau_s]
            ic_sig_val = ic_sig_val(:,1:negl_idx(1,1));
            oc_sig_val = oc_sig_val(:,1:negl_idx(1,1));

            % the following codes have bugs on ACC2_spec2_dataset1_index_8, WT
            % % localize tau_s using online monitoring. Thank you, Zhenya
            % tau = 0;
            % time_ = Br.P.traj{1}.time;
            % time = round(time_, 2);
            % sig_list = cellstr(split(this.sig_str, ","));
            % X = sigMatch(Br, sig_list);
            % trace = [time; X];
            % [rob_up, rob_low] = stl_eval_mex_pw(this.sig_str, this.phi_str, trace, tau);
            % rob_neg = find(rob_up < 0);
            % tau_s = (rob_neg(1, 1) - 1) * this.Ts;
            %
            % % extract the input and output signals of nn controller in the interval [0, tau_s]
            % ic_sig_val = ic_sig_val(:, 1:rob_neg(1, 1));
            % oc_sig_val = oc_sig_val(:, 1:rob_neg(1, 1));
        end

        function var = initializeCell(this, str)
            % initializeCell function can initialize the corresponding cell
            % according to str.
            %
            % Inputs:
            %   str: cell name
            % Outputs:
            %   intialized variable

            switch str
                case 'neuronOutSnapShot'
                    var = {};
                    for li = 1:this.layer_num
                        var{end+1} = zeros(this.nn_stru(1,li), 1);
                    end
                otherwise
                    error('Check cell name!');
            end
        end

        function [neuronOut] = calNeuronOut(this, input)
            % weightToNeuron function calculates the first part of the feed forward impact.
            % (a weight of a neuron to the output of this neuron)
            %
            % Inputs:
            %   input: input signal of nn at a timestamp (at a certain frame)
            % Outputs:
            %   weight2Neuron: the impact of a weight of a neuron to the output of this neuron during a system execution
            %   neuronOut: the output of each neuron during a system execution

            % calculate the frame number in a system execution (tau_s/Ts)
            [row, column] = size(input);
            % initialize the neuronOut to store neuronOutSnapShot
            neuronOut = {};

            for j = 1: column
                % record the output of each neuron at a certain frame
                neuronOutSnapShot = this.initializeCell('neuronOutSnapShot');
            
                for i = 1: this.layer_num
                    if i == 1
                        transFcn = this.net.layers{i}.transferFcn();
                        transFcn = str2func(transFcn);
                        neuronOutSnapShot{1,1} = transFcn(this.weight{1,1} * input(:,j) + this.bias{1,i});
                    else
                        transFcn = this.net.layers{i}.transferFcn();
                        transFcn = str2func(transFcn);
                        neuronOutSnapShot{1,i} = transFcn(this.weight{1,i} * neuronOutSnapShot{1,i-1} + this.bias{1,i});
                    end
                end
                neuronOut{end+1} = neuronOutSnapShot;
              
            end
        end

        function diagInfo = sysDiagnosis(this, mdl, in_sig, spec_i)
            % Given an external input signal of the system, forwardImpactAnalysis function returns the diagnostic information
            % of current system execution.
            % Inputs:
            %   mdl: simulink modelthe diagnostic information of the forward  impactof current system execution.
            %   in_sig: the input signals of the simulink model
            %   spec_i: specification index
            % Outputs:
            %   diagInfo: a struct store the diagnostic information of the forward impact of current system execution.
            %   diagInfo.Br: BreachSimulinkSystem of the Simulink model
            %   diagInfo.rob: robustness value
            %   diagInfo.state: the state of current system execution, passed or failed
            %   diagInfo.tau_s: the first timestamp at which the robustness value turns negative
            %   diagInfo.ic_sig_val: input signal of controller during [0, tau_s]
            %   diagInfo.oc_sig_val: output signal of controller during [0, tau_s]
            %   diagInfo.weight2Neuron: the impact of a weight of a neuron to the output of this neuron during a system execution
            %   diagInfo.neuronOut: the output of each neuron during a system execution
            %   diagInfo.neuron2Out: the gradient of the final output to the output of a neuron at a certain layer, i.e., the impact of the output of a neuron to the final output during a system execution
            %   diagInfo.weight2Out: the forward impact of each weight to the final output during a system execution

            % perform signal diagnosis
            [rob, Br, tau_s, ic_sig_val, oc_sig_val] = this.signalDiagnosis(mdl, in_sig, spec_i);
            % calculate the frame number in a system execution (tau_s/Ts)
            [~, frame_num] = size(ic_sig_val);
            % obtain the state of the system execution
            if rob >= 0
                state = 1;
            else
                state = 0;
            end
            % calculate the neuron output
            neuronOut = this.calNeuronOut(ic_sig_val);

            diagInfo.Br = Br;
            diagInfo.rob = rob;
            diagInfo.state = state;
            diagInfo.tau_s = tau_s;
            diagInfo.ic_sig_val = ic_sig_val;
            diagInfo.oc_sig_val = oc_sig_val;
            diagInfo.neuronOut = neuronOut;
        end

        function [bestInd, bestFitness, bestRobList, X_log, fitnessX_log, roblistX_log, mutant] = altWeightSearch(this, mutant, fl_mode, rep_mode, pert_level, DE_param, repLog_ID, core_num, res_name_prefix)
            % altWeightSearch function performs a differential evolution algorithm to search for a set of alternative weights.
            %
            % Inputs:
            %   mutant:
            %   fl_mode: 'tactical' or 'random'
            %   rep_mode: repair mode, i.e., the type of fitness function
            %   pert_level: perturbation level
            %   DE_param: the parameters of differential evoluation algorithm, including pop_size, max_gen, F, CR, mut_sgy, cr_sgy
            %   core_num: core number
            % The Unused Inputs:
            %   tol_rate: if y(1:subts_sz) < tol_rate * subts_sz, we discard current alt_weight. The default tol_rate is set as 0.7
            % Outputs:
            %   bestInd: the best individual of each generation
            %   bestFitness: the best fitness value of each generation
            %   bestRobList: the best robustness list of each generation
            %   X_log: record the population of each generation
            %   fitnessX_log: record the fitness of each generation
            %   roblistX_log: record the robustness list of each generation

            % there is a bug! if the safety rate of two individuals are the same, we need to select the one with a greater robustness value sum
            % reset this.dim
            % set a timer
            tic;

            if strcmp(fl_mode, 'tactical')
                sps_weight = mutant.sps_weight;
            elseif strcmp(fl_mode, 'random')
                % obtain all the weights
                weight_pool = [];
                for li = this.roil_l:this.roil_r
                    [rows, cols] = size(this.weight{1, li});
                    for i = 1:rows
                        for j = 1:cols
                            if ~ismember([li,i,j], mutant.mut_op(:,1:3),'rows')
                                weight_pool = [weight_pool; 0, this.weight{1,li}(i,j), li, i, j];
                            end
                        end
                    end
                end
                weight_num = size(weight_pool, 1);

                seli = randperm(weight_num, size(mutant.sps_weight,1));
                % randomly select sps_weight
                sps_weight = weight_pool(seli, :);
                % Overwrite the previous sps_weight stored in mutant
                mutant.sps_weight = sps_weight;
            else
                error('check fl_mode!');
            end

            if strcmp(DE_param.init_type, 'internal')
                this.dim = size(sps_weight,1);

                this.mdl_b = mutant.bug_mdl;
                % sort TS according to the robustness value of each test case
                sortedTS = sortTestSuite(mutant.mut_TS);

                % no change, no score, so the initial fitness value of the test suite is 0
                pre_fitnessbestX = 0;

                % obtain the lb and ub of the suspicious weights
                [lb, ub] = obtainBound(sps_weight);
                % initialization, X represents an individual, the column represents the dim of X
                % X = (ub - lb) .* rand(this.DE_param.pop_size, this.dim) + lb;
                % apply gaussian distribution to generate the initial population
                % In practice, the individuals generated by globalGauWeightInit underperform the original weight
                [X_g, ~] = globalGauWeightInit(this.weight, sps_weight, ceil(DE_param.pop_size/2));
                % generate
                X_l = localGauWeightInit(sps_weight, floor(DE_param.pop_size/2), pert_level);
                raw_X = [X_g; X_l];
                % shuffle the order of each row in raw_X
                rowrank = randperm(size(raw_X, 1));
                X = raw_X(rowrank,:);
            elseif strcmp(DE_param.init_type, 'external')
                % start repair from the break-point
                repLogDir = dir(fullfile(DE_param.init_folder,'*.mat'));
                repLogNames = {repLogDir.name};

                str1 = ['_K_', num2str(repLog_ID.K), '_'];
                str2 = ['_', repLog_ID.mut_op_mode, '_', num2str(repLog_ID.mut_op_noise), '_'];
                str3 = ['_repmode_', rep_mode, '_'];
                str4 = ['_mutant_', num2str(repLog_ID.mut_i), '_'];
                str5 = ['_budget_', num2str(repLog_ID.bi)];
                target_str = repLogNames(contains(repLogNames, str1) & contains(repLogNames, str2) & contains(repLogNames, str3) & contains(repLogNames, str4) & contains(repLogNames, str5));
                if numel(target_str) > 1 || numel(target_str) == 0
                    error('Check the DE_param.init_type and DE_param.init_folder!');
                end

                rep_log = load(target_str{1,1});
                this.dim = size(sps_weight,1);
                this.mdl_b = mutant.bug_mdl;
                % sort TS according to the robustness value of each test case
                sortedTS = sortTestSuite(mutant.mut_TS);
                % reset this.weight (has been corrected in mutation phase)
                [this.weight, ~, ~] = weightSync(this.net, mutant.mut_mode, mutant.mut_op);

                % no change, no score, so the initial fitness value of the test suite is 0
                pre_fitnessbestX = 0;

                % obtain the lb and ub of the suspicious weights
                [lb, ub] = obtainBound(sps_weight);
                nonempty_idx = find(~cellfun('isempty', rep_log.X_log));
                if ~isempty(nonempty_idx)
                    X = rep_log.X_log{1, nonempty_idx(end)};
                else
                    error('Check the rep_log!');
                end
            else
                error('Check the DE_param.init_type!');
            end

            % record the best individual and its fitness value of each generation
            % this.dim should be initialized in advanced! (the dim of the individual)
            bestInd = zeros(DE_param.max_gen, this.dim);
            bestFitness = zeros(DE_param.max_gen, 1);
            % record the robustness value list from the best individual
            bestRobList = cell(1, DE_param.max_gen);
            % record the fitness value of each individual X
            fitnessX = zeros(1, DE_param.pop_size);
            % store the fitness value of each individual U obtained by crossover
            fitnessU = zeros(1, DE_param.pop_size);
            % record the robustness value list of each individual X
            roblistX = zeros(mutant.mut_TS.size, DE_param.pop_size);
            % record the robustness value list of each individual U obtained by crossover
            roblistU = zeros(mutant.mut_TS.size, DE_param.pop_size);
            % X log
            X_log = cell(1, DE_param.max_gen);
            % fitnessX log
            fitnessX_log = cell(1, DE_param.max_gen);
            % roblistX log
            roblistX_log = cell(1, DE_param.max_gen);

            if strcmp(rep_mode, 'cheapincrdis')
                calib_bestInd = [];
                calib_bestFitness = 0;
                calib_bestRobList = [];
                calib_fitnessX = fitnessX;
                calib_roblistX = roblistX;
            end
            %% lookup table

            % check if a parallel pool exists
            pool = gcp('nocreate');
            % if a parallel pool exists, delete it
            if ~isempty(pool)
                delete(pool);
            end
            % create a new parallel pool with the desired configuration
            pool = parpool(core_num);

            % initialize the dictionary to store the ind-fitness maps
            indFitRobDict = [];

            mdl_bb = this.mdl_b;
            TT = this.T;
            Tss = this.Ts;
            in_namee = this.in_name;
            in_spann = this.in_span;
            phii = this.phi;
            weightt = this.weight;

            % obtain the fitness value of the initial population
            gen_i = 1;
            parfor pop_i = 1:DE_param.pop_size
                de_flag = 'or';
                alt_weight = X(pop_i, :);
                [fitnessX(1, pop_i), roblistX(:, pop_i)] = fitness(mdl_bb, TT, Tss, in_namee, in_spann, phii, weightt, alt_weight, sps_weight, sortedTS, gen_i, pop_i, rep_mode, de_flag);
            end

            % generate the ind-fitness maps based on current generation
            curIndFitRobDict = [X, fitnessX', roblistX'];
            % merge
            indFitRobDict = [indFitRobDict; curIndFitRobDict];
            % determine whether an X coresponds to a unique fitness value
            uniInd = unique(indFitRobDict(:,1:this.dim), 'rows');
            uniIndFitRobDict = unique(indFitRobDict, 'rows');
            if uniInd ~= uniIndFitRobDict(:, 1:this.dim)
                error('An X may correspond one more fitness!');
            end
            indFitRobDict = uniIndFitRobDict;

            % iterate
            for gen_i = 1:DE_param.max_gen

                % get the best fitness value
                [fitnessbestX, indexbestX] = max(fitnessX);
                % the individual whose fitness value is best in the first generation
                bestX = X(indexbestX, :);
                roblistbestX = roblistX(:,indexbestX);

                % mutation
                V = mutation(X, bestX, DE_param.F, DE_param.mut_sgy);
                % crossover
                U = crossover(X, V, lb, ub, DE_param.CR, DE_param.cr_sgy);

                parfor pop_i = 1:DE_param.pop_size
                    de_flag = 'cr';
                    alt_weight = U(pop_i, :);

                    logicIdx = ismember(indFitRobDict(:, 1:numel(alt_weight)), alt_weight, 'rows');
                    foundRows = indFitRobDict(logicIdx, :);
                    if size(foundRows, 1) > 1
                        error('duplicated items!');
                    end
                    if ~isempty(foundRows)
                        disp('Match successful!');
                        fitnessU(1, pop_i) = foundRows(1,numel(alt_weight)+1);
                        roblistU(:, pop_i) = foundRows(1,numel(alt_weight)+2:end)';
                    else
                        [fitnessU(1, pop_i), roblistU(:, pop_i)] = fitness(mdl_bb, TT, Tss, in_namee, in_spann, phii, weightt, alt_weight, sps_weight, sortedTS, gen_i, pop_i, rep_mode, de_flag);
                    end
                end

                % generate the ind-fitness maps based on current generation
                curIndFitRobDict = [U, fitnessU', roblistU'];
                indFitRobDict = [indFitRobDict; curIndFitRobDict];
                % determine whether an X coresponds to a unique fitness value
                uniInd = unique(indFitRobDict(:,1:this.dim), 'rows');
                uniIndFitRobDict = unique(indFitRobDict, 'rows');
                if uniInd ~= uniIndFitRobDict(:, 1:this.dim)
                    error('An X may correspond one more fitness!');
                end

                indFitRobDict = uniIndFitRobDict;

                % select next generation
                for pop_i = 1: DE_param.pop_size
                    if fitnessU(1, pop_i) > fitnessX(1, pop_i)
                        % find better individuals from U, then update X, fitnessX and roblistX
                        X(pop_i, :) = U(pop_i, :);
                        fitnessX(1, pop_i) = fitnessU(1, pop_i);
                        roblistX(:, pop_i) = roblistU(:, pop_i);
                        if fitnessU(1, pop_i) > fitnessbestX
                            % update bestX, fitnessbestX and roblistbestX based on fitnessU
                            bestX = U(pop_i, :);
                            fitnessbestX = fitnessU(1, pop_i);
                            roblistbestX = roblistU(:, pop_i);
                        end
                    end
                end

                % record the X, fitnessX and roblistX of current generation
                X_log{1,gen_i} = X;
                fitnessX_log{1,gen_i} = fitnessX;
                roblistX_log{1,gen_i} = roblistX;

                % display the best fitness value of the current generation
                fprintf('%d      %f\n', gen_i, fitnessbestX);
                % record the best individual and its fitness value
                bestInd(gen_i, :) = bestX;
                bestFitness(gen_i, 1) = fitnessbestX;
                bestRobList{1, gen_i} = roblistbestX;

                % although we can select and apply different fitness functions,
                % the loop exit condition is unique, i.e., #Repaired - #Broken > 0.
                if sum(bestRobList{1, gen_i} > 0) == sortedTS.size || gen_i == DE_param.max_gen
                    % when repair mode is cheapincrdis
                    if strcmp(rep_mode, 'cheapincrdis')
                        disp('cheapincrdis calibration!');
                        parfor pop_i = 1:DE_param.pop_size
                            de_flag = 'calib';
                            rep_mode = 'expdis';
                            alt_weight = X(pop_i, :);
                            [calib_fitnessX(1, pop_i), calib_roblistX(:, pop_i)] = fitness(mdl_bb, TT, Tss, in_namee, in_spann, phii, weightt, alt_weight, sps_weight, sortedTS, gen_i, pop_i, rep_mode, de_flag);
                        end
                        % get the best fitness value
                        [calib_bestFitness, calib_indexbestX] = max(calib_fitnessX);
                        % the individual whose fitness value is best in the first generation
                        calib_bestInd = X(calib_indexbestX, :);
                        calib_bestRobList = calib_roblistX(:,calib_indexbestX);
                        rep_time = toc;
                        calib_res_name = [res_name_prefix, '_geni_', num2str(gen_i), '_calib.mat'];
                        save(calib_res_name, "mutant", "DE_param", "bestInd", "bestFitness", "bestRobList", "X_log", "fitnessX_log", "roblistX_log", "calib_bestInd", "calib_bestFitness", "calib_bestRobList", "calib_fitnessX", "calib_roblistX", "rep_time");
                    end
                    return;
                end
                % if fitnessbestX reaches to 1.0 and ts_size reaches to D_size, exit DE.
                % if fitnessbestX == max_fitnessbestX && ts_size == this.D_size
                %     return;
                % elseif fitnessbestX > pre_fitnessbestX
                %     ts_size = 5 * ts_size;
                %     [test_suite, pre_fitnessbestX] = generateTestSuite(this.D, ts_size, mode);
                % end

                % save the result every DE_param.gen_span iterations
                if mod(gen_i, DE_param.gen_span) == 0
                    rep_time = toc;
                    res_name = [res_name_prefix, '_geni_', num2str(gen_i), '.mat'];
                    save(res_name, "mutant", "DE_param", "bestInd", "bestFitness", "bestRobList", "X_log", "fitnessX_log", "roblistX_log", "rep_time");
                end
            end
            % if DE fails
            % if (fitnessbestX < pre_fitnessbestX) || (fitnessbestX > pre_fitnessbestX && ts_size < this.D_size)
            %     disp('DE fails within max_gen!');
            % end

            % if the number of the safe traces is less than that before applying the alternative weights
            if sum(bestRobList{1, gen_i} > 0) <= sortedTS.sn
                disp('DE fails within max_gen!');
            end
        end
    end
end
