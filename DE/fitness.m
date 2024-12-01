function [fitnessVal, robList] = fitness(mdl_b, T, Ts, in_name, in_span, phi, weight, alt_weight, sps_weight, sortedTS, gen_i, pop_i, rep_mode, de_flag)
% The fitness function defines different approaches to calculate the final fitness value of the given alt_weight.
%
% Inputs:
%   mdl_b: the copy of the buggy model for alternative weights search
%   T: simulation time
%   Ts: sampling interval
%   in_name: the input signals of the system
%   in_span: the time span of the input signals
%   phi: specification in the form of STL
%   weight: the original weights of the nn controller
%   alt_weight: a set of alternative weights to replace
%   sps_weight: includes layer_idx, right_endpoint_idx, left_endpoint_idx and weight value and suspiciousness score of each suspicious weight
%   sortedTS: the sorted test suite
%   gen_i: generation idx
%   pop_i: individual idx in the whole population, for namingidx is used in the combination (this.mdl_m + idx) to name a model for fitness value calculation
%   rep_mode: repair mode
%   de_flag:
% Outputs:
%   fitnessVal: fitness value
%   robList: the robustness values of current test suite

% evalin can extract information from base workspace. But it may cause performance degradation.
% in_span = evalin('base', 'in_span');

% The following six terminologies of the possible repair results.
% Worsen: decrease robustness
% Worsen_negative_negtive (-1 -> -2)
% Worsen_positive_positive (2 -> 1)
% Worsen_positive_negative (2 -> -1)(Broken)
% Improve: increase robustness
% Improve_positive_positive (1 -> 2)
% Improve_negative_negtive (-2 -> -1)
% Improve_negative_positive (-1 -> 1)(Repaired)

% apply the alternative weights to mdl_new
mdl_new = [mdl_b, '_g_', num2str(gen_i), '_p_', num2str(pop_i), '_', de_flag];
weightReplace(mdl_b, weight, alt_weight, sps_weight, mdl_new);

% Here, we modify the mdl_m directly (replace 'T' and 'Ts' with constants)
% set_param(mdl_new, 'StopTime', num2str(T));
% set_param(mdl_new, 'FixedStep', num2str(Ts));

% reorganize sortedTS to adapt different repair modes
if strcmp(rep_mode, 'cheapdis') || strcmp(rep_mode, 'cheapcon1') || strcmp(rep_mode, 'cheapcon2')
    reoTS = reorganizeTS(TS, subts_sz, rep_mode);
else
    % rep_mode == 'expdis', 'expcon1', 'expcon2'
    reoTS = sortedTS;
end

%% decide the way to go through the given test suite
if strcmp(rep_mode, 'expdis') || strcmp(rep_mode, 'expcon1')
    robList = zeros(reoTS.size, 1);
    % perform simulations based on reoTS and obtain necessary info for fitness calculation
    for in_i = 1:reoTS.size
        in_sig = reoTS.tr_in_cell{1,in_i};
        [rob, ~] = calRob(mdl_new, T, Ts, in_name, in_span, phi, in_sig);
        robList(in_i, 1) = rob;
    end
elseif strcmp(rep_mode, 'cheapincrdis')
    % Approximate fitness, version IncrementalTests
    robList = zeros(reoTS.size, 1);
    pos_end_idx = reoTS.sn;
    neg_start_idx = reoTS.sn + 1;

    robList(1:pos_end_idx,1) = Inf;
    robList(neg_start_idx:end, 1) = -Inf;
    % perform simulations based on reoTS and obtain necessary info for fitness calculation
    % for negative
    for in_i = neg_start_idx:reoTS.size
        in_sig = reoTS.tr_in_cell{1,in_i};
        [rob, ~] = calRob(mdl_new, T, Ts, in_name, in_span, phi, in_sig);
        robList(in_i,1) = rob;
        if rob < 0
            break;
        end
    end
    % for positive
    for in_i = pos_end_idx:-1:1
        in_sig = reoTS.tr_in_cell{1,in_i};
        [rob, ~] = calRob(mdl_new, T, Ts, in_name, in_span, phi, in_sig);
        robList(in_i,1) = rob;
        if rob > 0
            break;
        end
    end
end

% obtain the indexes of state transition
% improve_neg_neg_idx = ;
% improve_neg_pos_idx = ;
% improve_pos_pos_idx = ;
%
% worsen_pos_pos_idx = ;
% worsen_pos_neg_idx = ;
% worsen_pos_pos_idx = ;
%% expensive (the entire test suite) and discrete (this fitness function is directly from Arachne)
% #Repaired - #Broken
if strcmp(rep_mode, 'expdis')
    fitnessVal = sum(robList > 0) - sortedTS.sn;
end
%% expensive (the entire test suite) and continuous_1 (in this function, we accumulate the delta rob rather than simply counting the safe traces)
% Sum of ΔRob for Improve_negative + sum of ΔRob for Worsen + 10K * (#Repaired - #Broken)
if strcmp(rep_mode, 'expcon1')
    fitnessVal = sum(robList) - sum(sortedTS.tr_rob_set) + 10000 * (sum(robList > 0) - sortedTS.sn);
end
%% expensive (the entire test suite) and continuous_2
if strcmp(rep_mode, 'expcon2')
    % Sum of ΔRob for Improve_negative + sum of ΔRob for (Worsen_negative + Worsen_positive_negative) + 10K * (#Repaired - #Broken)
    % Compared to expcon1, expcon2 doesn't consider the Worsen_positive_positive.
    fitnessVal = 0;
end
%% approximate subset_1 and discrete (in this function, we use fewer samples to improve the speed of alternative weights search process)
if strcmp(rep_mode, 'cheapsub1dis')
    fitnessVal = 0;
end
% approximate subset_2 and continuous_1
if strcmp(rep_mode, 'cheapsub1con1')
    fitnessVal = 0;
end
if strcmp(rep_mode, 'cheapsub1con2')
    fitnessVal = 0;
end

if strcmp(rep_mode, 'cheapsub2dis')
    fitnessVal = 0;
end

if strcmp(rep_mode, 'cheapsub2con1')
    fitnessVal = 0;
end

if strcmp(rep_mode, 'cheapsub2con2')
    fitnessVal = 0;
end
if strcmp(rep_mode, 'cheapincrdis')
    fitnessVal = sum(robList > 0) - sortedTS.sn;
end
end