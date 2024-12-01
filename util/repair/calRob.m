function [rob, Br] = calRob(mdl, T, Ts, in_name, in_span, phi, in_sig)
% Given a model and an input signal, calRob function returns the robustness and Br, according to the specification.
%
% Inputs:
%   mdl: a simulink model
%   in_sig: an external input signal of the system
% Outputs:
%   rob: robustness value
%   Br: an instance of the class BreachSimulinkSystem

% signal diagnosis
Br = BreachSimulinkSystem(mdl);
Br.Sys.tspan = 0:Ts:T;

if length(unique(in_span)) == 1
    con_type = 'UniStep';
    input_gen.cp = T/in_span(1,1);
else
    con_type = 'VarStep';
    input_gen.cp = T./in_span;
end

input_gen.type = con_type;
Br.SetInputGen(input_gen);

if strcmp(con_type, 'UniStep')
    for i = 1:numel(in_name)
        for cpi = 0:input_gen.cp - 1
            eval(['Br.SetParam({''',in_name{1,i},'_u',num2str(cpi),'''}, in_sig{1,i}(1,cpi+1));']);
        end
    end
elseif strcmp(con_type, 'VarStep')
    for i = 1:numel(in_name)
        for cpi = 0: input_gen.cp(1,i) - 1
            eval(['Br.SetParam({''',in_name{1,i},'_u',num2str(cpi),'''}, in_sig{1,i}(1,cpi+1));']);
            if cpi ~= input_gen.cp(1,i) - 1
                eval(['Br.SetParam({''',in_name{1,i},'_dt',num2str(cpi),'''},in_span(1,i));']);
            end
        end
    end
end

Br.Sim(0:Ts:T);
rob = Br.CheckSpec(phi);
end