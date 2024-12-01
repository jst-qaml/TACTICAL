function [hyCon] = isActivated(colNeuronOut, hyCon, nnStru, state, minOut, maxOut)
% isActivated function returns hyCon with updated actSpectrum.
%
% Inputs:
%   colNeuronOut: the values of the neurons in the hidden layer
%   hyCon: hyper config
%   nnStru: nn structure
%   state: 1: pass; 0: fail
%   minOut: the minimum value of all neurons of current mutant across all available test cases
%   maxOut: the maximum value of all neurons of current mutant across all available test cases
% Outputs:
%   hyCon: hyper config with updated actSpectrum

for hci = 1:numel(hyCon)
    curCri = hyCon(1,hci).criterion;
    hyperval1 = hyCon(1,hci).hyperval1;
    hyperval2 = hyCon(1,hci).hyperval2;
    switch curCri
        case 'INA'
            % INA (h) sml  {0,5%,10%}
            th = minOut + hyperval1 * (maxOut - minOut);
            actState = INA(colNeuronOut, th);
        case 'ITK'
            % ITK (k) sml  {1,2,3}/{2,3,4}/{3,4,5}
            actState = ITK(colNeuronOut, hyperval1);
        case 'PNA'
            % PNA (h, deltaT) ss sm asl ms mm ml ls lm ll  {0,5%,10%} * {0.5,1.0,1.5}
            th = minOut + hyperval1 * (maxOut - minOut);
            actState = PNA(colNeuronOut, th, hyperval2 * 10); % time to length
        case 'PTK'
            % PTK (deltaT, k) ss sm sl ms mm ml ls lm ll  {0.5,1.0,1.5} * {1,2,3}/{2,3,4}/{3,4,5}
            actState = PTK(colNeuronOut, hyperval1 * 10, hyperval2); 
        case 'PD'
            % PD  (h, deltaT) ss sm sl ms mm ml ls lm ll  {0,5%,10%} * {0.5,1.0,1.5}
            th = minOut + hyperval1 * (maxOut - minOut);
            actState = PD(colNeuronOut, th, hyperval2 * 10);
        case 'ND'
            % ND  (h, deltaT) ss sm sl ms mm ml ls lm ll  {0,5%,10%} * {0.5,1.0,1.5}
            th = minOut + hyperval1 * (maxOut - minOut);
            actState = ND(colNeuronOut, th, hyperval2 * 10);
        case 'MI'
            % MI  (deltaT) sml  {0.5,1.0,1.5}
            actState = MI(colNeuronOut, hyperval1 * 10);
        case 'MD'
            % MD  (deltaT) sml  {0.5,1.0,1.5}
            actState = MD(colNeuronOut, hyperval1 * 10);
        otherwise
            error("Please check the selected criterion!")
    end

    hLNum = numel(nnStru);
    newActState = cell(1, hLNum);
    for li = 2:hLNum
        layerActSpectrum = actState{1,li};
        layerActSpectrum(layerActSpectrum ~= 0) = 1;
        newActState{1,li} = layerActSpectrum;
    end

    % update activation spectra based on activation state of each neuron
    for li = 2:hLNum
        if state == 1
            % update the activation spectrum of the system executions For each neuron, [ap, np, af, nf]
            hyCon(1,hci).actSpectrum{1,li}(:,1) = hyCon(1,hci).actSpectrum{1,li}(:,1) + newActState{1,li};
            hyCon(1,hci).actSpectrum{1,li}(:,2) = hyCon(1,hci).actSpectrum{1,li}(:,2) + 1 - newActState{1,li};
        else
            hyCon(1,hci).actSpectrum{1,li}(:,3) = hyCon(1,hci).actSpectrum{1,li}(:,3) + newActState{1,li};
            hyCon(1,hci).actSpectrum{1,li}(:,4) = hyCon(1,hci).actSpectrum{1,li}(:,4) + 1 - newActState{1,li};
        end
    end
end
end