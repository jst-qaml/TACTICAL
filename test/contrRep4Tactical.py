import sys
# path
addpath = []
# dataset for simulation
is_nor = ''
dataset_run = ''
# dataset
dataset = []
D_size = ''
# model parameters
bench = ''
model = ''
model_m = ''
parameters = []
# io parameters
in_name = []
in_range = []
in_span = []
icc_name = []
ic_const = []
ics_name = []
oc_name = []
oc_span = []
# nn parameters
nn = ''
nn_stru = ''
# specification
phi_str = []
sig_str = []
# fl config (mi hci tops4rep)
fl_file = ''
fl_config = []
fl_mode = ''
# DE
pop_K_ratio = ''
max_gen = ''
gen_span = ''
F = ''
CR = ''
mut_sgy = ''
cr_sgy = ''
init_type = ''
init_folder = ''
# repair parameters
tsf_size = ''
pert_level = ''
search_mode = ''
rep_mode_list = []
budget = ''
# parallel computing 
core_num = ''
# script parameters
status = 0
arg = ''
linenum = 0

with open(sys.argv[1],'r') as conf:
	for line in conf.readlines():
		if not line.startswith('#'):
			argu = line.strip().split()
			if status == 0:
				status = 1
				arg = argu[0]
				linenum = int(argu[1])
			elif status == 1:
				linenum = linenum - 1
				if arg == 'addpath':
					addpath.append(argu[0])
				elif arg == 'dataset_run':
					dataset_run = argu[0]
				elif arg == 'dataset':
					dataset.append(argu[0])
				elif arg == 'D_size':
					D_size = argu[0]
				elif arg == 'is_nor':
					is_nor = argu[0]
				elif arg == 'bench':
					bench = argu[0]
				elif arg == 'model':
					model = argu[0]
				elif arg == 'model_m':
					model_m = argu[0]
				elif arg == 'parameters':
					parameters.append(argu[0])
				elif arg == 'in_name':
					in_name.append(argu[0])
				elif arg == 'in_range':
					in_range.append([float(argu[0]),float(argu[1])])
				elif arg == 'in_span':
					# print(len(argu))
					in_span = '{'
					for idx in range(len(argu)-1):
						in_span = in_span + argu[idx]+','
						# print(argu[idx])
					in_span = in_span + argu[len(argu)-1] + '}'
				elif arg == 'icc_name':
					icc_name.append(argu[0])
				elif arg == 'ic_const':
					ic_const.append(argu[0])
				elif arg == 'ics_name':
					ics_name.append(argu[0])
				elif arg == 'oc_name':
					oc_name.append(argu[0])
				elif arg == 'oc_span':
					# print(len(argu))
					oc_span = '{'
					for idx in range(len(argu)-1):
						oc_span = oc_span + argu[idx] + ','
						# print(argu[idx])
					oc_span = oc_span + argu[len(argu)-1] + '}'
				elif arg == 'nn':
					nn = argu[0]
				elif arg == 'nn_stru':
					nn_stru = argu[0]
				elif arg == 'phi_str':
					complete_phi = argu[0] + ';' + argu[1]
					for a in argu[2:]:
						complete_phi = complete_phi + ' ' + a
					phi_str.append(complete_phi)
				elif arg == 'sig_str':
					sig_str.append(argu[0])
				elif arg == 'fl_file':
					fl_file = argu[0]
				elif arg == 'fl_config':
					fl_config.append([float(argu[0]),float(argu[1]),float(argu[2])])
				elif arg == 'fl_mode':
					fl_mode = argu[0]
				elif arg == 'pop_K_ratio':
					pop_K_ratio = argu[0]
				elif arg == 'max_gen':
					max_gen = argu[0]
				elif arg == 'gen_span':
					gen_span = argu[0]
				elif arg == 'F':
					F = argu[0]
				elif arg == 'CR':
					CR = argu[0]
				elif arg == 'mut_sgy':
					mut_sgy = argu[0]
				elif arg == 'cr_sgy':
					cr_sgy = argu[0]
				elif arg == 'init_type':
					init_type = argu[0]
				elif arg == 'init_folder':
					init_folder = argu[0]
				elif arg == 'tsf_size':
					tsf_size = argu[0]
				elif arg == 'pert_level':
					pert_level = argu[0]
				elif arg == 'budget':
					budget = argu[0]
				elif arg == 'search_mode':
					search_mode = argu[0]
				elif arg == 'rep_mode':
					rep_mode_list.append(argu[0])
				elif arg == 'core_num':
					core_num = argu[0]
				else:
					continue
				if linenum == 0:
					status = 0
	# script name 
	for phi_i in range(len(phi_str)):
			ds = dataset[phi_i]
			sig_st = sig_str[phi_i]
			property = phi_str[phi_i].split(';')
			filename = model + '_spec_' + str(phi_i + 1) + '_searchmode_' + search_mode + '_budget_' + budget + '_' + fl_mode + '_Repair'
			print(filename)
			logname = '\'./output/' + filename + '.log\''
			param = '\n'.join(parameters)
			with open('scripts/'+filename,'w') as bm:
				bm.write('#!/bin/sh\n')
				bm.write('matlab -nodesktop -nosplash <<EOF\n')
				bm.write('clear;\n')
				bm.write('close all;\n')
				bm.write('clc;\n')
				bm.write('bdclose(\'all\');\n')
				# addpath
				for ap in addpath:
					bm.write('addpath(genpath(\'' + ap + '\'));\n')
				# load dataset for simulation 
				bm.write('dataset_run = \'' + dataset_run + '\';\n')
				bm.write('if strcmp(dataset_run, \'Null.mat\')\n')
				bm.write('\t D_run = \'\';\n')
				bm.write('else\n')
				bm.write('\t D_run = load(dataset_run);\n')
				bm.write('end\n')
				# load dataset
				if ds != '':
					bm.write('D = load(\'' + ds + '\');\n')
				bm.write('D_size = ' + D_size + ';\n')
				# model parameters info
				bm.write('bm = \'' + bench + '\';\n')
				bm.write('mdl = \''+ model + '\';\n')
				bm.write('mdl_m = \''+ model_m + '\';\n')
				bm.write(param + '\n')

				bm.write('is_nor = ' + is_nor + ';\n')
				bm.write('if is_nor == 1\n')
				bm.write('\t x_gain = D_run.ps_x.gain;\n')
				bm.write('\t x_gain = diag(x_gain);\n')
				bm.write('\t x_offset = D_run.ps_x.xoffset;\n')
				bm.write('\t y_gain = D_run.ps_y.gain;\n')
				bm.write('\t y_offset = D_run.ps_y.xoffset;\n')
				bm.write('end\n')

				# io parameters
				bm.write('in_name = {\'' + in_name[0] + '\'')
				for inname in in_name[1:]:
					bm.write(',')
					bm.write('\'' + inname + '\'')
				bm.write('};\n')

				bm.write('in_range = {[' + str(in_range[0][0]) + ' ' + str(in_range[0][1]) + ']')
				for ir in in_range[1:]:
					bm.write(',[' + str(ir[0]) + ' ' + str(ir[1]) + ']')
				bm.write('};\n')
				
				bm.write('in_span = ' + in_span + ';\n')

				if icc_name == []:
					bm.write('icc_name = {};\n')
				else:
					bm.write('icc_name = {\'' + icc_name[0] + '\'')
					for iccname in icc_name[1:]:
						bm.write(',')
						bm.write('\'' + iccname + '\'')
					bm.write('};\n')
					
				if ic_const == []:
					bm.write('ic_const = {};\n')
				else:
					bm.write('ic_const = {' + ic_const[0])
					for iccon in ic_const[1:]:
						bm.write(',')
						bm.write('' + iccon)
					bm.write('};\n')

				bm.write('ics_name = {\'' + ics_name[0] + '\'')
				for icsname in ics_name[1:]:
					bm.write(',')
					bm.write('\'' + icsname + '\'')
				bm.write('};\n')

				bm.write('oc_name = {\'' + oc_name[0] + '\'')
				for ocname in oc_name[1:]:
					bm.write(',')
					bm.write('\'' + ocname + '\'')
				bm.write('};\n')

				bm.write('oc_span = ' + oc_span + ';\n')

				# nn parameters
				bm.write('nn = load(\'' + nn + '\');\n')
				bm.write('net = nn.net;\n')
				bm.write('nn_stru = ' + nn_stru + ';\n')

				# specification
				bm.write('phi_str = \'' + property[1] + '\';\n')
				bm.write('spec_i = ' + str(phi_i + 1) + ';\n')
				bm.write('sig_str = ' + sig_st + ';\n')

				# fl config
				bm.write('fl_config = {[' + str(fl_config[0][0]) + ' ' + str(fl_config[0][1]) + ' ' + str(fl_config[0][2]) + ']')
				for flconf in fl_config[1:]:
					bm.write(',[' + str(flconf[0]) + ' ' + str(flconf[1]) + ' ' + str(flconf[0][2]) + ']')
				bm.write('};\n')

				bm.write('mi = fl_config{1,1}(1,1);\n')
				bm.write('hci = fl_config{1,1}(1,2);\n')
				bm.write('tops4rep = fl_config{1,1}(1,3);\n')
				bm.write('flRes = load(\'' + fl_file + '\');\n')
				bm.write('fl_mode = \'' + fl_mode + '\';\n')
				# DE
				bm.write('DE_param.pop_K_ratio = ' + pop_K_ratio + ';\n')
				bm.write('DE_param.max_gen = ' + max_gen + ';\n')
				bm.write('DE_param.K = tops4rep * nn_stru(1,1);\n')
				bm.write('DE_param.pop_size = DE_param.pop_K_ratio * DE_param.K;\n')
				bm.write('DE_param.gen_span = ' + gen_span + ';\n')
				bm.write('DE_param.F = ' + F + ';\n')
				bm.write('DE_param.CR = ' + CR + ';\n')
				bm.write('DE_param.mut_sgy = ' + mut_sgy + ';\n')
				bm.write('DE_param.cr_sgy = ' + cr_sgy + ';\n')
				bm.write('DE_param.init_type = \'' + init_type + '\';\n')
				bm.write('DE_param.init_folder = \'' + init_folder + '\';\n')
				
				# repair parameters
				bm.write('tsf_size = ' + tsf_size + ';\n')
				bm.write('pert_level = ' + pert_level + ';\n')
				bm.write('budget = ' + budget + ';\n')
				bm.write('search_mode = \'' + search_mode + '\';\n')
				bm.write('rep_mode_list = {\'' + rep_mode_list[0] + '\'')
				for repmode in rep_mode_list[1:]:
					bm.write(',')
					bm.write('\'' + repmode + '\'')
				bm.write('};\n')
				
				# parallel computing
				bm.write('core_num = ' + core_num + ';\n')
				bm.write('rng(round(rem(now, 1)*1000000));\n')
				bm.write('bm_cur = TACTICAL(bm, mdl, mdl_m, D_run, is_nor, D, D_size, net, nn_stru, T, Ts, in_name, in_range, in_span, icc_name, ic_const, ics_name, oc_name, oc_span, phi_str, sig_str, core_num);\n')
				
				bm.write('spsNeuronList = flRes.hyMCon(mi,hci).spsNeuronList(1:tops4rep,:);\n')
				bm.write('spsWeightList = [repelem(spsNeuronList(:,1:2), nn_stru(1,1), 1), repmat((1:nn_stru(1,1))\', tops4rep, 1)];\n')
				bm.write('spsWeightValList = zeros(size(spsWeightList,1),1);\n')
				bm.write('for i = 1:size(spsWeightList,1)\n')
				bm.write('\t spsWeightValList(i,1) = bm_cur.weight{1,spsWeightList(i,1)}(spsWeightList(i,2),spsWeightList(i,3));\n')
				bm.write('end\n')
				bm.write('spsWeightList = [repelem(spsNeuronList(:,3), nn_stru(1,1), 1), spsWeightValList, repelem(spsNeuronList(:,1:2), nn_stru(1,1), 1), repmat((1:nn_stru(1,1))\', tops4rep, 1)];\n')
				
				# convert current CPS to a mutant
				bm.write('mutant = struct();\n')
				bm.write('mut_TS = struct();\n')
				bm.write('mut_TS.size = tsf_size;\n')
				bm.write('mut_TS.tr_in_cell = D.tr_in_cell(1,75:tsf_size+74);\n')
				bm.write('mut_TS.tr_rob_set = D.tr_rob_set(75:tsf_size+74,:);\n')
				bm.write('mut_TS.sn = sum(D.tr_rob_set(75:tsf_size+74,:) > 0);\n')
				bm.write('mutant.mut_TS = mut_TS;\n')
				bm.write('mutant.phi = bm_cur.phi;\n')
				bm.write('mutant.bug_lb = bm_cur.min_weight;\n')
				bm.write('mutant.bug_ub = bm_cur.max_weight;\n')
				bm.write('mutant.sps_weight = spsWeightList;\n')
				bm.write('mutant.bug_mdl = bm_cur.mdl_m;\n')
				
				bm.write('for rmi = 1:numel(rep_mode_list)\n')
				bm.write('\t cur_rep_mode = rep_mode_list{1, rmi};\n')
				bm.write('\t for bi = 1:budget\n')
				bm.write('\t\t repLog_ID = struct;\n')
				bm.write('\t\t repLog_ID.K = tops4rep * nn_stru(1,1);\n')
				bm.write('\t\t repLog_ID.bi = bi;\n')
				
				bm.write('\t\t res_name = [mdl_m, \'_spec_\', num2str(spec_i), \'_K_\', num2str(DE_param.K), \'_popsize_\', num2str(DE_param.pop_size), \'_tsfsize_\', num2str(tsf_size), \'_searchmode_\', search_mode, \'_repmode_\', cur_rep_mode, \'_budget_\', num2str(bi), \'_\', fl_mode, \'.mat\'];\n')
				bm.write('\t\t res_name_prefix = [mdl_m, \'_spec_\', num2str(spec_i), \'_K_\', num2str(DE_param.K), \'_popsize_\', num2str(DE_param.pop_size), \'_tsfsize_\', num2str(tsf_size), \'_searchmode_\', search_mode, \'_repmode_\', cur_rep_mode, \'_budget_\', num2str(bi), \'_\', fl_mode];\n')
				bm.write('\t\t tic;\n')
				
				bm.write('\t\t [bestInd, bestFitness, bestRobList, X_log, fitnessX_log, roblistX_log, mutant] = bm_cur.altWeightSearch(mutant, fl_mode, cur_rep_mode, pert_level, DE_param, repLog_ID, core_num, res_name_prefix);\n')

				bm.write('\t\t rep_time = toc;\n')
				bm.write('\t\t save(res_name, "mutant", "DE_param", "bestInd", "bestFitness", "bestRobList", "X_log", "fitnessX_log", "roblistX_log", "rep_time");\n')
				bm.write('\t\t delete(fullfile([\'breach/Ext/ModelsData/\', \'*_breach.slx\']));\n')
				bm.write('\t\t delete(fullfile([\'*.slx\']));\n')
				bm.write('\t\t delete(fullfile([\'*.slx.autosave\']));\n')
				bm.write('\t\t delete(fullfile([\'*.slx.r202*\']));\n')
				bm.write('\t\t delete(fullfile([\'*_breach.slxc\']));\n')
				bm.write('\t end\n')
				bm.write('end\n')
		
				bm.write('logname = ' + logname + ';\n')
				bm.write('sendEmail(logname);\n')
				bm.write('quit force\n')
				bm.write('EOF\n')