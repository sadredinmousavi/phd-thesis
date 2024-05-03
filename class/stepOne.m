classdef stepOne < definition
    %STEPONE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        args
        configs
        handles
        savingData
        plotOptions
        symbolicOutput = {}
        designPsaiOutput = {}
        solverInputs = {}
        solverOutputs = {}
        
    end
    
    methods
        function obj = stepOne(handles, args, configs, plotOptions)
            %STEPONE Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                tmp = matlab.desktop.editor.getActive;
                cd(fileparts(tmp.Filename));
                cd('..');
                addpath(genpath('utils'));
                addpath(genpath('functions'));
                addpath(genpath('inputFiles'));
                addpath(genpath('data'));
                rmpath(genpath('backups'));
                rmpath(genpath('backup_v2'));
                obj = loadFromFile(obj, handles);
            else
                obj.handles = handles;
                obj.args = args;
                obj.configs = configs;
                obj.plotOptions = plotOptions;

                % Define filing configs
                tmp = matlab.desktop.editor.getActive;
                cd(fileparts(tmp.Filename));
                cd('..');
                addpath(genpath('utils'));
                addpath(genpath('functions'));
                addpath(genpath('inputFiles'));
                addpath(genpath('data'));
                rmpath(genpath('backups'));
                rmpath(genpath('backup_v2'));
                counterName = 1;
                cd('data')
                fileName1 = '';
                while strcmp(fileName1, '')
                    if isfile([num2str(counterName) '_object.mat'])
                        counterName = counterName + 1;
                    else
                        fileName1 = [num2str(counterName) '_input.mat'];
                        obj.savingData.fileName1 = [num2str(counterName) '_input.mat'];
                        obj.savingData.fileName2 = [num2str(counterName) '_movie.gif'];
                        obj.savingData.fileName3 = [num2str(counterName) '_datas.mat'];
                        obj.savingData.fileName4 = [num2str(counterName) '_object.mat'];
                    end
                end
                cd('..')

                prepareBeforeRun(obj)
                saveToFile(obj)
                disp(obj.savingData.fileName4)
%                 runSimulations(obj)
%                 plotOutput(obj)
            end
        end
        
        function saveToFile(obj, newName)
            if nargin > 1
                obj.savingData.fileName4 = newName;
            end
            tmp = matlab.desktop.editor.getActive;
            cd(fileparts(tmp.Filename));
            cd('..\data');
            save(obj.savingData.fileName4, "obj")
            cd('..')
        end
        
        function obj = loadFromFile(obj, fileName)
            cd('data')
            load(fileName);
            cd('..')
        end
        
        function [x_eqPoints, y_eqPoints] = calcEqPointByDynamics(obj, Psai)
            [x_eqPoints, y_eqPoints] = findEqPoints_Dynamics(obj.configs.vars.x_space ,obj.configs.vars.y_space, Psai, 50, 1);
        end
        
        function obj = prepareBeforeRun(obj)
            obj = symbolicCalculations(obj);
            obj = designPsai(obj);
            obj = seedAndParamsForMRs(obj);
            obj = seedAndParamsForFPs(obj);
        end
        
        function obj = writeFilesBeforeRun(obj)
            makeFunctionsFromSymbolic(obj.symbolicOutput.f, obj.symbolicOutput.b, obj.symbolicOutput.b_m, obj.symbolicOutput.str_psai, obj.symbolicOutput.str_psai_sym)
            makeFunctionsFromDesignPsai(obj.designPsaiOutput.time, obj.designPsaiOutput.eqP, obj.designPsaiOutput.Psai_t, obj.designPsaiOutput.x_eqPoints, obj.designPsaiOutput.y_eqPoints)
        end
        
        function obj = setEqPoints(obj, eqPoints)
            obj.configs.simulations.endTime       = max(max(max(eqPoints{:}))) + 20;
            obj.configs.simulations.tspan         = obj.configs.simulations.startTime:obj.configs.simulations.stepForOutput:obj.configs.simulations.endTime;
            obj.configs.eqPoints                  = eqPoints;
        end
        
        function obj = setWalls(obj, walls)
            obj.configs.walls = walls;
        end
        
        function obj = setMRsLoc(obj, seedData)
            obj.configs.mrs.seedData = seedData;
            seedAndParamsForMRs(obj);
        end
        
        function obj = seedAndParamsForMRs(obj)
            r_mr = obj.args.mr.D/2;
            m_mr = obj.args.mr.mass;
            k_mr = obj.args.mr.k;
            groups = obj.configs.mrs.seedData;
            x_mr_0 = [];
            y_mr_0 = [];
            for i=1:length(groups)
                x_serie_ = linspace(groups{i}.center(1)-groups{i}.offset, groups{i}.center(1)+groups{i}.offset, groups{i}.num);
                y_serie_ = linspace(groups{i}.center(2)-groups{i}.offset, groups{i}.center(2)+groups{i}.offset, groups{i}.num);
                [x_mr__,y_mr__] = meshgrid(x_serie_,y_serie_);
                x_mr_0 = [x_mr_0 reshape(x_mr__, 1, [])];
                y_mr_0 = [y_mr_0 reshape(y_mr__, 1, [])];
            end
            r_mr_0 = r_mr * ones(1, length(x_mr_0));
            m_mr_0 = m_mr * ones(1, length(x_mr_0));
            k_mr_0 = k_mr * ones(2, length(x_mr_0));
            %
            obj.configs.mrs.x_mr_0 = x_mr_0;
            obj.configs.mrs.y_mr_0 = y_mr_0;
            obj.configs.mrs.r_mr_0 = r_mr_0;
            obj.configs.mrs.m_mr_0 = m_mr_0;
            obj.configs.mrs.k_mr_0 = k_mr_0;
        end
        
        function obj = setFPsLoc(obj, fps)
            obj.configs.fps.fps = fps;
            seedAndParamsForFPs(obj);
        end
        
        function obj = seedAndParamsForFPs(obj)
            rho_fp = obj.args.fp.rho;
            fpsSeed = obj.configs.fps.fps;
            x_fp_0 = [];
            y_fp_0 = [];
            t_fp_0 = [];
            m_fp_0 = [];
            i_fp_0 = [];
            k_fp_0 = [];
            for i=1:length(fpsSeed)
                x_fp_0 = [x_fp_0 fpsSeed{i}.center(1)];
                y_fp_0 = [y_fp_0 fpsSeed{i}.center(2)];
                t_fp_0 = [t_fp_0 fpsSeed{i}.theta0];
                if fpsSeed{i}.type == 1
                    V = pi * fpsSeed{i}.height * fpsSeed{i}.radius^2;
                    m = rho_fp * V;
                    I = 1/2 * m * fpsSeed{i}.radius^2;
                    k = [16/3; 16/3] * fpsSeed{i}.radius;
                    m_fp_0 = [m_fp_0 m];
                    i_fp_0 = [i_fp_0 I];
                    k_fp_0 = [k_fp_0 k];
                elseif fpsSeed{i}.type == 2
                    a=1;
                elseif fpsSeed{i}.type == 3
                    a=1;
                end
            end
            %
            obj.configs.fps.x_fp_0 = x_fp_0;
            obj.configs.fps.y_fp_0 = y_fp_0;
            obj.configs.fps.t_fp_0 = t_fp_0;
            obj.configs.fps.m_fp_0 = m_fp_0;
            obj.configs.fps.i_fp_0 = i_fp_0;
            obj.configs.fps.k_fp_0 = k_fp_0;
        end
        
        function obj = setEpmsLoc(obj, magNum, a, z)
            for i=1:magNum
                phi(i) = (0 + (i-1)*(360/magNum) )*(pi/180);
                psai(i) = 0;
                MagPos(i,:) = [a*cos(phi(i)) a*sin(phi(i)) z 1 0 0];
            %     MagPos(i,4:6) = round(MagPos(i,4:6) ./ norm(MagPos(i,4:6)), 5);
            end

            obj.configs.epms.magNum = magNum;
            obj.configs.epms.a = a;
            obj.configs.epms.z = z;
            obj.configs.epms.MagPos = MagPos;
            obj.configs.epms.phi = phi;
            obj.configs.epms.psai = psai;
            %
            obj = prepareBeforeRun(obj);
        end
        
        function obj = setConfigVars(obj, domain, step, plotDomain)
            if nargin < 4
                plotDomain = max(a) + 0.02;
            end
            spaceRegion = -domain:step:domain;
            x_space = spaceRegion; z_space = x_space; y_space = x_space;
            Npoints_space = length(x_space);


            obj.configs.vars.domain = domain;
            obj.configs.vars.step = step;
            obj.configs.vars.plotDomain = plotDomain;
            obj.configs.vars.Npoints_space = Npoints_space;
            obj.configs.vars.x_space = x_space;
            obj.configs.vars.y_space = y_space;
            obj.configs.vars.z_space = z_space;
        end
        
        function obj = symbolicCalculations(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            dlgTitle    = 'PreCalculations Symbolic';
            dlgQuestion = 'Do you want symbolic calculations to be done ?';
            choice = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes');
            if strcmp(choice, 'Yes') == 1
                obj.handles.DoPreCalc = 1;
                obj.symbolicOutput = obj.handles.symbolicFunctionHandle(obj.configs.epms, obj.args);
                %
                makeFunctionsFromSymbolic(obj.symbolicOutput.f, obj.symbolicOutput.b, obj.symbolicOutput.b_m, obj.symbolicOutput.str_psai, obj.symbolicOutput.str_psai_sym)
            else
                obj.handles.DoPreCalc = 0;
            end
            
        end
        
        function obj = designPsai(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            dlgTitle    = 'PreCalculations Psai';
            dlgQuestion = 'Do you want Design Psai calculations to be done ?';
            choice = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes');
            if strcmp(choice, 'Yes') == 1
                obj.handles.DoDesignPsai = 1;
                obj.designPsaiOutput = designPsaiController(obj.configs, obj.handles);
                %
                makeFunctionsFromDesignPsai(obj.designPsaiOutput.time, obj.designPsaiOutput.eqP, obj.designPsaiOutput.Psai_t, obj.designPsaiOutput.x_eqPoints, obj.designPsaiOutput.y_eqPoints)
            else
                obj.handles.DoDesignPsai = 0;
            end
            
        end
        
        function obj = runSimulations(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            dlgTitle    = 'Simulations';
            dlgQuestion = 'Do you want to run the simulations ?';
            choice = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes');
            if strcmp(choice, 'Yes') == 1
                writeFilesBeforeRun(obj)
                [obj.solverInputs, obj.solverOutputs] = solveTheSystem(obj.configs, obj.args, obj.savingData);
            end
        end
        
        function obj = plotOutput(obj, plotOptionsFromInput)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if nargin > 1
                plotOptions1 = plotOptionsFromInput;
            else
                plotOptions1 = obj.plotOptions.dynamic;
            end
            plotOutputToGif(obj.solverInputs, obj.solverOutputs, obj.configs, obj.savingData, plotOptions1)
        end
        
        
    end
end

