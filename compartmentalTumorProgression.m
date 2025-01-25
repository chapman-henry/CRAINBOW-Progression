classdef compartmentalTumorProgression < handle
    %COMPARTMENTALTUMORPROGRESSION : methods for estimating parameters for
    %compartmental ODE models describing tumor progression in CRAINBOW
    %HER2+ mouse survey data


    properties
        mouse_survey
        size_threshold
        n_mice
        mouse_list
        ages_real
        ages_jitter
        d16_counts
        p95_counts
        d16_Field_params
        p95_Field_params
        param_estimate
        sum_of_squares
        edges
        t_start
        delta_t = 0.01
        t_final
        t_max_F_growth
    end

    methods
        function obj = compartmentalTumorProgression(cleanedTumorSurvey,P_d16,P_p95,small_large_cutoff)
            % TODO: 2) max F growth rate time

            %CONSTRUCTOR
            % INPUT:
            % cleanedTumorSurvey -- table of cleaned HERBOW tumor survey
            % data collected by Josh Ginzel
            % colums:
            %   ID (categorical) --- unique mouse id
            %   Age (double) --- mouse age (weeks)
            %   Gland (str?) ---
            %   Genotype (categorical) --- one of 'wt', 'd16', 'p95'
            %   Area (double) --- 2 dim imaged area of tumor (micron^2)
            %   Areamm2 (double) --- 2 dim imaged area of tumor (mm^2)
            %   Mm2 (double?) --- large flag?
            % P_d16 (struct) with 2 fields: F0, K --- gives parameters for F compartment logistic growth
            % P_p95 (struct) with 2 fields: "" ditto
            %
            obj.mouse_survey = cleanedTumorSurvey;
            obj.mouse_list = unique(obj.mouse_survey.ID);
            obj.n_mice = length(obj.mouse_list);
            % default parameter values for logistic growth of transformed
            % field cells
            % P_d16.F0 = 19040;
            % P_d16.K = 76950;
            % P_d16.r = 1.1123;
            %
            % P_p95.F0 = 4760;
            % P_p95.K = 16200;
            % P_p95.r = 0.8769;
            %obj.t_max_F_growth = max_puberty_time;
            %P_d16.r = obj.calcGrowthRate(P_d16.F0,P_d16.K,max_puberty_time);
            %P_p95.r = obj.calcGrowthRate(P_p95.F0,P_p95.K,max_puberty_time);

            obj.d16_Field_params = P_d16;
            obj.p95_Field_params = P_p95;

            obj.t_start = 3; %first observations of Field at week 3
            %obj.delta_t = 0.01;
            obj.t_final = 26; % max observed age at 25

            obj.size_threshold = small_large_cutoff;
            obj.edges = [min(obj.mouse_survey.Areamm2), obj.size_threshold, max(obj.mouse_survey.Areamm2)];

        end %constructor




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DATA PRE-PROCESSING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function mouseSummary(obj)
            % compiles number of Small and Large tumors, as defined by the scalar
            % 'threshold' for both d16 and p95 genotypes
            % relies on matlab function datasample
            % obj.size_threshold = threshold;
            % edges = [min(obj.mouse_survey.Areamm2), threshold, max(obj.mouse_survey.Areamm2)];
            % mouseList = unique(obj.mouse_survey.ID);
            varTypes = {'string', 'double', ...
                'double', 'double',	 ...
                'double', 'double'	};
            varNames = {'ID', 'AGE', ...
                'd16_1', 'd16_2'	 ...
                'p95_1', 'p95_2'	};
            sumData = table('Size', [obj.n_mice, length(varNames)], ...
                'VariableTypes', varTypes, 'VariableNames', varNames);
            for ii = 1:obj.n_mice
                tempMat = obj.mouse_survey(obj.mouse_survey.ID == obj.mouse_list(ii), :);
                ID = obj.mouse_list(ii);
                AGE = tempMat.Age(1);
                R = {ID, AGE};
                hasd16 = strcmpi(string(tempMat.Genotype),'d16');
                hasp95 = strcmpi(string(tempMat.Genotype),'p95');
                b = horzcat(hasd16, hasp95);
                for jj = 1:2
                    a = obj.genoCalc(tempMat(b(:,jj),:),obj.edges);
                    R = horzcat(R,a);
                end
                sumData(ii,:) = R;
            end


            %jitter the ages for the optimization function
            jitter_age = randn(obj.n_mice,1)*0.05;

            real_ages = sumData.AGE;
            jitter_ages = real_ages + jitter_age;

            [~,idx] = sort(jitter_ages);

            obj.ages_real = real_ages(idx);
            obj.ages_jitter = jitter_ages(idx);

            % writing data
            obj.d16_counts = horzcat(sumData.d16_1(idx),sumData.d16_2(idx));
            obj.p95_counts = horzcat(sumData.p95_1(idx),sumData.p95_2(idx));


            % we DON'T WANT mouse jittered ages to be the same
            % it will break call to lsqcurvefit()
            for ii = 2:obj.n_mice
                if obj.ages_jitter(ii-1) == obj.ages_jitter(ii)
                    % make ages different
                    obj.ages_jitter(ii,1) = obj.ages_jitter(ii,1) + 0.01;
                end
            end

        end %mouseSummary method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function bootstrapMouseSummary(obj)
            % creates 2-level bootstrap of mice and tumors, then
            % compiles number of Small and Large tumors, as defined by the scalar
            % 'threshold' for both d16 and p95 genotypes

            % move up one level.  Write to obj properties
            %obj.size_threshold = threshold;
            %obj.edges = [min(obj.mouse_survey.Areamm2), threshold, max(obj.mouse_survey.Areamm2)];

            %mouseList = unique(obj.mouse_survey.ID);
            resam_mouseList = datasample(obj.mouse_list, height(obj.mouse_list));
            varTypes = {'string', 'double', ...
                'double', 'double',	 ...
                'double', 'double'	};
            varNames = {'ID', 'AGE', ...
                'd16_1', 'd16_2'	 ...
                'p95_1', 'p95_2'	};
            sumData = table('Size', [length(resam_mouseList), length(varNames)], ...
                'VariableTypes', varTypes, 'VariableNames', varNames);
            %loop over mice
            for ii = 1:length(resam_mouseList)

                myMat = obj.mouse_survey(obj.mouse_survey.ID == resam_mouseList(ii), :);
                ID = resam_mouseList(ii);
                %bootstrap 2/2: sample tumors from mouse ii
                tempMat = datasample(myMat,height(myMat),1);
                AGE = tempMat.Age(1);
                R = {ID, AGE};
                hasd16 = strcmpi(string(tempMat.Genotype),'d16');
                hasp95 = strcmpi(string(tempMat.Genotype),'p95');
                b = horzcat(hasd16, hasp95);
                for jj = 1:2
                    a = obj.genoCalc(tempMat(b(:,jj),:),obj.edges);
                    R = horzcat(R,a);
                end
                sumData(ii,:) = R;
            end %loop over mice

            obj.n_mice = height(sumData);

            %jitter the ages for the optimization function
            jitter_age = randn(obj.n_mice,1)*0.05;

            real_ages = sumData.AGE;
            jitter_ages = real_ages + jitter_age;

            [~,idx] = sort(jitter_ages);

            obj.ages_real = real_ages(idx);
            obj.ages_jitter = jitter_ages(idx);

            % writing data
            obj.d16_counts = horzcat(sumData.d16_1(idx),sumData.d16_2(idx));
            obj.p95_counts = horzcat(sumData.p95_1(idx),sumData.p95_2(idx));


            % we DON'T WANT mouse jittered ages to be the same
            % it will break call to lsqcurvefit()
            for ii = 2:obj.n_mice
                if obj.ages_jitter(ii-1) == obj.ages_jitter(ii)
                    % make ages different
                    obj.ages_jitter(ii,1) = obj.ages_jitter(ii,1) + 0.01;
                end
            end

        end %bootstrapMouseSummary
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function a = genoCalc(obj,genoData,edges)
            % helper function for mouseSummary method
            [a,~] = histcounts(genoData.Areamm2,edges);
            a = num2cell(a);
        end %function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PARAMETER ESTIMATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function estimate_bootstrap_parameters(obj,genotype,model,n_samples)
            %INPUT:
            % genotype (str) = one of 'd16' or 'p95'
            % model (str) = one of '2C', '3C_partialObs', '3C_completeObs'
            % n_samples (int) = number of bootstrap iterations

            % allocating space for parameter estimates
            % K, F0, r, <p0>, p1, p2
            if strcmp(model,'2C')
                obj.param_estimate = zeros(n_samples,5);
            else
                obj.param_estimate = zeros(n_samples,6);
            end

            obj.sum_of_squares = zeros(n_samples,1);

            for ii = 1:n_samples
                obj.bootstrapMouseSummary;
                [obj.param_estimate(ii,:),obj.sum_of_squares(ii)] = obj.estimate_parameters(genotype,model);
            end

        end %estimate_boostrap_params
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [param_estimate,sum_of_squares] = estimate_parameters(obj,genotype,model)
            %INPUT:
            % genotype (str) = one of 'd16' or 'p95'
            % model (str) = one of '2C', '3C_partialObs', '3C_completeObs'


            switch genotype
                case 'd16'
                    locK = round(normrnd(obj.d16_Field_params.K,0.05*obj.d16_Field_params.K));
                    locF = round(normrnd(obj.d16_Field_params.F0,0.05*obj.d16_Field_params.F0));
                    %                    locr = obj.calcGrowthRate(locF,locK,obj.t_max_F_growth); % ASSUME MAX GROWTH RATE AT WEEK 5 (1 = 5 - 4)
                    locr = normrnd(obj.d16_Field_params.r,0.05*obj.d16_Field_params.r);
                    params = [locr,locK,locF,obj.t_start,obj.delta_t,obj.t_final,obj.ages_jitter'];
                    obs_data = obj.d16_counts;
                case 'p95'
                    locK = round(normrnd(obj.p95_Field_params.K,0.05*obj.p95_Field_params.K));
                    locF = round(normrnd(obj.p95_Field_params.F0,0.05*obj.p95_Field_params.F0));
                    locr = normrnd(obj.p95_Field_params.r,0.05*obj.p95_Field_params.r);
                    %locr = obj.calcGrowthRate(locF,locK,obj.t_max_F_growth); % ASSUME MAX GROWTH RATE AT WEEK 5 (1 = 5 - 4)
                    params = [locr,locK,locF,obj.t_start,obj.delta_t,obj.t_final,obj.ages_jitter'];
                    obs_data = obj.p95_counts;
            end
            switch model
                case '2C'
                    fun = @(theta,params)obj.LSQ_2compartments(theta,params);
                    xMin = [0,0];
                    x_guess = [0.01,0.1];
                    xMax = [1,1];
                case '3C_partialObs'
                    %passing in complete_obs_flag to function call
                    params = [false,params];
                    fun = @(theta,params)obj.LSQ_3compartments(theta,params);
                    xMin = [0,0,0];
                    x_guess = [0.001,0.5,0.1];
                    xMax = [1,100,3];
                case '3C_completeObs'
                    %passing in complete_obs_flag to function call
                    params = [true,params];
                    fun = @(theta,params)obj.LSQ_3compartments(theta,params);
                    xMin = [0,0,0];
                    x_guess = [0.001,0.5,0.1];
                    xMax = [1,100,3];
            end


            options = optimoptions('lsqcurvefit','Display','off');
            [param_estimate,sum_of_squares] = lsqcurvefit(fun,x_guess,params,obs_data,xMin,xMax,options);
            param_estimate = [locK,locF,locr,param_estimate];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LSQCURVEFIT input
        function pred_data = LSQ_2compartments(obj,theta,params)
            P.r = params(1);
            P.K = params(2);
            P.F0 = params(3);
            % to be estimated
            P.p1 = theta(1);
            P.p2 = theta(2);

            loc_t_start = params(4);
            loc_delta_t = params(5);
            loc_t_final = params(6);
            % shift observations to same coordinates as loc_t_start
            obs_t_steps = params(7:end) - loc_t_start;

            t_steps = loc_t_start:loc_delta_t:loc_t_final;

            modelState = obj.F2C(t_steps,P);


            obs_t_index = round(obs_t_steps./loc_delta_t);

            pred_data = modelState(obs_t_index,:);
        end % LSQ_2compartments
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function pred_data = LSQ_3compartments(obj,theta,params) %,complete_obs_flag)
            %INPUT:
            %params(end): complete_obs_flag (boolean) if true, use F3C_completeObs
            %   else, use F3C_partialObs
            complete_obs_flag = params(1);

            P.r = params(2);
            P.K = params(3);
            P.F0 = params(4);
            % to be estimated
            P.p0 = theta(1);
            P.p1 = theta(2);
            P.p2 = theta(3);

            loc_t_start = params(5);
            loc_delta_t = params(6);
            loc_t_final = params(7);
            % shift observations to same coordinates as loc_t_start
            obs_t_steps = params(8:end)-loc_t_start;

            t_steps = loc_t_start:loc_delta_t:loc_t_final;

            if complete_obs_flag
                modelState = obj.F3C_completeObs(t_steps,P);
            else
                modelState = obj.F3C_partialObs(t_steps,P);
            end

            obs_t_index = round(obs_t_steps./loc_delta_t);

            pred_data = modelState(obs_t_index,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DYNAMICS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function S = F2C(obj,t,P)
            % F2C = Field y(1), 2 compartments C1 y(2), C1 y(3)
            % INPUT:
            % P (struct) with fields
            % .F0 -- initial Field size
            % .r -- intrinsic growth rate of Field
            % .K -- carrying capacity of Field
            % .p1 -- transition rate from F to C1
            % .p2 -- transition rate from C1 to C2
            y0 = [P.F0, 0, 0];
            [~,Sv] = ode15s(@(t,y) growth(t,y,P),t,y0);
            function dydt = growth(t,y,P)

                dydt = zeros(3, 1);
                dydt(1) = P.r*y(1)*(1-y(1)/(P.K))-P.p1*y(1); %F
                dydt(2) = P.p1*y(1) - P.p2*y(2); %Small
                dydt(3) = P.p2*y(2);             %Large
            end
            % only return observed variables
            S = Sv(:,2:3);
        end % method F2C
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function S = F2C_all(obj,t,P)
            % F2C = Field y(1), 2 compartments C1 y(2), C1 y(3)
            % INPUT:
            % P (struct) with fields
            % .F0 -- initial Field size
            % .r -- intrinsic growth rate of Field
            % .K -- carrying capacity of Field
            % .p1 -- transition rate from F to C1
            % .p2 -- transition rate from C1 to C2
            y0 = [P.F0, 0, 0];
            [~,S] = ode15s(@(t,y) growth(t,y,P),t,y0);
            function dydt = growth(t,y,P)

                dydt = zeros(3, 1);
                dydt(1) = P.r*y(1)*(1-y(1)/(P.K))-P.p1*y(1); %F
                dydt(2) = P.p1*y(1) - P.p2*y(2); %Small
                dydt(3) = P.p2*y(2);             %Large
            end

        end % method F2C_all
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function S = F3C_partialObs(obj,t,P)
            % F3C = Field y(1), 3 compartments C1 y(2), C2 y(3)
            % INPUT:
            % P (struct) with fields
            % .F0 -- initial Field size
            % .r -- intrinsic growth rate of Field
            % .K -- carrying capacity of Field
            % .p0 -- transition rate from F to C1
            % .p1 -- transition rate from C1 to C2
            % .p2 -- transition rate from C2 to C3
            % OUTPUT: partial observations, only compartments C2 and C3

            y0 = [P.F0, 0, 0, 0];
            [~,Sv] = ode15s(@(t,y) growth(t,y,P),t,y0);
            function dydt = growth(t,y,P)

                dydt = zeros(4, 1);
                dydt(1) = P.r*y(1)*(1-y(1)/(P.K))-P.p0*y(1);
                dydt(2) = P.p0*y(1) - P.p1*y(2);
                dydt(3) = P.p1*y(2) - P.p2*y(3);
                dydt(4) = P.p2*y(3);
            end
            S = Sv(:,3:4);
        end % function F3C_partialObs

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function S = F3C_completeObs(obj,t,P)
            % F3C = Field y(1), 3 compartments C1 y(2), C2 y(3)
            % INPUT:
            % P (struct) with fields
            % .F0 -- initial Field size
            % .r -- intrinsic growth rate of Field
            % .K -- carrying capacity of Field
            % .p0 -- transition rate from F to C1
            % .p1 -- transition rate from C1 to C2
            % .p2 -- transition rate from C2 to C3
            % OUTPUT: complete observations compartmentsC1, C2 and C3

            y0 = [P.F0, 0, 0, 0];
            [~,Sv] = ode15s(@(t,y) growth(t,y,P),t,y0);
            function dydt = growth(t,y,P)

                dydt = zeros(4, 1);
                dydt(1) = P.r*y(1)*(1-y(1)/(P.K))-P.p0*y(1);
                dydt(2) = P.p0*y(1) - P.p1*y(2);
                dydt(3) = P.p1*y(2) - P.p2*y(3);
                dydt(4) = P.p2*y(3);
            end
            S = Sv(:,2:4);
        end % function F3C_completeObs

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Sv = F3C_all(obj,t,P)
            % F3C = Field y(1), 3 compartments C1 y(2), C2 y(3)
            % INPUT:
            % P (struct) with fields
            % .F0 -- initial Field size
            % .r -- intrinsic growth rate of Field
            % .K -- carrying capacity of Field
            % .p0 -- transition rate from F to C1
            % .p1 -- transition rate from C1 to C2
            % .p2 -- transition rate from C2 to C3
            % OUTPUT: complete observations compartmentsC1, C2 and C3

            y0 = [P.F0, 0, 0, 0];
            [~,Sv] = ode15s(@(t,y) growth(t,y,P),t,y0);
            function dydt = growth(t,y,P)

                dydt = zeros(4, 1);
                dydt(1) = P.r*y(1)*(1-y(1)/(P.K))-P.p0*y(1);
                dydt(2) = P.p0*y(1) - P.p1*y(2);
                dydt(3) = P.p1*y(2) - P.p2*y(3);
                dydt(4) = P.p2*y(3);
            end
           
        end % function F3C_completeObs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function r = calcGrowthRate(obj,F0,K,t)
        %     %INPUT:
        %     % F0 : (float) initial conditions
        %     % K : (float) carrying capacity
        %     % t : (float) time of max growth rate
        %
        %     ratio = F0 / (K - F0);
        %
        %     r = -log(ratio)/t;
        % end
    end % methods block
end % class