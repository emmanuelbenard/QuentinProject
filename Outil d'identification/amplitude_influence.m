run('datas.m')
Error_list=[];
for ii=1:length(Amp_list)
%% elevator
Ae = Amp_list(ii);
%% rudder
Ar = Amp_list(ii);
%% aileron
Aa = Amp_list(ii);
%% run the simulink model
open_system('model_plane_identification_script_v1')
sim('model_plane_identification_script_v1')

%%
    if (strcmp(coefficient,'Cl'))
    

            %Cl analysis
            Mx = Ix*x_dot(:,4)-Jxz*x_dot(:,6)+(Iz-Iy)*q.*r-Jxz*p.*q;
            Coeff_tot = Mx/(S*mean(yad1(:,3))*b);
            X = [linspace(1,1,length(time))' beta p r delta_a delta_r];
            theta_hat =inv(X'*X)*X'*Coeff_tot;

            % calculation of Cl coefficients
            Cl_0 = theta_hat(1);
            Cl_b = theta_hat(2);
            Cl_p = 2*V*theta_hat(3)/b;
            Cl_r = 2*V*theta_hat(4)/b;
            Cl_deltaa = theta_hat(5);
            Cl_deltar = theta_hat(6);

    elseif (strcmp(coefficient,'Cd'))

            % drag coefficients

            dynamic_pressure = mean(yad1(:,3));
            % calculation of Cz coefficient
            delta_az = -(Acc(:,3)-Acc(1,3)*linspace(1,1,length(time))');
            Cz_tot= delta_az*m/(S*dynamic_pressure);

            % calculation of Cx coefficient
            delta_ax = -m*(Acc(:,1)-Acc(1,1)*linspace(1,1,length(time))');
            Cx_tot= (delta_ax-thrust(:,1))/(S*dynamic_pressure);
            X = [linspace(1,1,length(time))' alpha q delta_e];
            Coeff_tot = -Cx_tot.*cos(alpha)-Cz_tot.*sin(alpha);
            theta_hat =inv(X'*X)*X'*Coeff_tot;

            % calculation of Cd coefficients
            Cd_0 = -theta_hat(1);
            Cd_alpha = theta_hat(2);
            Cd_q = 2*V*theta_hat(3)/cbar;
            Cd_deltae =  theta_hat(4);
        
              
        elseif (strcmp(coefficient,'CL'))
        
           % drag coefficients

            dynamic_pressure = mean(yad1(:,3));
            % calculation of Cz coefficient
            delta_az = -(Acc(:,3)-Acc(1,3)*linspace(1,1,length(time))');
            Cz_tot= delta_az*m/(S*dynamic_pressure);

            % calculation of Cx coefficient
            delta_ax = -m*(Acc(:,1)-Acc(1,1)*linspace(1,1,length(time))');
            Cx_tot= (delta_ax-thrust(:,1))/(S*dynamic_pressure);
            X = [linspace(1,1,length(time))' alpha q delta_e];
            Coeff_tot = -Cz_tot.*cos(alpha)+Cx_tot.*sin(alpha);
            theta_hat =inv(X'*X)*X'*Coeff_tot;

            % calculation of CL coefficients
            CL_0 = -theta_hat(1);
            CL_alpha = theta_hat(2);
            CL_q = 2*V*theta_hat(3)/cbar;
            CL_deltae = theta_hat(4);
        
        elseif (strcmp(coefficient,'Cm'))
            
            %Cm analysis
            My = Iy*x_dot(:,5)+(Ix-Iz)*r.*p -Jxz*(p.*p-r.*r);
            Coeff_tot = My/(S*mean(yad1(:,3))*cbar);
            X = [linspace(1,1,length(time))' alpha q delta_e];
            theta_hat =inv(X'*X)*X'*Coeff_tot;
            
            % calculation of Cm coefficients
            Cm_0 = -theta_hat(1);
            Cm_alpha = theta_hat(2);
            Cm_q = 2*V*theta_hat(3)/cbar;
            Cm_deltae = theta_hat(4);
            
        elseif (strcmp(coefficient,'Cn'))
            
            % Cn analysis
            Mz = Iz*x_dot(:,6)-Jxz*x_dot(:,4)+(Iy-Ix)*p.*q+Jxz*r.*q;
            Coeff_tot = Mz/(S*mean(yad1(:,3))*b);
            X = [linspace(1,1,length(time))' beta p r delta_a delta_r];
            theta_hat =inv(X'*X)*X'*Coeff_tot;
            
            % calculation of Cn coefficients
            Cn_0 = theta_hat(1);
            Cn_b = theta_hat(2);
            Cn_p = 2*V*theta_hat(3)/b;
            Cn_r = 2*V*theta_hat(4)/b;
            Cn_deltaa = theta_hat(5);
            Cn_deltar = theta_hat(6);

    elseif (strcmp(coefficient,'Cy'))
            
            % Cy analysis
            delta_ay = Acc(:,2)-Acc(1,2)*linspace(1,1,length(time))';
            Coeff_tot= delta_ay*m/(S*mean(yad1(:,3)));
            X = [linspace(1,1,length(time))' beta p r delta_a delta_r];
            theta_hat =inv(X'*X)*X'*Coeff_tot;
            theta_Cy = [Cy_0_true Cy_b_true Cy_p_true Cy_r_true Cy_deltaa_true Cy_deltar_true]';

            % calculation of Cy coefficients
            Cy_0 = theta_hat(1);
            Cy_b = theta_hat(2);
            Cy_p = 2*V*theta_hat(3)/b;
            Cy_r = 2*V*theta_hat(4)/b;
            Cy_deltaa = theta_hat(5);
            Cy_deltar = theta_hat(6);

        else 
        disp('Ce coefficient n''est pas disponible')
        break
                                
    end


[time nb_parameters] = size(X);
% matrix of prediction
K_coeff = X*inv(X'*X)*X';
I_coeff = eye(length(K_coeff));
% residuals 
residuals = (I_coeff - K_coeff)*Coeff_tot;

% measurement error variance
sigma2_hat = (residuals'*residuals)/(time-nb_parameters);

% fit_error = sqrt(sigma2_hat);

% The estimated parameter covariance matrix
cov_theta_hat = sigma2_hat*inv(X'*X);

% The standard errors of the estimated parameters
s_theta_hat = sqrt(diag(cov_theta_hat));

% residuals sum of squares:
SSe = Coeff_tot'*Coeff_tot-theta_hat'*X'*Coeff_tot;

% regression sum of sqares:
SSr = theta_hat'*X'*Coeff_tot-time*(mean(Coeff_tot))^2;

% total sum of squares
SSt = Coeff_tot'*Coeff_tot-time*mean(Coeff_tot)^2;

% coefficient of determination
R2  = 100*SSr/SSt; % indicator of the proportion of the measured output that can be explained by the model



% t paramter
t_parameter = abs(theta_hat./s_theta_hat);

% The estimated parameter correlation matrix
cor_theta_hat = (diag(1./s_theta_hat))*cov_theta_hat*(diag(1./s_theta_hat));


    if (strcmp(coefficient,'Cl'))
    
            paramter_Cl = {'Parameters'          'true'                   'estimated'         's parameter'                   't parameter'
                           'Cl_0_true'           num2str(Cl_0_true)       num2str(Cl_0)       num2str(s_theta_hat(1))      num2str(t_parameter(1))
                           'Cl_b_true'           num2str(Cl_b_true)       num2str(Cl_b)       num2str(s_theta_hat(2))      num2str(t_parameter(2))              
                           'Cl_p_true'           num2str(Cl_p_true)       num2str(Cl_p)       num2str(s_theta_hat(3))      num2str(t_parameter(3))   
                           'Cl_r_true'           num2str(Cl_r_true)       num2str(Cl_r)       num2str(s_theta_hat(4))      num2str(t_parameter(4))
                           'Cl_deltaa_true'      num2str(Cl_deltaa_true)  num2str(Cl_deltaa)  num2str(s_theta_hat(5))      num2str(t_parameter(5))
                           'Cl_deltar_true'      num2str(Cl_deltar_true)  num2str(Cl_deltar)  num2str(s_theta_hat(6))      num2str(t_parameter(6))
                           'coefficient of determination : ' num2str(R2)  ''                  ''                              ''
                            };

                % purcent of error
                ECl_beta = 100*abs((Cl_b_true-Cl_b)/Cl_b_true);
                ECl_p = 100*abs((Cl_p_true-Cl_p)/Cl_p_true);
                ECl_r = 100*abs((Cl_r_true-Cl_r)/Cl_r_true);
                ECl_deltaa = 100*abs((Cl_deltaa_true-Cl_deltaa)/Cl_deltaa_true);
                ECl_deltar = 100*abs((Cl_deltar_true-Cl_deltar)/Cl_deltar_true);
                ECl = [ECl_beta ECl_p ECl_r ECl_deltaa ECl_deltar];
                theta_hat_coefficient = {'Cl_{beta}','Cl_{p}','Cl_{r}','Cl_{deltaa}', 'Cl_{deltar}'};
                Error_list=[Error_list; ECl];

    elseif (strcmp(coefficient,'Cd'))

            paramter_Cd = {'Parameters'   'true'                   'estimated'         's parameter'                   't parameter'
                           'Cd_0'         num2str(Cd_0_true)       num2str(Cd_0)       num2str(s_theta_hat(1))      num2str(t_parameter(1))
                           'Cd_alpha'     num2str(Cd_alpha_true)   num2str(Cd_alpha)   num2str(s_theta_hat(2))      num2str(t_parameter(2))              
                           'Cd_q'         num2str(Cd_q_true)       num2str(Cd_q)       num2str(s_theta_hat(3))      num2str(t_parameter(3))   
                           'Cd_deltae'    num2str(Cd_deltae_true)  num2str(Cd_deltae)  num2str(s_theta_hat(4))      num2str(t_parameter(4))           
                           'coefficient of determination : ' num2str(R2) '' '' ''
                           };

                % purcent of error
                ECd_0 = 100*abs(Cd_0_true-Cd_0)/Cd_0_true;
                ECd_alpha = 100*abs(Cd_alpha_true-Cd_alpha)/Cd_alpha_true;
                ECd_q = 100*abs(Cd_q_true-Cd_q)/Cd_q_true;
                ECd_deltae = 100*abs(Cd_deltae_true-Cd_deltae)/Cd_deltae_true;
                ECd = [ECd_0 ECd_alpha];
                theta_hat_coefficient = {'Cd_{0}', 'Cd_{alpha}'};
                Error_list=[Error_list; ECd];

    elseif (strcmp(coefficient,'CL'))
            
            paramter_CL = {'Parameters'   'true'                   'estimated'         's parameter'                   't parameter'
                 'CL_0'         num2str(CL_0_true)       num2str(CL_0)       num2str(s_theta_hat(1))      num2str(t_parameter(1))
                 'CL_alpha'     num2str(CL_alpha_true)   num2str(CL_alpha)   num2str(s_theta_hat(2))      num2str(t_parameter(2))              
                 'CL_q'         num2str(CL_q_true)       num2str(CL_q)       num2str(s_theta_hat(3))      num2str(t_parameter(3))   
                 'CL_deltae'    num2str(CL_deltae_true)  num2str(CL_deltae)  num2str(s_theta_hat(4))      num2str(t_parameter(4))
                 'coefficient of determination : ' num2str(R2) '' '' ''
                 };
 
                % purcent of error
                ECL_0 = 100*abs(CL_0_true-CL_0)/CL_0_true;
                ECL_alpha = 100*abs(CL_alpha_true-CL_alpha)/CL_alpha_true;
                ECL_q = 100*abs(CL_q_true-CL_q)/CL_q_true;
                ECL_deltae = 100*abs(CL_deltae_true-CL_deltae)/CL_deltae_true;
                ECL = [ECL_0 ECL_alpha ECL_q ECL_deltae];
                theta_hat_coefficient = {'CL_{0}', 'CL_{alpha}', 'CL_{q}', 'CL_{deltae}'};
                Error_list=[Error_list; ECL];
            
    elseif (strcmp(coefficient,'Cm'))
                
                paramter_Cm = {'Parameters'   'true'                   'estimated'         's parameter'                   't parameter'
                 'Cm_0'         num2str(Cm_0_true)       num2str(Cm_0)       num2str(s_theta_hat(1))      num2str(t_parameter(1))
                 'Cm_alpha'     num2str(Cm_alpha_true)   num2str(Cm_alpha)   num2str(s_theta_hat(2))      num2str(t_parameter(2))              
                 'Cm_q'         num2str(Cm_q_true)       num2str(Cm_q)       num2str(s_theta_hat(3))      num2str(t_parameter(3))   
                 'Cm_deltae'    num2str(Cm_deltae_true)  num2str(Cm_deltae)  num2str(s_theta_hat(4))      num2str(t_parameter(4))
                 'coefficient of determination : ' num2str(R2) '' '' ''
                  };
 
                % purcent of error
                ECm_alpha = 100*abs((Cm_alpha_true-Cm_alpha)/Cm_alpha_true);
                ECm_q = 100*abs((Cm_q_true-Cm_q)/Cm_q_true);
                ECm_deltae = 100*abs((Cm_deltae_true-Cm_deltae)/Cm_deltae_true);
                ECm = [ECm_alpha ECm_q ECm_deltae];
                theta_hat_coefficient = {'Cm_{alpha}', 'Cm_{q}', 'Cm_{deltae}'};
                Error_list=[Error_list; ECm];
            
    elseif (strcmp(coefficient,'Cn'))
                 
               paramter_Cn = {'Parameters'          'true'                   'estimated'         's parameter'                   't parameter'
               'Cn_0_true'           num2str(Cn_0_true)       num2str(Cn_0)       num2str(s_theta_hat(1))      num2str(t_parameter(1))
               'Cn_b_true'           num2str(Cn_b_true)       num2str(Cn_b)       num2str(s_theta_hat(2))      num2str(t_parameter(2))              
               'Cn_p_true'           num2str(Cn_p_true)       num2str(Cn_p)       num2str(s_theta_hat(3))      num2str(t_parameter(3))   
               'Cn_r_true'           num2str(Cn_r_true)       num2str(Cn_r)       num2str(s_theta_hat(4))      num2str(t_parameter(4))
               'Cn_deltaa_true'      num2str(Cn_deltaa_true)  num2str(Cn_deltaa)  num2str(s_theta_hat(5))      num2str(t_parameter(5))
               'Cn_deltar_true'      num2str(Cn_deltar_true)  num2str(Cn_deltar)  num2str(s_theta_hat(6))      num2str(t_parameter(6))
               'coefficient of determination : ' num2str(R2)  ''                  ''                              ''
                };
            
                % purcent of error
                ECn_beta = 100*abs((Cn_b_true-Cn_b)/Cn_b_true);
                ECn_p = 100*abs((Cn_p_true-Cn_p)/Cn_p_true);
                ECn_r = 100*abs((Cn_r_true-Cn_r)/Cn_r_true);
                ECn_deltaa = 100*abs((Cn_deltaa_true-Cn_deltaa)/Cn_deltaa_true);
                ECn_deltar = 100*abs((Cn_deltar_true-Cn_deltar)/Cn_deltar_true);
                ECn = [ECn_beta ECn_p ECn_r ECn_deltaa ECn_deltar];
                theta_hat_coefficient = {'Cn_{beta}','Cn_{p}','Cn_{r}','Cn_{deltaa}', 'Cn_{deltar}'};
                Error_list=[Error_list; ECn];
      
    elseif (strcmp(coefficient,'Cy'))
                                 
                                 
                paramter_Cy = {'Parameters'          'true'                   'estimated'         's parameter'                   't parameter'
               'Cy_0_true'           num2str(Cy_0_true)       num2str(Cy_0)       num2str(s_theta_hat(1))      num2str(t_parameter(1))
               'Cy_b_true'           num2str(Cy_b_true)       num2str(Cy_b)       num2str(s_theta_hat(2))      num2str(t_parameter(2))              
               'Cy_p_true'           num2str(Cy_p_true)       num2str(Cy_p)       num2str(s_theta_hat(3))      num2str(t_parameter(3))   
               'Cy_r_true'           num2str(Cy_r_true)       num2str(Cy_r)       num2str(s_theta_hat(4))      num2str(t_parameter(4))
               'Cy_deltaa_true'      num2str(Cy_deltaa_true)  num2str(Cy_deltaa)  num2str(s_theta_hat(5))      num2str(t_parameter(5))
               'Cy_deltar_true'      num2str(Cy_deltar_true)  num2str(Cy_deltar)  num2str(s_theta_hat(6))      num2str(t_parameter(6))
               'coefficient of determination : ' num2str(R2)  ''                  ''                              ''
                 };
             
                % purcent of error
                ECy_b = 100*abs((Cy_b_true-Cy_b)/Cy_b_true);
                ECy_deltar = 100*abs((Cy_deltar_true-Cy_deltar)/Cy_deltar_true);
                ECy = [ECy_b ECy_deltar];
                theta_hat_coefficient = {'Cy_{b}', 'Cy_{deltar}'};
                Error_list=[Error_list; ECy];
                
 
    end
end

%% plot of the influence of the amplitude
figure(1)
title(strcat('influence de l''amplitude sur l''erreur relative pour le coeffcient: ', coefficient))
xlabel('Amplitude (�)')
ylabel('Erreur relative')
for i=1:length(theta_hat_coefficient)
    hold on 
    plot(Amp_list,Error_list(:,i),'Displayname',theta_hat_coefficient{i})
end
legend('show')
