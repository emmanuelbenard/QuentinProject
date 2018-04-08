run('datas.m')
% run 1
% elevator
Ae = Amp_elevator1;
test_elevator = elevator1;
% rudder
Ar = Amp_rudder1;
test_rudder = rudder1;
% aileron
Aa = Amp_aileron1;
test_aileron = aileron1;
%% run of the first simulation 
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
          return                      
    end

%% Run 2
% elevator
Ae = Amp_elevator2;
test_elevator = elevator2;
% rudder
Ar = Amp_rudder2;
test_rudder = rudder2;
% aileron
Aa = Amp_aileron2;
test_aileron = aileron2;
%% run of the second simulation 
open_system('model_plane_identification_script_v1')
sim('model_plane_identification_script_v1')

    if (strcmp(coefficient,'Cl'))
        Mx = Ix*x_dot(:,4)-Jxz*x_dot(:,6)+(Iz-Iy)*q.*r-Jxz*p.*q;
        C_tot_output = Mx/(S*mean(yad1(:,3))*b);

        % Cl_tot calculated with the estimated coefficient of the last part
        C_tot_model = Cl_0 + Cl_b*beta + Cl_p*p*b/(2*V) + Cl_r*r*b/(2*V) + Cl_deltaa * delta_a + Cl_deltar * delta_r;
        
    elseif (strcmp(coefficient,'Cd'))
        % drag coefficients
        dynamic_pressure = mean(yad1(:,3));
        delta_az = -(Acc(:,3)-Acc(1,3)*linspace(1,1,length(time))');
        Cz_tot= delta_az*m/(S*dynamic_pressure);
        delta_ax = -m*(Acc(:,1)-Acc(1,1)*linspace(1,1,length(time))');
        Cx_tot= (delta_ax-thrust(:,1))/(S*dynamic_pressure);
        C_tot_output = -Cx_tot.*cos(alpha)-Cz_tot.*sin(alpha);
        % Cd_tot calculated with the estimated coefficient of the last part
        C_tot_model =-Cd_0 + Cd_alpha*alpha + Cd_q*q*cbar/(2*V) + Cd_deltae * delta_e;
        
    elseif (strcmp(coefficient,'Cm'))
        
        My = Iy*x_dot(:,5)+(Ix-Iz)*r.*p -Jxz*(p.*p-r.*r);
        % Cm_tot fom datas
        C_tot_output = My/(S*mean(yad1(:,3))*cbar);
        % Cm_tot calculated with the estimated coefficient of the last part
        C_tot_model = Cm_0 + Cm_alpha*alpha + Cm_q*q*cbar/(2*V) + Cm_deltae * delta_e;
        
    elseif (strcmp(coefficient,'Cn'))

        Mz = Iz*x_dot(:,6)-Jxz*x_dot(:,4)+(Iy-Ix)*p.*q+Jxz*r.*q;
        C_tot_output = Mz/(S*mean(yad1(:,3))*b);
        % Cn_tot calculated with the estimated coefficient of the last part
        C_tot_model = Cn_0 + Cn_b*beta + Cn_p*p*b/(2*V) + Cn_r*r*b/(2*V) + Cn_deltaa * delta_a + Cn_deltar * delta_r;

    elseif (strcmp(coefficient,'CL'))
        
        dynamic_pressure = mean(yad1(:,3));
        delta_az = (Acc(:,3)-Acc(1,3)*linspace(1,1,length(time))');
        Cz_tot= delta_az*m/(S*dynamic_pressure);
        delta_ax = -m*(Acc(:,1)-Acc(1,1)*linspace(1,1,length(time))');
        Cx_tot= (delta_ax-thrust(:,1))/(S*dynamic_pressure);
        C_tot_output = -Cz_tot.*cos(alpha)+Cx_tot.*sin(alpha);
        % CL_tot calculated with the estimated coefficient of the last part
        C_tot_model = -CL_0 + CL_alpha*alpha + CL_q*q*cbar/(2*V) + CL_deltae * delta_e;

    elseif (strcmp(coefficient,'Cy'))

        delta_ay = Acc(:,2)-Acc(1,2)*linspace(1,1,length(time))';
        C_tot_output= delta_ay*m/(S*mean(yad1(:,3)));

        % Cy_tot calculated with the estimated coefficient of the last part
        C_tot_model = Cy_0 + Cy_b*beta + Cy_p*p*b/(2*V) + Cy_r*r*b/(2*V) + Cy_deltaa * delta_a + Cy_deltar * delta_r;

    else 
        disp('Ce coefficient n''est pas disponible')
          return                      
    
    end
figure
plot(time, C_tot_output,time,C_tot_model)
title('Comparaison entre la mesure et le modèle')
xlabel('temps (s)')
ylabel(strcat('coeffcient:',coefficient))
legend('Mesures','Modèle')
