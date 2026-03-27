clear; close all;

%% Chemin de sauvegarde des figures
fpath = 'C:\Users\moust\Documents\MATLAB\Tapha_Simu\Graphe2\Graphe3';
LETTERS='A':'Z'; 
%% Paramètres généraux
xmin = -3; xmax =3; Nx = 5 ; x = linspace(xmin, xmax, Nx);
AgeMax = 190*52; da = 1; VectAge = 0:da:AgeMax; Na = length(VectAge);
dt = 1; T = 170*52; 

time = 0:dt:T; Nt = length(time); 
scenario = 'FV';
%% Paramètres biologiques
a_0 = 10; alpha = 0.1; l_0 = 5; g = 471512.612; y = 0.3; kappa = 0.005;
b = 213.032; c = 2.535; coef =0.5*((1/35)*(log(1e-3)/log(1-0.05)));
%Initialgraine = 45*exp(-(x+0.04).^2 / (0.06/5)^2);
Initialgraine = 25*ones(1,Nx);
%% Génération de 3 vecteurs d_0 différents
d0_min = 1; d0_max = 10; d0_vectors = zeros(3, Nx);
% Valeurs de l'étendue du patch
mean_dispersion = b*sqrt(pi)/2*gamma(c-3/2)/gamma(c-1);
d0_values =  3; %[0.5 1 3 10];
qmax_list = [ 0.30 0.6 0.75];

%% Boucle sur chaque vecteur d_0
%for d_0 = d0_values
    k0 = (c-1)/(pi*b^2).*(1+(mean_dispersion.^2)/(b^2)).^(-c);
    %Q = exp(-x.^2/d_0);
    %% Initialisation des variables
    s = zeros(Nt, Nx);
    n = zeros(Nt, Na, Nx);
    s(1,:) = Initialgraine;
    %peupler qlques arbres adultes dès le départ pour fécondité immédiate
    for ix = 1:Nx
      for ia = 1:Na
         if VectAge(ia) >= a_0 && VectAge(ia) <= a_0+5
           n(1,ia,ix)   = coef;
         end
      end
    end
    %% Calcul de la fécondité phi(a)
  % q = Q;
    for i = 1:length(qmax_list)
     qmax = qmax_list(i);
        q = qmax*ones(1, Nx); 
      phi = zeros(Nx, Na);
    for a = 1:Na
        for ix = 1:Nx
          if VectAge(a) < a_0
            phi(ix,a) = 0;
          else
            phi(ix,a)= q(ix)*g*alpha*(alpha*(VectAge(a)-a_0))^(l_0-1)*...
                exp(-alpha*(VectAge(a)-a_0))/gamma(l_0);
          end
        end
    end
%% Fonctions temporaires
    tau = 52;         % période (ex : 12 mois ou 1 an)
    %A_beta  = 0.8;   % amplitude saisonnière
    %A_delta = 0.3;
    %beta0  = 0.1*(1+q) ; delta0 = 0.03*(1-q);
    beta = 0.8*(1+q); delta = 0.3*(1-q); r= 0.6*(1+q);
    mu_vec = zeros(1,Na); 
    mu_vec(VectAge <=180) = 0.005928; mu_vec(VectAge > 180) = 1000;
    mu = zeros(Na,Nx);
    for ia=1:Na
       mu(ia,:) = 2*mu_vec(ia)*(1-q); 
    end
    Mu =sum(mu,2); 
    
    for t = 1:(Nt-1)
      %beta_t  = beta0.*(1 + A_beta*sin(2*pi*time(t)/tau));
      %delta_t = delta0.*(1 + A_delta*sin(2*pi*time(t)/tau));
      %r_t     = 0.6.*(1+0.75*sin(2*pi*time(t)/tau)) ;
      %Intrn   = k0.*r_t.*trapz(VectAge, squeeze(n(t,:,:)).*phi(:,:).', 1); 
       Intrn   = k0.*r.*trapz(VectAge, squeeze(n(t,:,:)).*phi(:,:).', 1);
      %s(t+1,:) =  (Intrn + s(t,:)/dt)/(1/dt + delta_t) ;
      s(t+1,:) =  (Intrn + s(t,:)/dt)/(1/dt+ delta) ; 
        for a = 2:(Na-1)
            for ix = 1:Nx
               Fonctionr = trapz(VectAge, n(t,:,ix));
               H_val = max(1-kappa*Fonctionr, 0);
               %n(t+1,1,ix) = H_val*beta_t(ix) * s(t,ix);
               n(t+1,1,ix) =  H_val.*beta(ix)*s(t,ix);
               n(t+1,a,ix) =  (n(t,a,ix)/dt + n(t+1,a-1,ix)/da)/(1/dt +...
                   1/da + mu(a,ix));
            end
        end
    end  
    %% Calcul des totaux pour le tracé
    nxtPlot=zeros(Nt,Nx); 
    for t = 1:Nt
       for ix = 1:Nx
         nxtPlot(t,ix) = trapz(VectAge, n(t,:,ix));
       end
    end
    ntPlot = sum(nxtPlot,2);
    tsPlot = sum(s,2);
    %% Plot probabilité de survie
    Pi = zeros(Na,Nx);
    for a=1:Na
        Pi(a,:) = exp(-trapz(VectAge(1:a), mu(1:a,:)));
    end
    plotPi =  sum(Pi,2); 
    %% Plot fécondité intégrée 
    tphi = sum(phi,1);
    qphi=zeros (1, Nx) ;
    for ix= 1:Nx 
      qphi(ix) = trapz(VectAge, phi(ix,:).*Pi(:,ix)');
    end  
    %% Plot dynamique graines et arbres
    Pbar = mean(Pi, 2);
    nbar = mean(nxtPlot,2);
    sbar = mean(s,2);
    figure; 
    FigName = ['DynamicTreeSeedqmaarr170s_', num2str(qmax),'.jpeg'];
    %FigName= ['DynamicTreeSeeddt100_',   num2str(d_0),'.jpeg'];
    subplot(3,2,1)
    hold on;
    plot(time, tsPlot, 'r-', 'LineWidth', 1.6);
    plot(time, ntPlot, 'b-', 'LineWidth', 1.6);
    grid on        
    xlabel('Time (Weak)'); %ylabel('Density');
    %title(['q_{max}=', num2str(qmax)]);
    %title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(1) ')  }',...
    %'$ \quad d_0 = \quad $', num2str(d_0)],'interpreter','latex');
    title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(1) ')}',...
    '$\quad q_{max} = \quad $', num2str(qmax)],'interpreter','latex');
    legend('TS','TT', 'Location', 'eastoutside');
    subplot(3,2,2)
    hold on;
    plot(time, sbar, 'r-', 'LineWidth', 1.6);
    plot(time, nbar, 'b','LineWidth', 1.6) ;
    grid on        
    xlabel('Time (Weak)'); %ylabel('Density');
    title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(2) ')}'], ...
    'interpreter','latex');
    legend('MS','MT','Location', 'eastoutside');
    subplot(3,2,3)
    plot(VectAge, tphi,'g-','LineWidth',1.6);
    grid on        
    xlabel('Age (Weak)'); ylabel('$\psi(a)$','interpreter','latex'); 
    %title(['Fecundity function - d_0 = ', num2str(d_0(i_d))]);
    title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(3) ')}'], ...
    'interpreter','latex');
    subplot(3,2,4)
    plot(VectAge, Mu,'r-','LineWidth',1.6);
    grid on 
    xlabel('Age (Weak)'); ylabel('$\mu(a)$','interpreter','latex'); 
    %title('Mortality rate');
    title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(4) ')  }'],...
    'interpreter','latex');
    subplot(3,2,5)
    plot(VectAge, plotPi,'b-','LineWidth',1.6);
    grid on 
    xlabel('Age (Weak)','fontsize',11); 
    ylabel('$ \pi(a)$','interpreter','latex','fontsize',11); 
    %title('Survival probability');
    title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(5) ')}'],...
    'interpreter','latex');
    subplot(3,2,6)
    plot(VectAge, Pbar,'b-','LineWidth',1.6);
    grid on 
    xlabel('Age (Weak)','fontsize',11); 
    ylabel('$  \pi_m(a)$','interpreter','latex','fontsize',11); 
    %title('Survival probability');
    title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(6) ')}'],...
    'interpreter','latex');
    end 
   saveas(gcf,fullfile(fpath, FigName));
 %end



