function inversion

global mag2data

% Inversion of magnetic data
%
% reference "Compact Gravity Inversion" by Last and Kubick, Geophysics, 48,6,713-721 (1983)
%
% Written by Stefano Stocco on September, 2007
% Modified by Stefano Stocco on December, 2008

% Copyright (C) 2007 Stefano Stocco
% 
% This file is part of MAG2DATA.
% 
% MAG2DATA is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program; if not, see <http://www.gnu.org/licenses>.

Max_iter=mag2data.Max_iter;
l_0=mag2data.l_0;
target=mag2data.target/(4*pi);
F_obs_int=mag2data.F_obs_int;
kernelequalization=get(findobj('Tag','kernelequalization'),'Value');

if mag2data.controls.initialmodel==1    % se ho un modello iniziale
    S=mag2data.S;     % vettore delle suscettività derivante da un modello di primo tentativo
end

if mag2data.controls.campo_tot==1   % procede all'inversione considerando Fobs come campo totale
    A=mag2data.Ai;
else
    if mag2data.controls.grad==1
        A=mag2data.Ai_s;     % procede all'inversione considerando Fobs come gradiente verticale
    else
        A=mag2data.A;
    end
end

[m,n]=size(A);
I_n=eye(n);          % matrice identità n x n
I_m=eye(m);          % matrice identità m x m

eps_sus=[10^(-10) 10^(-11) 10^(-12) 10^(-13)];          % perturbazioni della funzione peso suscettività                         
S_iter=zeros(n,Max_iter,length(eps_sus));       % S_iter è una matrice tridimensionale
Fcalcolato=zeros(length(F_obs_int),Max_iter,length(eps_sus));   % Fcalcolato è una matrice tridimensionale
L2norm=zeros(length(eps_sus),Max_iter);
L2normuti=zeros(length(eps_sus),Max_iter);
weighted_resolution=zeros(n,Max_iter,length(eps_sus));

switch kernelequalization
    case 0
        Ws_inv=I_n;          % funzione peso suscettività per la prima iterazione
        D=(A*Ws_inv*A').*I_m;
        We_inv=l_0^(2)*D;    % funzione peso errori per la prima iterazione 
        K=1;
        err=l_0;
        
        if mag2data.controls.initialmodel==0    % se non ho un modello iniziale
            S=Ws_inv*A'*inv(A*Ws_inv*A'+We_inv)*F_obs_int;     % vettore delle suscettività calcolato per la prima iterazione (minimi quadrati)
        end
        
        fprintf('\nPerforming inversion...\n');
        tic;
        wait = waitbar(0,'Please wait... inversion computing');
        for l=1:length(eps_sus)              % primo ciclo = calcolo per diversi eps_sus

                waitbar(l/length(eps_sus), wait)
                set(wait,'Color',[0.5235    0.5785    0.7471],'Units','normalized','Position',[0.32 0.25 0.36 0.1])

            wait_iter = waitbar(0,'iteration');
            for k=1:Max_iter                % secondo ciclo: è il processo di inversione vero e proprio

                waitbar(k/Max_iter, wait_iter)
                set(wait_iter,'Color',[0.5235    0.5785    0.7471])

                cons=S./target;        
                H(find(cons<=1))=0;         % Heaviside step function
                H(find(cons>1))=1;          % constraint superiore (inferiore): valore dell'anomalia <= (>=) target
                Fobs_red=F_obs_int-target.*A*H';
                Ws_inv_red=((diag(S)).^(2)).*(I_n-diag(H))+(eps_sus(l))*I_n;    % funzione peso suscettività iterazioni successive alla prima

                if K==1
                    l_0=l_0;  % per la prima iterazione uso l_0 dato in input
                else
                    l_0=l_0*((max(abs(Fobs_red0-err)))/(max(abs(Fobs_red-err))));  % l_0 cambia ad ogni iterazione per le correzioni apportate
                end

                D=(A*Ws_inv_red*A').*I_m;
                We_inv_red=l_0^(2)*D;               % funzione peso errori iterazioni successive alla prima
                G_g=Ws_inv_red*A'*inv(A*Ws_inv_red*A'+We_inv_red);
                
                S=G_g*Fobs_red+target.*H';  % vettore suscettività iterazioni successive alla prima
                
%                 G_g=A'*((A*A'+We_inv_red)^(-1));
                res=G_g*A;            % model resolution matrix
                R=diag(res);    % diagonale della model resolution matrix (risoluzione del modello)
                weighted_resolution(:,k,l)=R;
                                
                if target>0                 % caso contrasto di suscettività positivo
                    cons=find(S<0);         % constraint inferiore
                else                        % caso contrasto di suscettività negativo
                    cons=find(S>0);         % constraint superiore
                end
                S(cons)=0;                  % se il contrasto è positivo azzero le celle con S negativo (e viceversa)
                Fobs_red0=Fobs_red;
                S_iter(:,k,l)=S;     % la matrice S_iter è 3D: i vettori delle suscettività vengono inseriti in colonne, 
                                     % e ogni colonna corrisponde ad un'iterazione;
                                     % la terza dimensione serve perché l'intero
                                     % ciclo viene computato per ogni eps_sus

                Fcal=A*S;          % forward per ogni inversione
                Fcalcolato(:,k,l)=Fcal;     % il forward viene posizionato nella matrice 3D
                residual=(F_obs_int-Fcal);    % residual
                L2norm(l,k)=norm(residual);        % norma2 su tutto il segnale
                siguti=find(mag2data.sig==1);              % calcolo la norma2 riferita alla parte utile del segnale
                resuti=(F_obs_int(siguti)-Fcal(siguti));
                L2normuti(l,k)=norm(resuti);  % norma2 della parte di segnale utile
            end
            close(wait_iter)
        end
        close(wait)
        toc
    case 1
        Wz=mag2data.Wz;
        Ws_inv=I_n;          % funzione peso suscettività per la prima iterazione
        D=(A*Ws_inv*Wz*A').*I_m;
        We_inv=l_0^(2)*D;    % funzione peso errori per la prima iterazione 
        K=1;
        err=l_0;
        
        if mag2data.controls.initialmodel==0
            S=Ws_inv*Wz*A'*inv(A*Ws_inv*Wz*A'+We_inv)*F_obs_int;     % vettore delle suscettività calcolato per la prima iterazione (minimi quadrati)
        end

        eps_sus=[10^(-10) 10^(-11) 10^(-12) 10^(-13)];          % perturbazioni della funzione peso suscettività                         
        S_iter=zeros(n,Max_iter,length(eps_sus));       % S_iter è una matrice tridimensionale
        Fcalcolato=zeros(length(F_obs_int),Max_iter,length(eps_sus));   % Fcalcolato è una matrice tridimensionale
        L2norm=zeros(length(eps_sus),Max_iter);
        L2normuti=zeros(length(eps_sus),Max_iter);
        
        fprintf('\nPerforming inversion...\n');
        tic;
        wait = waitbar(0,'Please wait... inversion computing');
        for l=1:length(eps_sus)              % primo ciclo = calcolo per diversi eps_sus

                waitbar(l/length(eps_sus), wait)
                set(wait,'Color',[0.5235    0.5785    0.7471],'Units','normalized','Position',[0.32 0.25 0.36 0.1])

            wait_iter = waitbar(0,'iteration');
            for k=1:Max_iter                % secondo ciclo: è il processo di inversione vero e proprio

                waitbar(k/Max_iter, wait_iter)
                set(wait_iter,'Color',[0.5235    0.5785    0.7471])

                cons=S./target;        
                H(find(cons<=1))=0;         % Heaviside step function
                H(find(cons>1))=1;          % constraint superiore (inferiore): valore dell'anomalia <= (>=) target
                Fobs_red=F_obs_int-target.*A*H';
                Ws_inv_red=((diag(S)).^(2)).*(I_n-diag(H))+(eps_sus(l))*I_n;    % funzione peso suscettività iterazioni successive alla prima

                if K==1
                    l_0=l_0;  % per la prima iterazione uso l_0 dato in input
                else
                    l_0=l_0*((max(abs(Fobs_red0-err)))/(max(abs(Fobs_red-err))));  % l_0 cambia ad ogni iterazione per le correzioni apportate
                end

                D=(A*Ws_inv_red*Wz*A').*I_m;
                We_inv_red=l_0^(2)*D;               % funzione peso errori iterazioni successive alla prima
                G_g=Ws_inv_red*Wz*A'*inv(A*Ws_inv_red*Wz*A'+We_inv_red);
                
                S=G_g*Fobs_red+target.*H';  % vettore suscettività iterazioni successice alla prima
                
%                 G_g=A'*((A*A'+We_inv_red)^(-1));
                res=G_g*A;            % model resolution matrix
                R=diag(res);    % diagonale della model resolution matrix (risoluzione del modello)
                weighted_resolution(:,k,l)=R;
                                
                if target>0                 % caso contrasto di suscettività positivo
                    cons=find(S<0);         % constraint inferiore
                else                        % caso contrasto di suscettività negativo
                    cons=find(S>0);         % constraint superiore
                end
                S(cons)=0;                  % se il contrasto è positivo azzero le celle con S negativo (e viceversa)
                Fobs_red0=Fobs_red;
                S_iter(:,k,l)=S;     % la matrice S_iter è 3D: i vettori delle suscettività vengono inseriti in colonne, 
                                     % e ogni colonna corrisponde ad un'iterazione;
                                     % la terza dimensione serve perché l'intero
                                     % ciclo viene computato per ogni eps_sus

                Fcal=A*S;          % forward per ogni inversione
                Fcalcolato(:,k,l)=Fcal;     % il forward viene posizionato nella matrice 3D
                residual=(F_obs_int-Fcal);    % residual
                L2norm(l,k)=norm(residual);        % norma2 su tutto il segnale

                siguti=find(mag2data.sig==1);              % calcolo la norma2 riferita alla parte utile del segnale
                resuti=(F_obs_int(siguti)-Fcal(siguti));
                L2normuti(l,k)=norm(resuti);  % norma2 della parte di segnale utile
            end
            close(wait_iter)
        end
        close(wait)
        toc
end

min_k=0;
while min_k<4           % si sceglie S tale che abbia norma2 min, ma si escludono le prime iterazioni
    [val_min,iter_min] = min(L2norm');
    min_l=find(val_min==min(val_min));
    min_k=iter_min(min_l);
        if min_k<4
            L2norm(min_l,min_k)=L2norm(min_l,min_k)+100;
        end
end

min_kuti=0;
while min_kuti<4        % si sceglie Suti tale che abbia norma2 min, ma si escludono le prime iterazioni
        [val_minuti,iter_minuti] = min(L2normuti');
        min_luti=find(val_minuti==min(val_minuti));
        min_kuti=iter_minuti(min_luti);
            if min_kuti<4
                L2normuti(min_luti,min_kuti)=L2normuti(min_luti,min_kuti)+100;
            end
end

resall=weighted_resolution(:,:,min_l);
mag2data.resall=resall;
if min_l~=min_luti
    resuti=weighted_resolution(:,:,min_luti);
    mag2data.resuti=resuti;
    mag2data.controls.differentresolutions=1;
end

% plotto l'andamento della norma2 
iter=[1:Max_iter];
figure('Name','L2norm Trend')
plot(iter,L2norm(min_l,:),'r', iter,L2normuti(min_luti,:),'b', 'Linewidth',1.5), grid on
legend('L_2norm entire signal','L_2norm "useful" signal')
xlabel('Iterations')
ylabel('L_2norm')
set(gca,'XTick',iter,'XLim',[1 Max_iter])


fprintf('\n%%%% Inversion result ("entire signal") %%%%\n');
fprintf('Perturbation number = %i\n',eps_sus(min_l));
fprintf('Iteration = %i\n',min_k);
L2norm=L2norm(min_l,min_k);
fprintf('L2norm = %i\n',L2norm);

mag2data.min_k=min_k;
mag2data.S=S_iter(:,min_k,min_l);
mag2data.S_iter=S_iter(:,:,min_l);
mag2data.Fcalcolato=Fcalcolato(:,min_k,min_l);
mag2data.eps_sus=eps_sus(min_l);
mag2data.min_k=min_k;
mag2data.L2norm=L2norm;


fprintf('\n%%%% Inversion result ("useful signal") %%%%\n');
fprintf('Perturbation number = %i\n',eps_sus(min_luti));
fprintf('Iteration = %i\n',min_kuti);
L2normuti=L2normuti(min_luti,min_kuti);
fprintf('L2norm = %i\n',L2normuti);

mag2data.min_kuti=min_kuti;
mag2data.Suti=S_iter(:,min_kuti,min_luti);
mag2data.S_iteruti=S_iter(:,:,min_luti);
mag2data.Fcalcolatouti=Fcalcolato(:,min_kuti,min_luti);
mag2data.eps_susuti=eps_sus(min_luti);
mag2data.min_kuti=min_kuti;
mag2data.L2normuti=L2normuti;


fprintf('\n%%%% Inversion parameters %%%%\n');
fprintf('Max iteration = %i\n',Max_iter);
fprintf('Noise/Signal = %i\n',l_0);
fprintf('Susceptibility contrast = %i\n\n',target*4*pi);