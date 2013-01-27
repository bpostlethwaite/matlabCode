function kernelFpseudo3D

global mag2data

%%
% kernelFpseudo3D computes the kernel for the case of total field anomaly (21/2D).
% 
% Written by Stefano Stocco on September, 2007

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

 
switch mag2data.controls.computation
    case 'computing_forward'
        %% Input
        N_dati_ori = mag2data.N_dati_ori;   % Number of original data
        lp = mag2data.lp;   % Length of profile
        order = mag2data.order; % Order of interpolation
        hsi = mag2data.hsi; % Height of bottom sensor
        hss = mag2data.hss; % Height of top sensor
        sus = mag2data.model.sus/(4*pi);    % Susceptibility distribution
        angoloBeta = mag2data.controls.angoloBeta;
        sequenzaBeta = mag2data.controls.sequenzaBeta;
        Beta = (90-mag2data.Beta)*pi/180;   % Angle north/profile
        I = mag2data.I*pi/180;  % Inclination of magnetic field
        F=mag2data.F;
        CX = mag2data.CX;   % Prisms along x
        CZ = mag2data.CZ;   % Prisms along z
        dz = mag2data.dz;   % thickness of prisms
        L = mag2data.L/2;
        alreadyforward = mag2data.controls.alreadyforward;
        % angoloBeta shows if, in a trial-and-error, Beta changes, so I don't
        % compute the kernel if it is not necessary
        angoloBeta(sequenzaBeta)=Beta;
        diffBeta=0;
        if sequenzaBeta>1
            diffBeta=angoloBeta(sequenzaBeta)-angoloBeta(sequenzaBeta-1);
        end
        Fobs=ones(N_dati_ori,1);
        steplp=lp/(N_dati_ori-1);
        xx=[0:steplp:lp];
        % Interpolation of experimental data
        if order~=1
            dstat=(xx(2)-xx(1))/order;
            x=[min(xx):dstat:max(xx)];
            F_obs_int=interp1(xx,Fobs,x,'linear','extrap');
        else
            x=xx;
            F_obs_int=Fobs;
        end
        N_dati=length(F_obs_int);
        h=[hsi,hss];
        dx=lp/CX;   % Prisms dimension along x
        CTOT=CX*CZ; % Total number of prisms

        if (alreadyforward==1 & diffBeta==0)

        else
            disp('Wait... kernel computing');
            tic
            for k=1:2;

                % Kernel definition

                Xi=x-min(x);
                Zi=[0:dz:CZ];

                xj=[dx/2:dx:CX*dx-dx/2];    % centre (x) of prisms
                zj=[dz/2:dz:CZ*dz-dz/2];    % centre (z) of prisms
                [XJ,ZJ]=meshgrid(xj,zj);
                Xj=reshape(XJ',1,size(XJ,1)*size(XJ,2));
                Zj=reshape(ZJ',1,size(ZJ,1)*size(ZJ,2));

                Zj=Zj+h(k);

                % Distance and angles
                A=ones(N_dati,CTOT);
                Ximat=diag(Xi)*A;
                Xjmat=A*diag(Xj);
                Zjmat=A*diag(Zj);
                shiftx=(dx/2)*A;
                shiftz=(dz/2)*A;

                Xplus=Ximat-Xjmat+shiftx+eps;
                Xminus=Ximat-Xjmat-shiftx+eps;
                Zplus=Zjmat+shiftz;
                Zminus=Zjmat-shiftz;
                Xplus2=Xplus.^2;
                Xminus2=Xminus.^2;
                Zplus2=Zplus.^2;
                Zminus2=Zminus.^2;

                % Distance and angles
                r1(:,:)=(Zminus2+Xplus2).^0.5;
                r2(:,:)=(Zplus2+Xplus2).^0.5;
                r3(:,:)=(Zminus2+Xminus2).^0.5;
                r4(:,:)=(Zplus2+Xminus2).^0.5;
                teta1(:,:)=atan(Zminus./Xplus);
                teta2(:,:)=atan(Zplus./Xplus);
                teta3(:,:)=atan(Zminus./Xminus);
                teta4(:,:)=atan(Zplus./Xminus);
                delta_X(:,:)=Ximat-Xjmat+2*eps;
                delta_Z(:,:)=Zjmat-shiftz;

                % Synthetic anomaly of prisms
                at=log((r1.^2+L^2).^0.5+L)-log((r1.^2+L^2).^0.5-L)+log((r3.^2+L^2).^0.5-L)-log((r3.^2+L^2).^0.5+L);
                bt=atan(L./delta_X)-atan(L./(delta_X-dx))-atan(L.*delta_Z./(delta_X.*((r1.^2+L^2).^0.5)))+atan(L.*delta_Z./((delta_X-dx).*((r3.^2+L^2).^0.5)));
                ct=cos(I)^2*sin(Beta)-sin(I)^2;
                Asup=F*(sin(2*I)*sin(Beta).*at-ct.*bt);          % matrice geometrica riferita alla parte superiore del prisma

                dt=log((r2.^2+L^2).^0.5+L)-log((r2.^2+L^2).^0.5-L)+log((r4.^2+L^2).^0.5-L)-log((r4.^2+L^2).^0.5+L);
                et=atan(L./delta_X)-atan(L./(delta_X-dx))-atan(L.*(delta_Z+dz)./(delta_X.*((r2.^2+L^2).^0.5)))+atan(L.*(delta_Z+dz)./((delta_X-dx).*((r4.^2+L^2).^0.5)));
                Ainf=F*(sin(2*I)*sin(Beta).*dt-ct.*et);          % matrice geometrica riferita alla parte inferiore del prisma

                if k==1
                    mag2data.Ai=Asup-Ainf;          % matrice geometrica (kernel) bottom sensor
                else    % k==2
                    mag2data.As=Asup-Ainf;          % matrice geometrica (kernel) top sensor
                end

            end
        end     % è la fine dell'if alreadyforward
        toc

        mag2data.Ai_s=mag2data.Ai-mag2data.As;

        if alreadyforward==1
        else
            mag2data.Xi=Xi;
            mag2data.controls.alreadyforward=1;
        end

        mag2data.xx=xx;
        mag2data.CTOT=CTOT;
        mag2data.controls.angoloBeta = angoloBeta;
        mag2data.controls.sequenzaBeta=sequenzaBeta+1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'computing_inversion'
        %% Input
        Fobs = mag2data.Fobs;   % Magnetic data
        xx = mag2data.xx;   % Coordinate of profile
        order = mag2data.order; % Order of interpolation
        hsi = mag2data.hsi; % Heigth of bottom sensor
        hss = mag2data.hss; % Heigth of top sensor
        F = mag2data.F; % Earth magnetic field intensity
        I = mag2data.I*pi/180;  % Inclination of magnetic field
        Beta = (90-mag2data.Beta)*pi/180;   % Angle north/profile
        CX = mag2data.CX;   % Prisms along x
        CZ = mag2data.CZ;   % Prisms along z
        dz = mag2data.dz;   % thickness of prisms
        L = mag2data.L/2;

        if mag2data.controls.Risoluzione~=1  % if I have already computed the kernel for the model resolution, I do not compute it again
            N_dati_ori=length(Fobs);
            lp=max(xx)-min(xx);
            % Interpolation of experimental data
            if order~=1
                dstat=(xx(2)-xx(1))/order;
                x=[min(xx):dstat:max(xx)];
                F_obs_int=interp1(xx,Fobs,x,'linear','extrap');
                F_obs_int=F_obs_int';
            else
                dstat=xx(2)-xx(1);
                x=xx;
                F_obs_int=Fobs;
            end
            N_dati=length(F_obs_int);
            h=[hsi,hss];
            dx=lp/CX;   % Prisms dimension along x
            CTOT=CX*CZ; % Total number of prisms
            
            disp('Wait... kernel computing');
            tic
            for k=1:2
                %%
                % Kernel definition

                Xi=x-min(x);
                Zi=[0:dz:CZ];

                xj=[dx/2:dx:CX*dx-dx/2];    % centre (x) of prisms
                zj=[dz/2:dz:CZ*dz-dz/2];    % centre (z) of prisms
                [XJ,ZJ]=meshgrid(xj,zj);
                Xj=reshape(XJ',1,size(XJ,1)*size(XJ,2));
                Zj=reshape(ZJ',1,size(ZJ,1)*size(ZJ,2));

                Zj=Zj+h(k);

                % Distance and angles
                A=ones(N_dati,CTOT);
                Ximat=diag(Xi)*A;
                Xjmat=A*diag(Xj);
                Zjmat=A*diag(Zj);
                shiftx=(dx/2)*A;
                shiftz=(dz/2)*A;

                Xplus=Ximat-Xjmat+shiftx+eps;
                Xminus=Ximat-Xjmat-shiftx+eps;
                Zplus=Zjmat+shiftz;
                Zminus=Zjmat-shiftz;
                Xplus2=Xplus.^2;
                Xminus2=Xminus.^2;
                Zplus2=Zplus.^2;
                Zminus2=Zminus.^2;

                % Distance and angles
                r1(:,:)=(Zminus2+Xplus2).^0.5;
                r2(:,:)=(Zplus2+Xplus2).^0.5;
                r3(:,:)=(Zminus2+Xminus2).^0.5;
                r4(:,:)=(Zplus2+Xminus2).^0.5;
                teta1(:,:)=atan(Zminus./Xplus);
                teta2(:,:)=atan(Zplus./Xplus);
                teta3(:,:)=atan(Zminus./Xminus);
                teta4(:,:)=atan(Zplus./Xminus);
                delta_X(:,:)=Ximat-Xjmat+2*eps;
                delta_Z(:,:)=Zjmat-shiftz;

                % Synthetic anomaly of prisms

                at=log((r1.^2+L^2).^0.5+L)-log((r1.^2+L^2).^0.5-L)+log((r3.^2+L^2).^0.5-L)-log((r3.^2+L^2).^0.5+L);
                bt=atan(L./delta_X)-atan(L./(delta_X-dx))-atan(L.*delta_Z./(delta_X.*((r1.^2+L^2).^0.5)))+atan(L.*delta_Z./((delta_X-dx).*((r3.^2+L^2).^0.5)));
                ct=cos(I)^2*sin(Beta)-sin(I)^2;
                Asup=F*(sin(2*I)*sin(Beta).*at-ct.*bt);          % matrice geometrica riferita alla parte superiore del prisma

                dt=log((r2.^2+L^2).^0.5+L)-log((r2.^2+L^2).^0.5-L)+log((r4.^2+L^2).^0.5-L)-log((r4.^2+L^2).^0.5+L);
                et=atan(L./delta_X)-atan(L./(delta_X-dx))-atan(L.*(delta_Z+dz)./(delta_X.*((r2.^2+L^2).^0.5)))+atan(L.*(delta_Z+dz)./((delta_X-dx).*((r4.^2+L^2).^0.5)));
                Ainf=F*(sin(2*I)*sin(Beta).*dt-ct.*et);          % matrice geometrica riferita alla parte inferiore del prisma

                if k==1
                    Ai=Asup-Ainf;          % kernel bottom sensor

                    % depth weighting matrix
                    vector=[1:CX:CTOT];
                    depth_factor=abs(Ai(1,vector)).^(-1/2);
                    for v=1:length(vector)
                        Wz(1,CX*(v-1)+1:CX*v)=depth_factor(v);
                        mag2data.Wz=diag(Wz,0);
                    end
                    mag2data.Ai=Ai;

                else    % k==2
                    mag2data.As=Asup-Ainf;          % kernel top sensor
                end

            end
            toc
            
            mag2data.Ai_s=mag2data.Ai-mag2data.As;    % kernel vertical gradient

            mag2data.Xi=Xi;
            mag2data.dstat=dstat;
            mag2data.N_dati=N_dati;
            mag2data.F_obs_int=F_obs_int;
            mag2data.CTOT=CTOT;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'computing_synthetic'
        %% Input
        N_dati_ori = mag2data.N_dati_ori;   % Magnetic data
        lp = mag2data.lp;   % Coordinate of profile
        order = mag2data.order; % Order of interpolation
        hsi = mag2data.hsi; % Heigth of bottom sensor
        hss = mag2data.hss; % Heigth of top sensor
        F = mag2data.F; % Earth magnetic field intensity
        I = mag2data.I*pi/180;  % Inclination of magnetic field
        Beta = (90-mag2data.Beta)*pi/180;   % Angle north/profile
        CX = mag2data.CX;   % Prisms along x
        CZ = mag2data.CZ;   % Prisms along z
        dz = mag2data.dz;   % thickness of prisms
        Qn = mag2data.Qn/100;   % Kroinesberger ratio
        I_R = mag2data.I_R*pi/180;  % Remanent I
        Beta_R = (90-mag2data.Beta_R)*pi/180;   % Remanent Beta
        L = mag2data.L/2;

        if mag2data.controls.Risoluzione~=1  % if I have already computed the kernel for the model resolution, I do not compute it again
            Fobs=ones(N_dati_ori,1);
            steplp=mag2data.lp/(N_dati_ori-1);
            xx=[0:steplp:lp];
            % Interpolation of experimental data
            if order~=1
                dstat=(xx(2)-xx(1))/order;
                x=[min(xx):dstat:max(xx)];
                F_obs_int= interp1(xx,Fobs,x,'linear','extrap');
            else
                dstat=xx(2)-xx(1);
                x=xx;
                F_obs_int=Fobs;
            end
            N_dati=length(F_obs_int);
            h=[hsi,hss];        % Height of sensors (bottom,top)
            dx=lp/CX;   % Prisms dimension along x
            CTOT=CX*CZ; % Total number of prisms
            
            disp('Wait... kernel computing');
            tic
            for k=1:2;
                %%
                % Kernel definition

                Xi=x-min(x);
                Zi=[0:dz:CZ];

                xj=[dx/2:dx:CX*dx-dx/2];    % centre (x) of prisms
                zj=[dz/2:dz:CZ*dz-dz/2];    % centre (z) of prisms
                [XJ,ZJ]=meshgrid(xj,zj);
                Xj=reshape(XJ',1,size(XJ,1)*size(XJ,2));
                Zj=reshape(ZJ',1,size(ZJ,1)*size(ZJ,2));

                Zj=Zj+h(k);

                % Distance and angles
                A=ones(N_dati,CTOT);
                Ximat=diag(Xi)*A;
                Xjmat=A*diag(Xj);
                Zjmat=A*diag(Zj);
                shiftx=(dx/2)*A;
                shiftz=(dz/2)*A;

                Xplus=Ximat-Xjmat+shiftx+eps;
                Xminus=Ximat-Xjmat-shiftx+eps;
                Zplus=Zjmat+shiftz;
                Zminus=Zjmat-shiftz;
                Xplus2=Xplus.^2;
                Xminus2=Xminus.^2;
                Zplus2=Zplus.^2;
                Zminus2=Zminus.^2;

                % Distance and angles
                r1(:,:)=(Zminus2+Xplus2).^0.5;
                r2(:,:)=(Zplus2+Xplus2).^0.5;
                r3(:,:)=(Zminus2+Xminus2).^0.5;
                r4(:,:)=(Zplus2+Xminus2).^0.5;
                teta1(:,:)=atan(Zminus./Xplus);
                teta2(:,:)=atan(Zplus./Xplus);
                teta3(:,:)=atan(Zminus./Xminus);
                teta4(:,:)=atan(Zplus./Xminus);
                delta_X(:,:)=Ximat-Xjmat+3*eps;
                delta_Z(:,:)=Zjmat-shiftz;

                % Synthetic anomaly of prisms
                at=log((r1.^2+L^2).^0.5+L)-log((r1.^2+L^2).^0.5-L)+log((r3.^2+L^2).^0.5-L)-log((r3.^2+L^2).^0.5+L);
                bt=atan(L./delta_X)-atan(L./(delta_X-dx))-atan(L.*delta_Z./(delta_X.*((r1.^2+L^2).^0.5)))+atan(L.*delta_Z./((delta_X-dx).*((r3.^2+L^2).^0.5)));
                ct_induced=cos(I)^2*sin(Beta)-sin(I)^2;
                ct_remanent=cos(I_R)^2*sin(Beta_R)-sin(I_R)^2;
                Asup_induced=F*(sin(2*I)*sin(Beta).*at-ct_induced.*bt);          % matrice geometrica riferita alla parte superiore del prisma
                Asup_remanent=Qn*F*(sin(2*I_R)*sin(Beta_R).*at-ct_remanent.*bt);
                Asup=Asup_induced+Asup_remanent;

                dt=log((r2.^2+L^2).^0.5+L)-log((r2.^2+L^2).^0.5-L)+log((r4.^2+L^2).^0.5-L)-log((r4.^2+L^2).^0.5+L);
                et=atan(L./delta_X)-atan(L./(delta_X-dx))-atan(L.*(delta_Z+dz)./(delta_X.*((r2.^2+L^2).^0.5)))+atan(L.*(delta_Z+dz)./((delta_X-dx).*((r4.^2+L^2).^0.5)));
                Ainf_induced=F*(sin(2*I)*sin(Beta).*dt-ct_induced.*et);          % matrice geometrica riferita alla parte inferiore del prisma
                Ainf_remanent=Qn*F*(sin(2*I_R)*sin(Beta_R).*dt-ct_remanent.*et);
                Ainf=Ainf_induced+Ainf_remanent;

                if k==1
                    Ai=Asup-Ainf;          % kernel bottom sensor

                    % depth weighting matrix
                    vector=[1:CX:CTOT];
                    depth_factor=abs(Ai(1,vector)).^(-1/2);
                    for v=1:length(vector)
                        Wz(1,CX*(v-1)+1:CX*v)=depth_factor(v);
                        mag2data.Wz=diag(Wz,0);
                    end
                    mag2data.Ai=Ai;

                else    % k==2
                    As=Asup-Ainf;          % kernel top sensor
                    mag2data.As=As;
                end

            end
            toc
            
            mag2data.Ai_s=Ai-As;

            mag2data.A=[Ai;As];
            mag2data.xx=xx;
            mag2data.dstat=dstat;
            mag2data.Xi=Xi;
            mag2data.N_dati=N_dati;
            mag2data.CTOT=CTOT;
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%