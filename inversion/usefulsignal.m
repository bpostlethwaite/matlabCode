function usefulsignal

global mag2data

%%
% usefulsignal searchs the useful part of the signal (anomaly).
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


switch mag2data.controls.conjugatekernels
    case 0
        % Dati in: Fobs, lp, dstat=step
        Fobs=mag2data.F_obs_int;
        lp=mag2data.lp;
        step=mag2data.dstat;
        Xi=mag2data.Xi;
        % parametri
        if round(length(Fobs)/2)==length(Fobs)/2
            zero_ped=length(Fobs)+1;
        else
            zero_ped=length(Fobs);
        end
        %zero_ped=512;
        kny=pi/step;
        deltak=2*pi/((zero_ped-1)*step);
        k=0:deltak:kny;   % vettore dei numeri d'onda
        finestra=blackmanharris(length(Fobs));
        Fobs=Fobs.*finestra;
        % trasformata di Fourier
        A=fft(Fobs,zero_ped);
        A=abs(A(1:length(k)));
        AdB=-20*log10(max(A)./A);
        
%         % spettro Fobs
%         figure
%         plot(k,AdB), grid on
%         xlabel('Wave number k')
%         ylabel('dB')
        
        % imposto la soglia e trovo lambda corrispondente
        soglia1=(-12);
        m=find(AdB==max(AdB));
        while AdB(m)>soglia1
            m=m+1;
        end
        win=((soglia1-AdB(m))*(k(m-1)-k(m))/(AdB(m-1)-AdB(m)))+k(m);
        win=2*pi/win;
        window=win/step;
        window=round(window);
        fprintf('\nRunning window length (m)= %.2f\n',win);
        % valore quadratico medio lungo la finestra
        evenodd=window/2;
        if round(evenodd)==evenodd
            pos=(window)/2;
        else 
            pos=(window+1)/2;
        end
        fine=length(Fobs)-window+1;
        varamp=zeros(length(Xi),1);
        F_quad=Fobs.^2;
        for i=1:fine
            varamp(pos)=(sum(F_quad(i:i+window-1)))/(window-1);
            pos=pos+1;
        end
        null=find(varamp==0);
        varamp(null)=NaN;
        ampdB=-10*log10(max(varamp)./varamp);
        % plotto la potenza del segnale lungo x e mostro le soglie
        soglie=[-84 -72 -60 -54 -48 -42 -36 -30 -24 -18 -12 -6];
        figure('Name','Select a Threshold for the "Mask" Signal','Tag','Potenza segnale');
        plot(Xi,ampdB), grid on, hold on
        line([0 max(Xi)],[-6 -6],'Color','g'), hold on
        line([0 max(Xi)],[-12 -12],'Color','r'), hold on
        line([0 max(Xi)],[-18 -18],'Color','black'), hold on
        line([0 max(Xi)],[-24 -24],'Color','m'), hold on
        line([0 max(Xi)],[-30 -30],'Color','g'), hold on
        line([0 max(Xi)],[-36 -36],'Color','r'), hold on
        line([0 max(Xi)],[-42 -42],'Color','black'), hold on
        line([0 max(Xi)],[-48 -48],'Color','m'), hold on
        line([0 max(Xi)],[-54 -54],'Color','g'), hold on
        line([0 max(Xi)],[-60 -60],'Color','c')
        title('Running Signal Power \gamma^2')
        xlabel('Distance (m)')
        ylabel('dB')
        set(gca,'XLim',[0 max(Xi)])
        set(gca,'YTick',soglie)
        pot_seg=findobj('Tag','Potenza segnale');
        % scelgo la soglia appropriata
        [asc,ord]=ginput(1);
        for lin=1:length(soglie)
            dst(lin)=abs(ord-soglie(lin));
        end
        mindst=find(dst==min(dst));
        soglia2=soglie(mindst);
        sup=find(ampdB>soglia2);
        sig=zeros(length(Xi),1);
        sig(sup)=1;
        
%         cacca=1;
%         if exist('cacca')
%         figure
%         subplot(2,1,1)
%         plot(Xi,Fobs), grid on
%         xlabel('Distanza (m)')
%         ylabel('Magnetic anomaly (nT)')
%         set(gca,'XLim',[0 max(Xi)])
%         
%         subplot(2,1,2)
%         stairs(Xi,sig), grid on
%         set(gca,'YLim',[0 8])
%         set(gca,'XLim',[0 max(Xi)])
%         end

        mag2data.sig=sig;
        close (pot_seg)
    case 1
        % Dati in: Fobs, lp, dstat=step
        lp=mag2data.lp;
        step=mag2data.dstat;
        Xi=mag2data.Xi;
        for i=1:2
            if i==1
                Fobs=mag2data.F_inf;
            else
                Fobs=mag2data.F_sup;
            end
            % parametri
            if round(length(Fobs)/2)==length(Fobs)/2
                zero_ped=length(Fobs)+1;
            else
                zero_ped=length(Fobs);
            end
            %zero_ped=512;
            kny=pi/step;
            deltak=2*pi/((zero_ped-1)*step);
            k=0:deltak:kny;   % vettore dei numeri d'onda
            finestra=blackmanharris(length(Fobs));
            Fobs=Fobs.*finestra;
            % trasformata di Fourier
            A=fft(Fobs,zero_ped);
            A=abs(A(1:length(k)));
            AdB=-20*log10(max(A)./A);
            % % spettro Fobs
            % figure
            % plot(k,AdB), grid on
            % xlabel('Wave number k')
            % ylabel('dB')
            % imposto la soglia e trovo lambda corrispondente
            soglia1=(-12);
            m=find(AdB==max(AdB));
            while AdB(m)>soglia1
                m=m+1;
            end
            win=((soglia1-AdB(m))*(k(m-1)-k(m))/(AdB(m-1)-AdB(m)))+k(m);
            win=2*pi/win;
            window=win/step;
            window=round(window);
            fprintf('\nRunning window length (m)= %.2f\n',win);
            % valore quadratico medio lungo la finestra
            evenodd=window/2;
            if round(evenodd)==evenodd
                pos=(window)/2;
            else 
                pos=(window+1)/2;
            end
            fine=length(Fobs)-window+1;
            varamp=zeros(length(Xi),1);
            F_quad=Fobs.^2;
            for j=1:fine
                varamp(pos)=(sum(F_quad(j:j+window-1)))/(window-1);
                pos=pos+1;
            end
            null=find(varamp==0);
            varamp(null)=NaN;
            ampdB=-10*log10(max(varamp)./varamp);
            % plotto la potenza del segnale lungo x e mostro le soglie
            soglie=[-84 -72 -60 -54 -48 -42 -36 -30 -24 -18 -12 -6];
            figure('Tag','Potenza segnale');
            plot(Xi,ampdB), grid on, hold on
            line([0 max(Xi)],[-6 -6],'Color','g'), hold on
            line([0 max(Xi)],[-12 -12],'Color','r'), hold on
            line([0 max(Xi)],[-18 -18],'Color','black'), hold on
            line([0 max(Xi)],[-24 -24],'Color','m'), hold on
            line([0 max(Xi)],[-30 -30],'Color','g'), hold on
            line([0 max(Xi)],[-36 -36],'Color','r'), hold on
            line([0 max(Xi)],[-42 -42],'Color','black'), hold on
            line([0 max(Xi)],[-48 -48],'Color','m'), hold on
            line([0 max(Xi)],[-54 -54],'Color','g'), hold on
            line([0 max(Xi)],[-60 -60],'Color','c')
            xlabel('Distance (m)')
            ylabel('dB')
            set(gca,'XLim',[0 max(Xi)])
            set(gca,'YTick',soglie)
            pot_seg=findobj('Tag','Potenza segnale');
            % scelgo la soglia appropriata
            [asc,ord]=ginput(1);
            for lin=1:length(soglie)
                dst(lin)=abs(ord-soglie(lin));
            end
            mindst=find(dst==min(dst));
            soglia2=soglie(mindst);
            sup=find(ampdB>soglia2);
            if i==1
                sig=zeros(2*length(Xi),1);
                sig(sup)=1;
            else
                sig(sup+length(Xi))=1;
            end
            % cacca=1;
            % if exist('cacca')
            % figure
            % subplot(2,1,1)
            % plot(Xi,Fobs), grid on
            % xlabel('Distanza (m)')
            % ylabel('Magnetic anomaly (nT)')
            % set(gca,'XLim',[0 max(Xi)])
            % 
            % subplot(2,1,2)
            % stairs(Xi,sig), grid on
            % set(gca,'YLim',[0 10])
            % set(gca,'XLim',[0 max(Xi)])
            % end
            close (pot_seg)
        end
        mag2data.sig=sig;
end