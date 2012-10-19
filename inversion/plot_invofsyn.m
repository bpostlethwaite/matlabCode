function plot_invofsyn

global mag2data

%%
% plot_invofsyn plots the result of the inversion process of synthetic data.
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


sig=mag2data.sig;
CX=mag2data.CX;
CZ=mag2data.CZ;
S=mag2data.S.*4*pi.*10^3;
Suti=mag2data.Suti.*4*pi.*10^3;
Xi=mag2data.Xi;
xx=mag2data.xx;
Fobs=mag2data.Fobs;
F_obs_int=mag2data.F_obs_int;
Fcalcolato=mag2data.Fcalcolato;
Fcalcolatouti=mag2data.Fcalcolatouti;
x_centrocella=mag2data.x_centrocella;
z_centrocella=mag2data.z_centrocella;
xspace=mag2data.xspace;
zspace=mag2data.zspace;
matrix_celle=mag2data.model.matrix_celle.*10^3;
matrix_celle_inv=reshape(S,CX,CZ);
matrix_celleuti_inv=reshape(Suti,CX,CZ);


xlab=' ';
xlabtop=min(Xi);
for i=1:round(max(Xi)/(xspace(6)-xspace(1)))
    for j=1:4
        xlab=strvcat(xlab,' ');
    end
    xlab=strvcat(xlab, num2str((xspace(6)-xspace(1))*i));
    xlabtop=[xlabtop, (xspace(6)-xspace(1))*i];
end

% xlab=' ';
% for i=1:round(max(Xi)/10)
%     for j=1:10/xspace(2)-1
%         xlab=strvcat(xlab,' ');
%     end
%     xlab=strvcat(xlab,num2str(10*i));
% end

if length(zspace)>5
    zlab=[];
    for i=1:length(zspace)
        if round(i/2)~=i/2
            zlab=strvcat(zlab, num2str(zspace(i)));
        else
            zlab=strvcat(zlab, ' ');
        end
    end
else
    zlab=zspace;
end


switch mag2data.controls.conjugatekernels
    case 0
        if mag2data.controls.campo_tot==1
            Foss='Observed total field';
            Fcalc='Computed total field';
            assey='Mag. anomaly (nT)';
        else
            Foss='Observed vertical gradient';
            Fcalc='Computed vertical gradient';
            assey='Mag. anomaly (nT/m)';
        end

        figure('Name','Inversion Result - Entire Signal')
        subplot(2,1,1)
        plot(xx,Fobs,'.', Xi,Fcalcolato,'r', 'MarkerSize',8, 'LineWidth',1.5);
        legend(Foss,Fcalc)
        xlabel('Distance (m)')
        ylabel(assey)
        set(gca,'XLim',[min(Xi) max(Xi)], 'XTick',xlabtop)
%         subplot(3,1,2)
%         imagesc(x_centrocella(1,:),z_centrocella(:,1),matrix_celle'), shading flat, grid on
%         set(gca, 'XTick',xspace, 'YTick',zspace, 'XTickLabel',[])
%         axis('ij')
%         ylabel('Depth (m)')
        subplot(2,1,2)
        imagesc(x_centrocella(1,:),z_centrocella(:,1),matrix_celle_inv'), shading flat, colorbar('horiz'), grid on
        set(gca, 'XTick',xspace, 'YTick',zspace, 'XAxisLocation','top', 'XTickLabel',xlab, 'YTickLabel',zlab)
        axis('ij')
        ylabel('Depth (m)')
        annotation('textbox',[0.89 0.025 0.1 0.05], 'String','x 10^-^3', 'Linestyle','none');

        figure('Name','Inversion Result - "Useful" Signal')
        subplot(2,1,1)
        plot(xx,Fobs,'.', Xi,Fcalcolatouti,'r', 'MarkerSize',8, 'LineWidth',1.5);
        legend(Foss,Fcalc)
        xlabel('Distance (m)')
        ylabel(assey)
        set(gca,'XLim',[min(Xi) max(Xi)], 'XTick',xlabtop)
        hold on
        rett=0;
        for i=2:length(Xi)
            diffsig=sig(i)-sig(i-1);
            if diffsig>0
                rett=rett+1;
                xirett(rett)=Xi(i);
            end
            if diffsig<0
                xfrett(rett)=Xi(i-1);
            end
        end
        for j=1:rett
            rectangle('Position',[xirett(j),min(F_obs_int),(xfrett(j)-xirett(j)),(max(F_obs_int)-min(F_obs_int))]), hold on
        end
        hold off
    case 1
        Fossbottom='Observed total field bottom sensor';
        Fosstop='Observed total field top sensor';
        Fcalcbottom='Computed total field bottom sensor';
        Fcalctop='Computed total field top sensor';
        assey='Mag. anomaly (nT)';
        figure('Name','Inversion Result - Entire Signal')
        subplot(2,1,1)
        plot(xx,Fobs(1:length(xx)),'.', xx,Fobs(length(xx)+1:end),'.', Xi,Fcalcolato(1:length(Xi)),'r', Xi,Fcalcolato(length(Xi)+1:end),'m', 'MarkerSize',8, 'LineWidth',1.5);
        legend(Fossbottom,Fosstop,Fcalcbottom,Fcalctop)
        xlabel('Distance (m)')
        ylabel(assey)
        set(gca,'XLim',[min(Xi) max(Xi)], 'XTick',xlabtop)
%         subplot(3,1,2)
%         imagesc(x_centrocella(1,:),z_centrocella(:,1),matrix_celle'), shading flat, grid on
%         set(gca, 'XTick',xspace, 'YTick',zspace, 'XTickLabel',[])
%         axis('ij')
%         ylabel('Depth (m)')
        subplot(2,1,2)
        imagesc(x_centrocella(1,:),z_centrocella(:,1),matrix_celle_inv'), shading flat, colorbar('horiz'), grid on
        set(gca, 'XTick',xspace, 'YTick',zspace, 'XAxisLocation','top', 'XTickLabel',xlab, 'YTickLabel',zlab)
        axis('ij')
        ylabel('Depth (m)')

        annotation('textbox',[0.89 0.025 0.1 0.05], 'String','x 10^-^3', 'Linestyle','none', 'FontSize',12);

        figure('Name','Inversion Result - "Useful" Signal')
        subplot(2,1,1)
        plot(xx,Fobs(1:length(xx)),'.', xx,Fobs(length(xx)+1:end),'.', Xi,Fcalcolatouti(1:length(Xi)),'r', Xi,Fcalcolatouti(length(Xi)+1:end),'m', 'MarkerSize',8, 'LineWidth',1.5);
        legend(Fossbottom,Fosstop,Fcalcbottom,Fcalctop)
        xlabel('Distance (m)')
        ylabel(assey)
        set(gca,'XLim',[min(Xi) max(Xi)])
        hold on
        for k=1:2
            rett=0;
            for i=2:length(Xi)
                diffsig=sig(i)-sig(i-1);
                if diffsig>0
                    rett=rett+1;
                    xirett(rett)=Xi(i);
                end
                if diffsig<0
                    xfrett(rett)=Xi(i-1);
                end
            end
            for j=1:rett
                in=(k-1)*(length(Xi))+1;
                fin=(k)*(length(Xi));
                rectangle('Position',[xirett(j),min(F_obs_int(in:fin)),(xfrett(j)-xirett(j)),max(F_obs_int(in:fin))-min(F_obs_int(in:fin))]), hold on
            end
            sig=sig(length(Xi)+1:end);
        end
        hold off
end        
        
% subplot(3,1,2)
% imagesc(x_centrocella(1,:),z_centrocella(:,1),matrix_celle'), shading flat, grid on
% set(gca, 'XTick',xspace, 'YTick',zspace, 'XTickLabel',[])
% axis('ij')
% ylabel('Depth (m)')
subplot(2,1,2)
imagesc(x_centrocella(1,:),z_centrocella(:,1),matrix_celleuti_inv'), shading flat, colorbar('horiz'), grid on
set(gca, 'XTick',xspace, 'YTick',zspace, 'XAxisLocation','top', 'XTickLabel',xlab, 'YTickLabel',zlab)
axis('ij')
ylabel('Depth (m)')

annotation('textbox',[0.89 0.025 0.1 0.05], 'String','x 10^-^3', 'Linestyle','none');