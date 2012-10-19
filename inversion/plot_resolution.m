function plot_resolution

global mag2data

%%
% plot_resolution plots the model resolution matrix and its main diagonal.
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


CX=mag2data.CX;
CZ=mag2data.CZ;
x_centrocella=mag2data.x_centrocella;
z_centrocella=mag2data.z_centrocella;

switch mag2data.controls.pesi
    case 0
        R=mag2data.R;
%         figure        
%         imagesc(mag2data.datares), colorbar('horiz')
%         axis('ij')
%         title('Data Resolution Matrix')
        figure('Name','Model Resolution Matrix')
        imagesc(mag2data.res), colorbar('horiz')
        axis('ij')
        title('Model Resolution Matrix')
        figure('Name','Main Diagonal of the Model Resolution Matrix')
        XX_R=reshape(R,CX,CZ);
        imagesc(x_centrocella(1,:),z_centrocella(:,1),XX_R'), shading flat
        axis('ij')
        title('Main Diagonal of the Model Resolution Matrix')
        xlabel('Distance (m)')
        ylabel('Depth (m)')
        set(gca, 'CLim',[0 1], 'XAxisLocation','top')
        colorbar('horiz')
        figure('Name','Covariance Matrix')
        imagesc(mag2data.covm), colorbar('horiz')
        axis('ij')
        title('Covariance Matrix')
    case 1
        Max_iter=mag2data.Max_iter;
        resall=mag2data.resall;
        min_k=mag2data.min_k;
        switch Max_iter
            case 10
                figure
                for i=1:Max_iter                    
                    XX_R=reshape(resall(:,i),CX,CZ);
                    subplot(5,2,i)
                    set(gca,'FontSize',8);
                    imagesc(x_centrocella(1,:),z_centrocella(:,1),XX_R'), shading flat 
                    axis('ij')
                    if i==min_k
                        title(['Iteration',int2str(i)],'Color','r')
                    else
                        title(['Iteration',int2str(i)])
                    end
                    set(gca, 'CLim',[0 1])
%                     colorbar('EastOutside','FontSize',8)
%                     xlabel('Distance (m)')
%                     ylabel('Depth (m)')
                end
            otherwise
                for i=1:Max_iter
                    figure
                    XX_R=reshape(resall(:,i),CX,CZ);
                    imagesc(x_centrocella(1,:),z_centrocella(:,1),XX_R'), shading flat, colorbar('horiz')
                    axis('ij')
                    if i==min_k
                        title(['Iteration',int2str(i)],'Color','r')
                    else
                        title(['Iteration',int2str(i)])
                    end
                    xlabel('Distance (m)')
                    ylabel('Depth (m)')
                end
        end
        if mag2data.controls.differentresolutions==1;
            resuti=mag2data.resuti;
            min_kuti=mag2data.min_kuti;
            switch Max_iter
                case 10
                    figure
                for i=1:Max_iter                    
                    XX_R=reshape(resuti(:,i),CX,CZ);
                    subplot(5,2,i)
                    set(gca,'FontSize',8);
                    imagesc(x_centrocella(1,:),z_centrocella(:,1),XX_R'), shading flat
                    axis('ij')
                    if i==min_kuti
                        title(['Iteration',int2str(i)],'Color','r')
                    else
                        title(['Iteration',int2str(i)])
                    end
                    set(gca, 'CLim',[0 1])
%                     colorbar('EastOutside','FontSize',8)
%                     xlabel('Distance (m)')
%                     ylabel('Depth (m)')
                end
                otherwise
                    for i=1:Max_iter
                        figure
                        XX_R=reshape(resuti(:,i),CX,CZ);
                        imagesc(x_centrocella(1,:),z_centrocella(:,1),XX_R'), shading flat, colorbar('horiz')
                        axis('ij')
                        if i==min_kuti
                            title(['Iteration',int2str(i)],'Color','r')
                        else
                            title(['Iteration',int2str(i)])
                        end
                        xlabel('Distance (m)')
                        ylabel('Depth (m)')
                    end
            end
        end
end