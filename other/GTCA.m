% Still to do: More fail safes for incorrect names and numbers
% Create a profile parse option, to cut out max amp signal


close all
clear all
%% STEP 1
% Decide whether to compile text data into gridded matlab data or not. If
% data is already compiled skip to step 2
h=msgbox(['Welcome to the gravity & topography compiler and analyzer.'...
    '  This program works with data from the Scripps Institution of Oceanography.'...
'  To download txt data go the website http://topex.ucsd.edu/cgi-bin/get_data.cgi'...
' and enter in the longitude and latitude coordinates (google earth helps with this)'...
' and save the webpage as a text file for use with this program. Please keep the selection'...
' to within a 10 degrees lon & lat for the sake of computing time. If formatted data already exists from previous runs'...
' in gridded .mat format the program gives you the option of skipping the compiling step.'...
'  After the compiling step the program will offer to finish, meaning this'...
' program can be used as a txt to gridded .mat converter. You will be asked to input the names of the two'...
' text files with the gravity and topography data from the Scripps Website'],'GRAVITY AND TOPOGRAPHY COMPILER AND ANALYZER','modal');
waitfor(h)

% Ask user to compile txt data, make sure that the data has the right
% format for names
comp{1}(1:20)=1;
name = [];
while comp{1}(1) == 1  %Loop to query user to either format txt data or continue
prompt{1} = {['Would you like to compile txt data into gridded .mat data? [Y/N]'...
    ' (enter ''N'' if this is not your first run and want to use existing gridded .mat format data):']};
dlg_title{1} = 'Compile txt data into usable gridded .mat file';
def = {'Y'};
num_lines{1} = 1;
P(1) = inputdlg(prompt{1},dlg_title{1},num_lines{1},def);
Q(1) = cell2mat(P(1));
    if isempty(Q(1))
        Q(1) = 'Y'; end % Set defaults
    if Q(1) == 'y' || Q(1) == 'Y'
            comp{1}(1)=0; % Move to next loop
    elseif Q(1) == 'n' || Q(1) == 'N';
        comp{1}(2)=0; % Pass key for NO to gridding txt data, skip next loop
            comp{1}(1)=0; % Move to next loop
    else
        warndlg('Please enter either the character ''Y'' or ''N'' in the selection box'); 
    end       
end

while comp{1}(2) == 1 %Loop for YES to gridding data, first query directory to see if files are
    %   in proper format
    prompt{2} = {'Enter the name of the text file holding the topography data (do not add .txt extension)', ['Enter the name'...
        ' of the text file holding the gravity data (do not add .txt extension)']};
    dlg_title{2} = 'Text Data from Scripps Website';
    name5=inputdlg(prompt{2},dlg_title{2},num_lines{1}); name3 = cell2mat(name5(1)); name4 = cell2mat(name5(2));
    
    Q(2) = exist(sprintf('%s.txt',name3));
    Q(3) = exist(sprintf('%s.txt',name4));
  if Q(2) == 2 && Q(3) == 2 %Prompt for name then run compiler
    prompt{2} = {['Enter the name which you would like to call the .mat file (geographic identifier),'...
        ' do not include the extension. This may take a minute']};
    dlg_title{2} = 'Name .mat File';
    name=inputdlg(prompt{2},dlg_title{2},num_lines{1}); name = cell2mat(name);
    compilR(name,name3,name4);
    comp{1}(2) = 0; % Move on to the next step        
  else
      h=warndlg(sprintf(['Matlab did not find the files ''%s.txt'' & ''%s.txt'' in the search path. '...
          'Please make sure the text files are in the search path and spelled correctly.'],name3,name4),'!! Warning !!');
      waitfor(h)
end
end

while comp{1}(3) == 1 % This loop will check to see if the user would like to work on the data
    % created with the last step, if it was completed, else ask if the user
    % has other .mat data to offer.
a=1;   
    if isempty(name) == 0 && a == 1
        button = questdlg(sprintf('would you like to work on the %s.mat file?',name),...
            'User Query','Yes','No, I want to work on a different file','Just wanted to convert the txt file, I''m all done','Yes');
        if strcmp(button,'Yes')
            comp{1}(3) = 0;
            a = 0;
            break
        elseif strcmp(button,'Just wanted to convert the txt files, I''m all done')
                h=msgbox('Exiting Program, Have a super day');
                waitfor(h)
                return
        else
            break
        end      
     else 
        prompt{2} = {'Please enter the name of the .mat file you would like to continue to work on'};
        dlg_title{2} = 'Select .mat data file';
        name=inputdlg(prompt{2},dlg_title{2},num_lines{1}); name = cell2mat(name);
            comp{1}(3)= 0;
            
     end
end    
%% STEP 2
% Load .mat file and prep (Resize into (m,m)) file for further analysis and plotting
load(name)
comp{1}(4)=0; %Bypass for save .mat file loop if .mat file already resized
h(2)=0;
while length(lon) ~= length(lat)
    comp{1}(4)=1; % Switch for Save loop once resize complete
    h(2)=figure(1);
    axis([min(lon) max(lon) min(lat) max(lat)])
    imagesc(lon,lat,topo)
    set(gca,'ydir','normal')
    xlabel('Longitude')
    ylabel('Latitude')
    title(sprintf('%s Topography (meters), Pending trimming',name));
    shading interp
    view(0,90)
    colorbar
diffL = abs(length(lon)-length(lat));  % Determine length difference between two spatial axis  
    if length(lon) > length(lat)
        h(1)=msgbox('Longitude is longer than latitude, we need to resize this, press OK to continue');
        waitfor(h(1))
        button = questdlg('Which side would you like to trim down?, See figure',...
            'User Query','West Side','East Side',...
            'I have no idea what''s going on, abort','West Side');
        if strcmp(button,'West Side')
            topo = topo(:,diffL+1:length(lon));
            grav = grav(:,diffL+1:length(lon));
            lon=lon(diffL+1:length(lon));
            break
        elseif strcmp(button,'East Side')
            topo = topo(:,1:length(lon)-diffL);
            grav = grav(:,1:length(lon)-diffL);
            lon=lon(1:length(lon)-diffL);
            break
        else
            h(1)=msgbox('Exiting Program, Have a super day');
                waitfor(h(1))
            return
        end
    elseif length(lat) > length(lon)
        h(1)=msgbox('Latitude is longer than longitude, we need to resize this, press OK to continue');
        waitfor(h(1))
        button = questdlg('Which side would you like to trim down?, See figure',...
            'User Query','North Side','South Side',...
            'I have no idea what''s going on, abort','North Side');
        if strcmp(button,'North Side')
            topo=topo(diffL+1:length(lat),:);
            grav=grav(diffL+1:length(lat),:);
            lat=lat(diffL+1:length(lat));
            break
        elseif strcmp(button,'South Side')
            topo=topo(1:length(lat)-diffL,:);
            grav=grav(1:length(lat)-diffL,:);
            lat=lat(1:length(lat)-diffL);
            break
        else
               h(1)=msgbox('Exiting Program, Have a super day');
                waitfor(h(1))
            return
        end
    end
        

end
        if comp{1}(4) == 1 % This provides an opportunity to save the resized data
            
        button = questdlg('Would you like to save the resized .mat data in order to skip the resizing step in the future?',...
            'Save Data','Yes, Save it and continue','No, Don''t save it',...
            'Save and  quit','Yes, Save it and continue');
         if strcmp(button,'Yes, Save it and continue')
            prompt{2} = {'Enter the name which you would like to call the resized .mat file, do not include the extension'};
            dlg_title{2} = 'Name resized .mat File';
            def = {name};
            name2=inputdlg(prompt{2},dlg_title{2},num_lines{1},def); name2 = cell2mat(name2);
            name = name2;
            save(sprintf('%s.mat',name2),'lon','lat','topo','grav')
     
        elseif strcmp(button,'No, Don''t save it')
            name2=name;
            else
               h(1)=msgbox('Exiting Program, Have a super day');
                waitfor(h(1))
            return
         end

        end
%% STEP 3
% Set up the filters, hanning window, and variables, mean depth etc.
% For now it will be automatic, perhaps add user defined variables later.
if h(2)
close 1
end

% Set variables
G = 6.673e-11;
rho_w = 1025;
rho_c = 2800;
rho_m = 3330;
delta_rho = rho_c-rho_w;
E = 6.5e10; v=0.25; g=9.82;
mean_sea = mean(topo(topo<0));
s = abs(mean_sea);
d=6000;
% Create 2D cosine window
width=length(topo(:,1));
hanV = 0.5*(1-cos(2*pi*[0:width-1]'/(0.5*width))); hanV(round(width/4):round(3*width/4),1)=1;
hanH = 0.5*(1-cos(2*pi*[0:width-1]/(0.5*width))); hanH(1,round(width/4):round(3*width/4))=1;
HAN = hanV*hanH;
gravw = grav.*HAN; 

%% Step 4 Select n Profiles and calculate the wave vectors
comp{1}(1:4)=1; % Reset flags


        
while  comp{1}(1) == 1; %Comp(1) flag will be used to allow user to continue
    % repeating the profile process until at some stage user says finish and
    % continue
    
    a =1;
if exist('lonX') == 1 && exist('latY') == 1 && a == 1;
    button = questdlg('Would you like to use your existing Profiles?',...
            'User Query','Yes, use last profiles','No, select new profiles',...
            'abort','No, select new profiles');
            if strcmp(button,'Yes, use last profiles')
                comp{1}(8)=0;
                a=0;
                
            elseif strcmp(button,'No, select new profiles')
                a=0;
                
            else
                return
            end
end

    if comp{1}(8) == 1       
 h(1) = msgbox(['To continue with the analysis we need to pick several profiles'...
     ' from the gravity map. Once this dialogue box is closed you will be asked'... 
     ' to specify how many longitudinal and latitudinal profiles you would like.'...
     ' A figure will be opened to assist in your decision.'...
     ' Please keep the total profiles to 6 and under. If you select profiles over'...
     ' a gravity high (load) and moat system (flexure form) you will be presented with'...
     ' an option to use the load & moat profile system. This is recommended.'],'Dialogue Box','modal');
     waitfor(h(1))
    figure(1); % Open figure to assist in # of longitude and latitude profiles
    axis([min(lon) max(lon) min(lat) max(lat)])
    imagesc(lon,lat,grav)
    set(gca,'ydir','normal')
    xlabel('Longitude')
    ylabel('Latitude')
    title(sprintf('%s Gravity [mgal], Selecting Profiles',name));
    shading interp
    view(0,90)
    colorbar 
    end
            % INPUT # of profiles. Need to check that they are numbers,
            % and they sum to 6 or less. Use flag 2 to achieve this
            while  comp{1}(2) == 1 && comp{1}(8) == 1
            prompt{2} = {'How many Longitude Profiles (North/South):','How many Latitude Profiles (East/West):'};
            dlg_title{2} = 'Pick # of profiles';
            def = {'2','2'};
            nprofile=inputdlg(prompt{2},dlg_title{2},num_lines{1},def); 
            nprofile = [str2num(nprofile{1}), str2num(nprofile{2})];
                if sum(nprofile) <= 6 %Check for max n profiles
                    comp{1}(2) = 0;
                
                else
                    h(1)=msgbox('Please enter 6 or less total profiles');
                    waitfor(h(1))
                end
            end

    while comp{1}(3) == 1 %This flag loop will break after BOTH long and lat
        % profiles have been chosen and verified by the user.
        
        if nprofile(1) > 0 && comp{1}(8) == 1% Select n longitudinal profiles
        h(1) = msgbox(sprintf(['Now we are going to present you with a gravity map'...
            ' and with the mouse cursor please select the locations of %1.0d '...
            ' longitudinal profiles (north/south) with the mouse cursor'],nprofile(1)));
        waitfor(h(1))
    figure(1); % Open figure to assist in # of longitude and latitude profiles
        axis([min(lon) max(lon) min(lat) max(lat)])
    imagesc(lon,lat,grav)
    set(gca,'ydir','normal')
    xlabel('Longitude')
    ylabel('Latitude')
    title(sprintf('%s Gravity [mgal], Selecting Longitudinal Profiles (North/South)',name));
    shading interp
    view(0,90)
    colorbar 
    [lonX,latY] = ginput(nprofile(1));
            if nprofile(2) == 0
            latY= [];
            end
        end
    
          if nprofile(2) > 0 && comp{1}(8) == 1 % Select n longitudinal profiles
        h(1) = msgbox(sprintf(['Now we are going to present you with a gravity map'...
            ' and with the mouse cursor please select the locations of %1.0d '...
            ' Latitude profiles (east/west) with the mouse cursor'],nprofile(2)));
        waitfor(h(1))
    figure(1); % Open figure to assist in # of longitude and latitude profiles
    title(sprintf('%s Gravity [mgal], Selecting Latitudinal Profiles (East/West)',name));
    [lonX2,latY] = ginput(nprofile(2));
    clear lonX2
            if nprofile(1) == 0
            lonX= [];
            end
          end
        
    % Now we need to check to make sure that the ginput matches the n asked
    % for profiles and verify that these profiles are what the user wants
    % to use
        if nprofile(1) == length(lonX) && nprofile(2) == length(latY)
          
            figure(1); % Open figure to assist in # of longitude and latitude profiles
            axis([min(lon) max(lon) min(lat) max(lat)])
            imagesc(lon,lat,grav)
            set(gca,'ydir','normal')
            xlabel('Longitude')
            ylabel('Latitude')
            title(sprintf('%s Gravity [mgal], Verify Your Profiles',name));
            shading interp
            view(0,90)
            colorbar 
            title(sprintf('%s Gravity [mgal], Verify Your Profiles',name));
            if nprofile(1) > 0
                for ii = 1:nprofile(1)
                    line([lonX(ii) lonX(ii)],[min(lat) max(lat)],'LineWidth',3,'Color',[.8 .8 .8]);
                end
            end
            if nprofile(2) > 0
                for ii = 1:nprofile(2)
                    line([min(lon) max(lon)],[latY(ii) latY(ii)],'LineWidth',3,'Color',[.8 .8 .8]);
                end
            end
            button = questdlg('Are you happy with these profile selections (See figure)?',...
            'User Query','Yes, lets continue','No, I would like to pick new ones',...
            'I''m sick and tired of this program, abort','Yes, lets continue');
            if strcmp(button,'Yes, lets continue')
                comp{1}(3) = 0;
                close 1
            elseif strcmp(button,'No, I would like to pick new ones')
                
            else
                   h(1)=msgbox('Exiting Program, Have a super day');
                waitfor(h(1))
                return
            end
        end
    
    end          
        
   %% Step 5: Setting up Wave Vectors, Applying Filters Windows and
   % performing Transforms
    h(1)= msgbox(['We are now going to begin setting up the appropriate wave '...
        'vectors and perform the Fourier transforms. This may take a moment']);
        waitfor(h(1))
    N = length(lat);          
    if isempty(latY)  %Creating Xdir East/West WaveVectors, Cosine dependent
        
    else
        for ii = 1:nprofile(2)
        Kx(ii) = 1./(abs(lat(7)-lat(6))*111*1000*cosd(latY(ii))); 
        NNx(ii,:)=[-N/2+Kx(ii)/N : N/2]*Kx(ii)/N;
        end
        
    end
    
    if isempty(lonX) % Create Ydir North/South Wavevector
        
    else
        Ky = 1./(abs(lon(7)-lon(6))*111*1000); 
        NNy(:,1)=[-N/2+Ky/N : N/2]*Ky/N;
    end
        % Ask about filtering and perform filtering.
      button = questdlg('Would you like to Nyquist Filter the Topo Data? This may take a moment.',...
            'User Query','Yes, Filter it','No, don''t filter but continue',...
            'abort','Yes, Filter it');
            if strcmp(button,'Yes, Filter it') 
               sze= [length(lat),length(lon)];
               if exist('Ky','var') == 1;
                   cutoff = 0.5*Ky;
               else
                   cutoff = 0.5*Kx(1);
               end
               n=7;
               butter2 = lowpassfilter(sze,cutoff,n);
               topo = filter2(butter2,topo);
            elseif strcmp(button,'abort') 
                return
            end
        
        Wtopo = topo.*HAN; % Apply Windowing   
        
        if isempty(latY)  %Creating Xdir Fourier Transforms
        
    else
        for ii = 1:nprofile(2)
       Xtopo(ii,:)= (fftshift(fft(Wtopo(find(lat<latY(ii),1,'first'),:))));
       gravobX(ii,:) = gravw(find(lat<latY(ii),1,'first'),:);
        end
        
    end
    
    if isempty(lonX) % Creating Ydir (N/S) Fourier Transforms
       
    else
        for ii = 1:nprofile(1)
             Ytopo(:,ii) = (fftshift(fft(Wtopo(:,find(lon>lonX(ii),1,'first')))));
             gravobY(:,ii) = gravw(:,find(lon>lonX(ii),1,'first'));
        end
    end
%% Load and Moat system:
% Want to  find the gravity high, use this to select a min to either side
% then somewhere inbetween cut out the load and stitch two together, then
% run through the minimizer in step 6. Want to apply the transfer function
% to the whole profile though, just minimize the L&M (load and moat) profiles.
a = 1;
while a == 1
 button = questdlg('Would you like to use the Load & Moat system to minimize amplitude contribution to the analysis?',...
            'User Query','Yes, Use L&M','Please give me more info on L&M first',...
            'No, continue with analysis without L&M','Yes, Use L&M');
            if strcmp(button,'Yes, Use L&M')
                a = 0;
                
            elseif strcmp(button,'Please give me more info on L&M first')
                h=msgbox(['The Load and Moat system removes a gravity high associated'...
    ' with a topographical high or load. While performing analysis this load can skew the '...
    ' minimization routine due to amplitude differences. The important signal required for'...
    ' accurate minimization is in the flexure. The Load and Moat system works best for profiles'...
    ' with a clear gravitational load and associated flexural moat. If the L&M system is activated'...
    ' the load will be removed from the profile at some chosen interval between the load high'...
    ' and the moat. This will help reduce the amplitude error in the minimization and provide'...
    ' more accurate results'],'LOAD & MOAT SYSTEM README','modal');
    waitfor(h)
    
            else
                comp{1}(10)=0;
                a = 0;
            end
end


if comp{1}(10) == 1

    LtoM = 0.5; % LtoM needs to be set by user command, it is the percentage
    % from the max value (load) to the min value (moat) thus 100% or 1 would be
    % splicing the two moats together and totally cutting out the load.

    if isempty(latY)  %Creating Xdir Fourier Transforms
        
    else
            for ii = 1:nprofile(2)
            mxX(ii) = find(gravobX(ii,:)==max(gravobX(ii,:)),1,'first');
            mnX(ii,:) = [find(gravobX(ii,:)==min(gravobX(ii,1:mxX(ii))),1,'first'),...
            find(gravobX(ii,:)==min(gravobX(ii,mxX(ii):length(gravobX(ii,:)))),1,'first')];
            mnmxX(ii,:) = [round(LtoM*(mxX(ii)-mnX(ii,1))),round(LtoM*(mnX(ii,2)-mxX(ii)))];
            gravobXLM{ii}(1,:) = [gravobX(ii,1:(mnX(ii,1)+mnmxX(ii,1))),gravobX(ii,(mnX(ii,2)-mnmxX(ii,2)):length(gravobX(ii,:)))];
            end
        
    end
    
    if isempty(lonX) % Creating Ydir (N/S) Fourier Transforms
       
    else
        for ii = 1:nprofile(1)
              mxY(ii) = find(gravobY(:,ii)==max(gravobY(:,ii)),1,'first');
              mnY(ii,:) = [find(gravobY(:,ii)==min(gravobY(1:mxY(ii),ii)),1,'first'),...
              find(gravobY(:,ii)==min(gravobY(mxY(ii):length(gravobY(:,ii)),ii)),1,'first')];
              mnmxY(ii,:) = [round(LtoM*(mxY(ii)-mnY(ii,1))),round(LtoM*(mnY(ii,2)-mxY(ii)))];
              gravobYLM{ii}(:,1) = [gravobY(1:(mnY(ii,1)+mnmxY(ii,1)),ii),gravobY((mnY(ii,2)-mnmxY(ii,2)):length(gravobY(:,ii)),ii)];
        end
    end

end    
  %%  Step 6: Create Transfer Function and minimize he
  
he=linspace(0,100,15).*1000;
            % Ask user to use either regular rms approach or the
            % differentiation approach. Also ask for rms fit precision.
            while comp{1}(7) == 1
            prompt{2} = {['Regular Profile rms fit (similar to Watts'' method)'...
                ' or Differentiation Method: Type [''R'' for regular, or ''D'' for differentiation]'],...
                'Percentage Precision of RMS fit [0.5% recommended]:'};
            dlg_title{2} = 'Method of Analysis';
            def = {'R','0.5'};
            Q=inputdlg(prompt{2},dlg_title{2},num_lines{1},def); 
            P = cell2mat(Q(1)); Prec = str2num(Q{2})/100;
            if strcmp(P,'D') || strcmp(P,'R')
                comp{1}(7) = 0;
                break
            else
                h(1)=msgbox('Please Enter either ''D'' or ''R'' into the analysis method selection box');
                waitfor(h(1));
            end
            end

while comp{1}(4) == 1

    for jj = 1:length(he)
D = E*he.^3/(12*(1-v));

if isempty(latY)  %Creating East/West Transfer Functions and doing Fourier realm multiplication
    % and calculating rms for given h values
        
    else
        for ii = 1:nprofile(2)
       grav_transX(ii,:) = 2*pi*G*(rho_c-rho_w).*exp(-2*pi.*abs(NNx(ii,:)).*s).*...
           (1-  ((1+(D(jj).*((2*pi.*abs(NNx(ii,:))).^4)/(g*(rho_m-rho_c)))).^-1).*(exp(-2*pi.*abs(NNx(ii,:)).*d)));
       grav_X(ii,:) = grav_transX(ii,:).*Xtopo(ii,:); 
       grav_X(ii,:) = 1/1e-5.*(real(ifft(ifftshift(grav_X(ii,:)))));
       grav_X(ii,:) = grav_X(ii,:).*hanH;
       if comp{1}(10) == 1
           grav_XLM{ii}(1,:) = [grav_X(ii,1:(mnX(ii,1)+mnmxX(ii,1))),grav_X(ii,(mnX(ii,2)-mnmxX(ii,2)):length(grav_X(ii,:)))]; 
       end
        
       if strcmp(P,'D') && comp{1}(10)==1
        Dgrav_XLM{ii}(1,:) = diff(grav_XLM{ii}(1,:));
        DgravobXLM{ii}(1,:) = diff(gravobXLM{ii}(1,:));
        rmsX(ii,jj) = sqrt(sum((DgravobXLM{ii}(1,:)-Dgrav_XLM{ii}(1,:)).^2)./length(DgravobXLM{ii}(1,:))); 
       elseif strcmp(P,'D') 
        Dgrav_X(ii,:) = diff(grav_X(ii,:));
        DgravobX(ii,:) = diff(gravobX(ii,:));
        rmsX(ii,jj) = sqrt(sum((DgravobX(ii,:)-Dgrav_X(ii,:)).^2)./length(DgravobX(ii,:)));
       
       elseif comp{1}(10)==1
           rmsX(ii,jj) = sqrt(sum((gravobXLM{ii}(1,:)-grav_XLM{ii}(1,:)).^2)./length(gravobXLM{ii}(1,:)));
           
       else
           rmsX(ii,jj) = sqrt(sum((gravobX(ii,:)-grav_X(ii,:)).^2)./length(gravobX(ii,:)));
       end
        end
        
end
    
if isempty(lonX) % Creating North/South Transfer Functions and doing Fourier realm multiplication
        % and calculating rms for given h values
        
    else
        for ii = 1:nprofile(1)
            grav_transY(:,ii) = 2*pi*G*(rho_c-rho_w).*exp(-2*pi.*abs(NNy).*s).*...
                (1-  ((1+(D(jj).*((2*pi.*abs(NNy)).^4)/(g*(rho_m-rho_c)))).^-1).*(exp(-2*pi.*abs(NNy).*d)));
            grav_Y(:,ii) = grav_transY(:,ii).*Ytopo(:,ii); 
            grav_Y(:,ii) = 1/1e-5.*(real(ifft(ifftshift(grav_Y(:,ii)))));
            grav_Y(:,ii) = grav_Y(:,ii).*hanV;
        if comp{1}(10) == 1
            grav_YLM{ii}(:,1) = [grav_Y(1:(mnY(ii,1)+mnmxY(ii,1)),ii),grav_Y((mnY(ii,2)-mnmxY(ii,2)):length(grav_Y(:,ii)),ii)];    
        end
        
        if strcmp(P,'D') && comp{1}(10)==1
            Dgrav_YLM{ii}(:,1) = diff(grav_YLM{ii}(:,1));
            DgravobYLM{ii}(:,1) = diff(gravobYLM{ii}(:,1));    
            rmsY(jj,ii) = sqrt(sum((DgravobYLM{ii}(:,1)-Dgrav_YLM{ii}(:,1)).^2)./length(DgravobYLM{ii}(:,1)));
        elseif strcmp(P,'D')
            Dgrav_Y(:,ii) = diff(grav_Y(:,ii));
            DgravobY(:,ii) = diff(gravobY(:,ii));    
            rmsY(jj,ii) = sqrt(sum((DgravobY(:,ii)-Dgrav_Y(:,ii)).^2)./length(DgravobY(:,ii)));
        elseif comp{1}(10)==1
            rmsY(jj,ii) = sqrt(sum((gravobYLM{ii}(:,1)-grav_YLM{ii}(:,1)).^2)./length(gravobYLM{ii}(:,1)));
        else
            rmsY(jj,ii) = sqrt(sum((gravobY(:,ii)-grav_Y(:,ii)).^2)./length(gravobY(:,ii)));
       end
        end
end
    end
    
    if length(he) == 1
       comp{1}(4) = 0;
       comp{1}(5) = 0;
       
    else
        comp{1}(5) = 1;
    end
    
    while comp{1}(5) == 1
    % Now sum rms along profile dimensions to get a range of rms for
    % different he's. Determine if rms is within the error tolerance to
    % break loop
    
            if isempty(latY)
                
            else
                rmsX = sum(rmsX,1)./nprofile(2);
            end
            if isempty(lonX)
                
            else
                rmsY = sum(rmsY,2)./nprofile(1);
            end
                if exist('rmsX','var') == 1
                else
                    rmsX=0;
                end
                
                if exist('rmsY','var') == 1
                else
                    rmsY=0;
                end
            
            rms = (rmsX+rmsY')./(1 - isempty(rmsX) + 1 - isempty(rmsY));
   % Set h to one value if we are within error range and recompute 
           if abs(rms(find(rms==min(rms),1,'first')) - rms(find(rms==min(rms),1,'first')-1))/mean(rms) < Prec    
            he = he(rms == min(rms));
            rmstol = 100*abs(rms(find(rms==min(rms),1,'first')) - rms(find(rms==min(rms),1,'first')-1))/mean(rms);
            comp{1}(5) = 0; 
        
        else
            he = linspace(0.75*he(rms == min(rms)),1.25*he(rms == min(rms)),15);
            comp{1}(5) = 0;
        end
   end

end

%% STEP 7: Save progress and query user for finish
                if exist('gravobX','var') == 1
                else
                    gravobX=[];
                end
                
                if exist('grav_X','var') == 1
                else
                    grav_X=[];
                end
                if exist('gravobY','var') == 1
                else
                    gravobY=[];
                end
                
                if exist('grav_Y','var') == 1
                else
                    grav_Y=[];
                end
    a=1;
     h(1) = msgbox(sprintf(['The program has minimized the elastic crustal thickness to within'...
         ' an root mean squared error tolerance of %3.2f %%. The elastic thickness averaged over your'...
         ' profiles is %3.2f km. You will be provided a chance to save and plot this information.'],rmstol,he/1000));
     waitfor(h(1))

     button = questdlg('Save and Plot Data',...
            'User Query','Save & Plot Data','Just Save my Data',...
            'Quit','Save & Plot Data');
        if strcmp(button,'Save & Plot Data') 
         
        save(sprintf('%s.mat',name),'he','name','grav_X','grav_Y','gravobX','gravobY','lonX','latY','nprofile','-append');
        
        elseif strcmp(button,'Just Save my Data') 
        save(sprintf('%s.mat',name),'he','grav_X','grav_Y','gravobX','gravobY','-append');
        comp{1}(6)=0;
        close all
        a=0;
        else
             h(1)=msgbox('Exiting Program, Have a super day');
                waitfor(h(1))
            return
        end
%% PLOT DATA
if comp{1}(6)==1
close all

GRAVOB=[gravobX ; gravobY'];
GRAVPR=[grav_X ; grav_Y'];
if nprofile(1) > 0 
   for qq = 1:nprofile(1)
       grid(qq,:) = lon;
       TITLE(qq) = lonX(qq);
   end
end

if nprofile(2) > 0
    for qq = 1: nprofile(2)
        grid(qq+nprofile(1),:)= lat';
        TITLE(qq+nprofile(1)) = latY(qq);
    end
end

figure(1)
for ii = 1: length(grid(:,1))
subplot(length(grid(:,1)),1,ii)
plot(grid(ii,:),GRAVOB(ii,:),grid(ii,:),GRAVPR(ii,:),'r')
title(sprintf('Profile for %3.2f',TITLE(ii)))  
end

end

button = questdlg('Perform another run',...
            'User Query','Let me try more Profiles, (close figure to proceed)','No, I am all done',...
            'Quit','Let me try more Profiles, (close figure to proceed)');
        if strcmp(button,'Let me try more Profiles, (close figure to proceed)') 
            comp{1}(1:20)=1;
            if a == 1;
            h(1)=figure(1);
            waitfor(h(1));
            end
            elseif strcmp(button,'No, I am all done') 
             h(1)=msgbox('Exiting Program, Have a super day');
                waitfor(h(1))
            return
        else
             h(1)=msgbox('Exiting Program, Have a super day');
                waitfor(h(1))
            return
        end
end

 