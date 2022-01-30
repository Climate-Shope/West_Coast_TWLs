%Second part in wave reconstruction code. For use after part M3a_A
%Useful for offloading to another matlab-running machine as code will take
%days potentially to complete depedning on the number of reconstruction
%points
load('Douglas_TS_Input');

dirout=[]; %Change This to desired output location 


%% TS RECONSTRUCTION 

global PosEscalares
global PosDireccionales 

VISIBILITY = 'on'; 
OVERWRITE = 0; 

% -------------------------------------------------------------------------

RUN_INDS = [ 1:numel(ii_grid_array)] ; % RUNNING ALL TRANSECTS

% -------------------------------------------------------------------------





for ii =RUN_INDS 
    
    idprofile = ii_grid_array(ii);

    OBJECTID = idprofile; 
    disp([num2str(ii),'/',num2str(Nact_trs), ' - prf ',num2str(idprofile),'/',num2str(Ntr)]) 
    
    %% select appropriate GOW point data to use
    appro=GRIDid(idprofile);
    eval(['ALLDATA=ALLDATA_' num2str(appro) ';']);
    eval(['classification=classification_' num2str(appro) ';']);
    eval(['positions=positions_' num2str(appro) ';']);
    
    %%%% identify conditons to remove

    dir_idx=find(ALLDATA.X(:,3) > 0 & ALLDATA.X(:,3) < 180);
    
    
    
    %% skip if mat file already exists 
    fname= [dirout,filesep,'WD4R_prf',dec2base(OBJECTID,10,4)]; 
    if exist([fname,'.mat'],'file')==2 && 1-OVERWRITE, continue, end 
    
    %% DO NOT CONSIDER NANs in propagations 
    IND = find(1-isnan(Hsinterp(:,idprofile))); 
    if isempty(IND), continue, end % all nans 
    
    wave_data       = struct(); 
    wave_data.Hs    = Hsinterp(IND,idprofile);  % m 
    wave_data.T     = Tpinterp(IND,idprofile);  % secs 
 
    
    if all(isnan(wave_data.Hs)), 
        warning(['ALL Hs data is NAN - ',num2str(OBJECTID)])
        disp('skipping...') 
        continue 
    end 
    
    %% reconstruct Hs at head of profile: 

    
    datos = ALLDATA.X; % data offshore
    if numel(classification.bmus)>534e3 
      datos(indsremove,:)=[]; 
    end
    
    % do not consider NAN - USE IND 
    casos_datos  = classification.final(IND,:);  % not normalized
    
    PosEscalares = classification.escalar; 
    PosDireccionales = classification.direccional; 
    
    t1 = now;    

    
    % Hs 
    Propagaciones   = [wave_data.Hs(:)]'; % propagated
    
    
    
    %%
    %Reconstruction functions from Camus et al. 2011
    [Hs]=fun_reconstruct_RBF(datos,casos_datos,Propagaciones);
    get(gcf), title(['Hs - Profile: ',dec2base(idprofile,10,4)])
    saveas(gcf, [dirout,filesep,'rbf_Hs_',dec2base(idprofile,10,4),'.png'])
    hsremove = Hs<0.1; 
    Hs(hsremove)=0.1; 

    % Tm %Reconstruction functions from Camus et al. 2011
    [Tp]=fun_reconstruct_RBF(datos,casos_datos,[wave_data.T(:)]'); 
    saveas(gcf, [dirout,filesep,'rbf_T_',dec2base(idprofile,10,4),'.png'])
    get(gcf), title(['T - Profile: ',dec2base(idprofile,10,4)])
    Tp(hsremove)=1.5;
    
    Hs(dir_idx)=0.1;
    Tp(dir_idx)=1.5;


    t2 = now; 
    duration_run = datevec(t2-t1); 
    
    Hsnnb  = wave_data.Hs(classification.bmus); 
    Tmnnb  = wave_data.T(classification.bmus); 
	
    if numel(classification.bmus)>534e3
      Hsnnb(indsremove,:)=[]; 
      Tmnnb(indsremove,:)=[];   
    end
    
     %%%%% Figure
    figure ('visible',VISIBILITY,'position',[ 560   500   867   368])
    hold on 
    plot(time, Hsnnb,'+','markersize',1,'color',[1 1 1].*0.2) 
    plot(time, datos(:,1),'b'), 
    plot(time, Hs,'r') 
    legend('Hs nnb', 'Hs offshore', 'Hs propagated')
    grid on, box on 
    title(['Prf. ',dec2base(idprofile,10,4)])
    ylabel('m')   
    tlabel('x') 
    axis tight 
    
    file=[dirout,filesep,'ts_',num2str(OBJECTID)]; 
    set(gcf,'PaperPositionMode','auto','InvertHardcopy','on')
    print(gcf,'-dpng','-r100',file)
        
    close all 
    
    %% save data

     Hs = Hs(:); 
    T  = Tp(:); 

    eval(['x0 = ' NAME_GRD_transects '_Contour_pts.X(idprofile);']);
    eval(['y0 = ' NAME_GRD_transects '_Contour_pts.Y(idprofile);']);
    eval(['z0 = ' NAME_GRD_transects '_Contour_pts.Z(idprofile);']);
    OBJECTID = dec2base(idprofile,10,4); 
    nname=['WD4R_prf',dec2base(idprofile,10,4)];
    eval([nname '=m;']);
    save(sprintf('%s.mat', fname),'Hs','T','x0','y0','z0','OBJECTID');
    eval(['clear ' nname]);

   
end

