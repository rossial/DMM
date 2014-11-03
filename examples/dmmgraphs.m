% ---------------------------------------------------
% DMMgraphs read DMM output and makes graphs of
% model parameters and unobsevables
% A. Rossi, November 2014
%

% ---------- to be set by the user  ---------------------------------------
path = 'h:\arossi\dmm\nile\';   % location of the .nml file e.g. nile.nml
file = 'nile';                  % name of the nml file

freq      = 1;                  % 1 annual, 4 quarterly, 12 monthly
startyear = 1985;               
startper  = 1;
space     = 8;                  % spacing for the tme labels 
NHIST     = 200;                % to compute percentiles
% -------------------------------------------------------------------------

% Load metadata, priors, observations
close all
fid = fopen([path,file,'.PRI']);
nt     = str2double(fscanf(fid,'%s',1));
np     = str2double(fscanf(fid,'%s',1));
maxvec = str2double(fscanf(fid,'%s',1));
maxhyp = str2double(fscanf(fid,'%s',1));
nf     = str2double(fscanf(fid,'%s',1));
nz     = str2double(fscanf(fid,'%s',1));
seed   = fscanf(fid,'%s',1);
nx     = str2double(fscanf(fid,'%s',1));
ny     = str2double(fscanf(fid,'%s',1));
nobs   = str2double(fscanf(fid,'%s',1));
nv     = str2double(fscanf(fid,'%s',1));

NS = zeros(nv,1);
for i = 1:nv
    NS(i) = str2double(fscanf(fid,'%s',1));
end
estimation = fscanf(fid,'%s',1);

prior    = zeros(nt,4);
psiprior = zeros(maxvec,maxhyp);
tipo     = cell(nt,1);
name     = cell(nt+np,1);
for i = 1:nt
    for j =1:4
        prior(i,j) = str2double(fscanf(fid,'%s',1));
    end
    tipo(i) = cellstr(fscanf(fid,'%s',1));
    name(i) = cellstr(['\theta_{' int2str(i) '}']);
end
nstate = zeros(maxvec,1);
Sdyn   = cell(maxvec,1);
for i = 1:maxvec
    nstate(i,1) = str2double(fscanf(fid,'%s',1));
    for j =1:maxhyp
        psiprior(i,j) = str2double(fscanf(fid,'%s',1));
    end
    Sdyn(i) = cellstr(fscanf(fid,'%s',1));
end

% data + forecasts
lab = plotlabsmall(startyear,startper,nobs+nf,freq);
ind = fliplr(nobs+nf:-space:1);
if nf > 0
    fore = load([path,file,seed,'.FST']);
end
yk  = zeros(nobs+nf,ny+nz);
IYK = zeros(nobs,ny+1);
for i = 1:nobs
    for j = 1:ny+nz
        yk(i,j) = str2double(fscanf(fid,'%s',1));
    end
end
fclose(fid);

% track missings
for i = 1:nobs
    K = 0;
    for j = 1:ny
        if yk(i,j) ~= -99999
            K = K+1;
            IYK(i,j) = j;
        end
        IYK(i,ny+1) = K;
    end
end
nmiss = ny*nobs-sum(IYK(1:nobs,ny+1));

% theta and psi
par = load([path,file,seed,'.PAR']);
theta = par(:,1:nt);
if np > 0
    psi   = par(:,nt+1:nt+np);
    for i = 1:np
        name(nt+i) = cellstr(['\psi_' int2str(i)]);
    end
end

% Beta marginal prior for psi
psimarg = zeros(np,2);
n = 0;
for k = 1:maxvec
    a0 = sum(psiprior(k,1:nstate(k)));
    for j = 1:nstate(k)-1
        psimarg(n+1,1) = psiprior(k,j);
        psimarg(n+1,2) = a0-psiprior(k,j);
        n = n+1;
    end
end

% Plot posterior and prior pdfs
[nr,nc] = numsubplot(nt+np);
for i = 1:nt
    if prior(i,3)-prior(i,4) ~= 0
        [prc,iqr] = myprc(theta(:,i),NHIST,.01,.99);        
        hh = 0.79*length(theta(:,i))^(-1/5)*iqr;
        step = (prc(2)-prc(1))/200;
        [x,fx] = ker1boun(theta(:,i),hh,step,prc(1),prc(2));        
        if tipo{i} == 'BE'
            s = prior(i,1);
            v = prior(i,2);
            a = prior(i,3);
            b = prior(i,4);
            y = [a:(b-a)/200:b]';            
            fy = dmmprior((y-a)/(b-a),s,v,tipo{i});            
            fy = fy/(b-a);
        elseif tipo{i} == 'IG'
            s = prior(i,1);
            v = prior(i,2);
            Med = s/(v-2);
            y  = [0:5*Med/200:5*Med]';
            fy = dmmprior(y,s,v,tipo{i});            
        elseif tipo{i} == 'NT'
            med = prior(i,1);
            va  = prior(i,2);
            lb  = prior(i,3);
            ub  = prior(i,4);
            y   = (lb:(ub-lb)/200:ub)';
            fy  = dmmprior(y,med,va,tipo{i});            
%           plb = normcdf(lb,med,sd);
%           pub = normcdf(ub,med,sd);
%           fy  = normpdf(y,med,sd)/(pub-plb);
        end
        subplot(nr,nc,i)
        plot(x,fx,'b',y,fy,'r--')
        axis([prc(1) prc(2) min(fx) max(fx)])
        title(name{i,:})
        grid on
    end
end
for i = 1:np
    [prc,iqr] = myprc(psi(:,i),NHIST,.01,.99);    
    hh = 0.79*length(psi(:,i))^(-1/5)*iqr;
    step = (prc(2)-prc(1))/200;
    [x,fx] = ker1boun(psi(:,i),hh,step,prc(1),prc(2));
    
    s  = psimarg(i,1);
    v  = psimarg(i,2);
    y  = (1/200:1/200:1-1/200)';
    fy = dmmprior(y,s,v,'BE');
    
    subplot(nr,nc,i+nt)
    plot(x,fx,'b',y,fy,'r--')
    axis([prc(1) prc(2) min(fx) max(fx)])
    grid on
    title(name{i+nt,:})
end

% Data + missings + forecasts
prc = nan*zeros(nobs+nf,ny,2);
if nmiss > 0
    miss = load([path,file,seed,'.MIS']);
    K = 0;
    for i = 1:nobs
        for j = 1:ny
            if yk(i,j) == -99999
                K = K+1;
                yk(i,j)    = mean(miss(:,K));
                [ppp,iqr]  = myprc(miss(:,K),NHIST,.10,.90);
                prc(i,j,1) = ppp(1);
                prc(i,j,2) = ppp(2);
            end
        end
    end
end
if nf > 0
    for j = 1:ny
        yk(nobs+1:nobs+nf,j) = mean(fore(:,nf*(j-1)+1:nf*j))';
        K = 0;
        for k=nf*(j-1)+1:nf*j
            K = K+1;
            [ppp,iqr] = myprc(fore(:,k),NHIST,.10,.90);
            prc(nobs+K,j,1) = ppp(1);
            prc(nobs+K,j,2) = ppp(2);
        end
    end
end
[nr,nc] = numsubplot(ny);
for i = 1:ny
    subplot(nr,nc,i)
    plot(1:nobs+nf,yk(:,i),'k',1:nobs+nf,prc(:,i,1),'b*',1:nobs+nf,prc(:,i,2),'b*')
    axis([1 nobs+nf min([yk(:,i); squeeze(prc(:,i,1))]) max([yk(:,i);squeeze(prc(:,i,2))])])
    grid on
    title(['Series No ' num2str(i)])
    if nf > 0
        hold on;
        plot(nobs:nobs+1,[-10^10;10^10],'k')
        hold off;
    end
    set(gca,'xtick',ind);
    set(gca,'xticklabel',lab(ind));
end

% State inference
bstate  = load([path,file,seed,'.UNB']);
G       = size(bstate,1);
[nr,nc] = numsubplot(nx);
for i = 1:nx
    stato = bstate(:,(i-1)*nobs+1:i*nobs);
    if nf > 0
        stato = [stato fore(:,nf*ny+nf*(i-1)+1:nf*ny+nf*i)];
    end 
    statomed = mean(stato);
    for j = 1:nobs+nf
        [prc(j,:),iqr] = myprc(stato(:,j),NHIST,.10,.90);
    end    
    subplot(nr,nc,i)
    plot(1:nobs+nf,statomed,'k',1:nobs+nf,prc(:,1),'b--',1:nobs+nf,prc(:,2),'b--')
    axis([1 nobs+nf min(prc(:,1)) max(prc(:,2))])
    grid on
    title(['State No ' num2str(i)])
    if nf > 0
        hold on;
        plot(nobs:nobs+1,[-10^10;10^10],'k')
        hold off;
    end
    set(gca,'xtick',ind);
    set(gca,'xticklabel',lab(ind));    
end

% Innovations
inn  = load([path,file,seed,'.INN']);  % G x nser*T
[nr,nc] = numsubplot(ny);
for i = 1:ny
    minn = mean(inn(:,nobs*(i-1)+1:nobs*i))';
    subplot(nr,nc,i)
    plot(1:nobs,minn,'k')
    axis([1 nobs min(minn(1:nobs,i)) max(minn(1:nobs,i))])
    set(gca,'xtick',ind);
    set(gca,'xticklabel',lab(ind));
    grid on
    title(['Innovations (full line) and series No ' num2str(i)])
end

% Discrete latent variables
if nv > 0
    Z = load([path,file,seed,'.DIS']);
    if nf > 0
        Z(:,nobs+1:nobs+nf) = fore(:,end-nf+1:end);
    end    
    medS = zeros(nobs+nf,nv);    
    for j = 1:nobs+nf
        for k=1:NS
            medS(j,k) = numel(find(Z(:,j)==k));
        end
    end
    medS = medS/G;
    [nr,nc] = numsubplot(NS);
    for i = 1:NS
        subplot(nr,nc,i)
        plot(1:nobs+nf,medS(:,i),'k')
        axis([1 nobs+nf 0 1])
        grid on
        title(['Posterior mean of S' num2str(i)])
        set(gca,'xtick',ind);
        set(gca,'xticklabel',lab(ind));
    end
end