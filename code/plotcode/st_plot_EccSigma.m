function data = st_plot_EccSigma(df,binsize)
% rmPlotEccSigma - plot sigma versus eccentricity in selected ROI
%
% data = rmPlotEccSigma(v, lr, [fwhm=0], [plot_position_variance], [plotFlag=1]);
%
%  INPUT
%   v: view
%   lr: left (1) or right (2) if FLAT view is used
%   fwhm: plot sigma (0) or fwhm (1)
%   plot_pv: plot 1D position variance (1), 2D position variance (2), sigma in mm (3) or not
%   (0)
%
% Uses data above view's co-threshold. It does makes sure that this
% threshold is above that used for the fine-search stage - so not to
% include blurred data. 
%
% 2007/08 SOD: ported from my scripts.
% 2008/11 KA: included the ability to plot position variance
% 2008/11 RAS: added a separate plot flag.
% if ~exist('v','var') || isempty(v),
%     v = getCurView;
% end
% 
% if strcmp(v.viewType,'Flat')
%     if notDefined('lr'), lr=1;  end
% end
% 
% if notDefined('fwhm')
%     fwhm=0;
% end
% 
% if notDefined('plotFlag'),	plotFlag = 1;		end
% if notDefined('plot_pv'),    plot_pv=0;			end
% 
% %% load stuff
% % load retModel information
% try
%     rmFile   = viewGet(v,'rmFile');
%     rmParams = viewGet(v,'rmParams');
% catch %#ok<CTCH>
%     error('Need retModel information (file)');
% end
% % load ROI information
% try
%     roi.coords = v.ROIs(viewGet(v,'currentROI')).coords;
%     titleName  = v.ROIs(viewGet(v,'currentROI')).name;
% catch %#ok<CTCH>
%     error('Need ROI');
% end
% 
% % load all coords
% if strcmp(v.viewType,'Flat')
%     if lr ==1
%         allCrds  = viewGet(v,'coords','left');
%     elseif lr==2
%         allCrds  = viewGet(v,'coords','right');
%     end
% else
%     allCrds  = viewGet(v,'coords');
% end


% get data
% rmData   = load(rmFile);
% [tmp, roi.iCrds] = intersectCols(allCrds,roi.coords);
% if plot_pv == 1 % 1D position variance plot
%     getStuff = {'ecc','pol','sigma','varexp','pv'};
% elseif plot_pv == 2 % 2D position variance plot
%     getStuff = {'ecc','pol','sigma','varexp','pv2'};
% elseif plot_pv == 3 % sigma in mm
%     getStuff = {'ecc','pol','sigma','varexp','s_mm'};
% else
%     getStuff = {'ecc','pol','sigma','varexp'};
% end
titleName = df.roiName2(1);
getStuff = {'ecc','pol','sigma','varexp'};
getStuff = {'ecc','sigma','cv_varexp'};

% binsize = 1;
sigmaMax = 15;
eccMax = 12;
MarkerSize = 8;
plot_pv = 0;
for m = 1:numel(getStuff)
    tmp = df.(getStuff{m});
    roi.(getStuff{m}) = tmp;
end
% if plot_pv == 2 % 2D position variance plot
%     roi.pv=roi.pv2;
% elseif plot_pv == 3 % plot sigma in mm
%     roi.pv = roi.s_mm;
% end
% if fwhm ==1
%     roi.sigma = sd2fwfm(roi.sigma);
% end

% % bin size (eccentricity range) of the data
% if max([rmParams.stim(:).stimSize])>4,
%     binsize = 1;
%     %binsize = 1;
% else
%     %binsize = 0.25;
%     binsize = 0.5;
% end

%%--- thresholds
% take all data that would be processed in the 'search/fine' stage
% do not take data below this threshold because this is from the
% 'grid/coarse' stage and involves blurring.
thresh.varexp  = 0;

% take all data within the stimulus range, and decrease it by a small 
% amount to be more conservative.
% thresh.ecc = [0 eccMax] + [.5 -2];% * binsize/2; 
% thresh.ecc = [0 eccMax] + [.5 -2];% * binsize/2; 
thresh.ecc = [0 eccMax] + ([1 -1]* binsize/2); 

% thresh.ecc = viewGet(v, 'mapclip') + [0.5 -0.5];
% thresh.ecc = [0 max([rmParams.stim(:).stimSize])] + [0.5 -0.5];% * binsize/2; 

% thresh.ecc = [0 30];
% basically no sigma threshold
thresh.sig = [0.01 sigmaMax-0.5]; 

%--- plotting parameters
xaxislim = [0 eccMax];

% find useful data given thresholds
ii = roi.cv_varexp > thresh.varexp & ...
     roi.ecc > thresh.ecc(1) & roi.ecc < thresh.ecc(2) & ...
     roi.sigma > thresh.sig(1) & roi.sigma < thresh.sig(2);

% weighted linear regression:
roi.p = linreg(roi.ecc(ii),roi.sigma(ii),roi.cv_varexp(ii));
roi.p = flipud(roi.p(:)); % switch to polyval format
xfit = thresh.ecc;
yfit = polyval(roi.p,xfit);

% bootstrap confidence intervals
if exist('bootstrp','file') 
    B = bootstrp(1000,@(x) localfit(x,roi.ecc(ii),roi.sigma(ii),roi.cv_varexp(ii)),(1:numel(roi.ecc(ii))));
    B = B';
    pct1 = 100*0.05/2;
    pct2 = 100-pct1;
    b_lower = prctile(B',pct1);
    b_upper = prctile(B',pct2);
    keep1 = B(1,:)>b_lower(1) &  B(1,:)<b_upper(1);
    keep2 = B(2,:)>b_lower(2) &  B(2,:)<b_upper(2);
    keep = keep1 & keep2;
    b_xfit = linspace(min(xfit),max(xfit),100)';
    fits = [ones(100,1) b_xfit]*B(:,keep);
    b_upper = max(fits,[],2);
    b_lower = min(fits,[],2);
end

% if plot_pv
%     roi.p2 = linreg(roi.ecc(ii),roi.pv(ii),roi.varexp(ii));
%     roi.p2 = flipud(roi.p2(:)); % switch to polyval format
%     y2fit = polyval(roi.p2,xfit);
% end

% output struct
data.xfit = xfit(:);
data.yfit = yfit(:);
data.x    = (thresh.ecc(1):binsize:thresh.ecc(2))';
data.y    = nan(size(data.x));
data.ysterr = nan(size(data.x));
data.roi = roi;

% if plot_pv
%     data.y2fit = y2fit(:);
%     data.y2    = nan(size(data.x));
%     data.y2sterr = nan(size(data.x));
% end

% plot averaged data
for b=thresh.ecc(1):binsize:thresh.ecc(2),
    bii = roi.ecc >  b-binsize./2 & ...
          roi.ecc <= b+binsize./2 & ii;
    if any(bii),
        % weighted mean of sigma
        s = wstat(roi.sigma(bii),roi.cv_varexp(bii));
        % store
        ii2 = find(data.x==b);
        data.y(ii2) = s.mean;
        data.ysterr(ii2) = s.sterr;
        
        if plot_pv
            % weighted mean of position variance
            s = wstat(roi.pv(bii),roi.varexp(bii));
            % store
            data.y2(ii2) = s.mean;
            data.y2sterr(ii2) = s.sterr;
        end
    else
       fprintf(1,'[%s]:Warning:No data in eccentricities %.1f to %.1f.\n',...
            mfilename,b-binsize./2,b+binsize./2);
    end;
end;

% plot if requested


return;


function B=localfit(ii,x,y,ve)
B = linreg(x(ii),y(ii),ve(ii));
B(:);
return

