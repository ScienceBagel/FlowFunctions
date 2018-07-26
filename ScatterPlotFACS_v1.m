function [h, LRx, LRy] = ScatterPlotFACS_v1(X,varargin)
%
% Creates scatter plot with colorcoded density. 
%
% h = ScatterPlotFACS(X) - if X is a vector, the program produces a
% Histogram plot. If X is a double column matrix then the Scatter plot is
% done of the second column vs. the first.
%
% h = ScatterPlotFACS(X, Y) - X and Y should be vectors for the Scatter
% plot
% By default the plot is done in the double arcsinh scale. The linear range
% is determined automatically from 95% of negative values. 0.02% of
% of outliers are placed at the borders of the plot.
%
% All of these and some other parameters can be modified in the most general use of the
% function:
%
% h = ScatterPlotFACS(X, Y, 'Param', ParamValue, ..., 'Param', ParamValue),
% or
% h = ScatterPlotFACS(X, [], 'Param', ParamValue, ..., 'Param', ParamValue),
%
% The allowed parameter pairs are:
%
% 'LogXscale' - 1 for arcsinh scale on X axis, 0 for linear plot
% 'LogYscale' - 1 for arcsinh scale on Y axis, 0 for linear plot
%
% 'LinearRangeFraction' - a number between 0 and 1 that determines a fraction
%                         of negative values used to determine the linear range.
%                         The default value is 0.05.
%
% 'LinearRangeX' and 'LinearRangeY' - real numbers that determine linear
%                                      range in asinh plots.
%
% 'OutlierFraction' - a number between 0 and 1 that determines a fraction
%                         of ouliers that are going to be plotted on the axes borders.
%                         The default value is 2e-4.
%
% 'Xlim' and 'Ylim' followed by two-number vector, e.g. [-100 1e5] set hard
%                   limits on the graph axes. Otherwise the limits are determined from
%                   OutlierFraction
%
% 'Nbins'        followed by a natural number: the number of bins for
%                  Histogram plot and also the number of colors for Color coding of density
%               (decreasing Nbins might speed up the plots)
%
% 'MaxPointsToPlot' followed by a natural number: the maximal number of points to plot. Reduce
%                   the number to speed up the plot. A random sample of
%                   MaxPointsToPlot of the data points will be shown.
%
% 'PlotGate'        followed by the two-column array of (X, Y) pairs of
%                   polygon vertices for polygonal or rectangular gate, or two pair vector
%                   for the range gate in histogram plot.
LRx = Inf;
LRy = Inf;


if size(varargin,2) == 1;
    varargin = varargin{:,1};
end;

HistPlot = 0; % flag for Histogram plot
LogXscale = 0; % flag for Asinh scale on X axis
LogYscale = 0; % flag for Asinh scale on Y axis
LinearRangeFraction = 0.05;
OutlierFraction = 5e-5;
Nbins = 30;
MaxPointsToPlot = 50000;
PlotGate =0;

%tranposing X array if improperly oriented
if (size(X, 2)>2) & (size(X, 1) > 3), 
    X = X';
end;

% The first parameter can be two column (two rows) array, in which case
% the first column is X and the second column Y. In this case, if more
% parameters are used then the second parameter should be an empty vector
% ([]). If the first two parameters are vectors then the first is X and the
% second Y. Finally, if X is a vector and there are no other parameters or
% the second parameter is empty ([]) then the program assumes Histogram
% plot.

if isempty(varargin) | isempty(varargin{1}),
    if size(X, 2) == 2,
        Y = X(:, 2);
        X = X(:, 1);       
    elseif size(X, 2) == 1,
        HistPlot = 1; % do histogram
    else
        error('X should contain no more than two columns!')
    end;
else % now to normal scenario: two parameters X and Y
    if isvector(X),
        X=X(:); % make it a column vector (already done at the beginning but just in case)
    
    if isnumeric(Y),
        if isvector(Y),
            Y=Y(:); % make it a column vector
            if (length(X)~=length(Y)),
                 error('X and Y should be vectors of same length!')
            end;
        else
            error('In two parameter scatter plot X and Y should be both vectors!')
        end;
    end;
    
    elseif size(X, 2) == 2,
        Y = X(:, 2);
        X = X(:, 1); 
    else
        error('In two parameter scatter plot X and Y should be both vectors!')
    end;
    

end;

if HistPlot
    varargin = Y;
end;
j1=find(strcmp(varargin, 'PlotGate'));
if length(j1)>1,
    error('Too many PlotGate inputs!')
elseif length(j1)==1,
    Gate = varargin{j1+1};
    PlotGate = 1;
end;


j1=find(strcmp(varargin, 'LogXscale'));
if length(j1)>1,
    error('Too many LogXscale inputs!')
elseif length(j1)==1,
    LogXscale = varargin{j1+1};
end;

j1=find(strcmp(varargin, 'LogYscale'));
if length(j1)>1,
    error('Too many LogYscale inputs!')
elseif length(j1)==1,
    LogYscale = varargin{j1+1};
end

j1=find(strcmp(varargin, 'LinearRangeFraction'));
if length(j1)>1,
    error('Too many LinearRangeFraction inputs!')
elseif length(j1)==1,
    LinearRangeFraction = varargin{j1+1};
end

j1=find(strcmp(varargin, 'LinearRangeX'));
if length(j1)>1,
    error('Too many LinearRangeX inputs!')
elseif length(j1)==1,
    LRx = varargin{j1+1};
elseif LogXscale %determine from data
    Xneg = X(X < 0);
    if isempty(Xneg),
        LRx = 100; % basically  a random number
    else
        Xneg = sort(Xneg, 'ascend');
        LRx = abs(Xneg(ceil(length(Xneg)*LinearRangeFraction)));
    end;
    clear Xneg;
end;

j1=find(strcmp(varargin, 'LinearRangeY'));
if length(j1)>1,
    error('Too many LinearRangeY inputs!')
elseif length(j1)==1,
    LRy = varargin{j1+1};
elseif ((~HistPlot)&LogYscale) %determine from data
    Yneg = Y(Y < 0);
    if isempty(Yneg),
        LRy = 100; % basically a random number
    else
        Yneg = sort(Yneg, 'ascend');
        LRy = abs(Yneg(ceil(length(Yneg)*LinearRangeFraction)));
    end;
    clear Yneg;
end;

j1=find(strcmp(varargin, 'OutlierFraction'));
if length(j1)>1,
    error('Too many OutlierFraction inputs!')
elseif length(j1)==1,
    OutlierFraction = varargin{j1+1};
end

j1=find(strcmp(varargin, 'Xlim'));
if length(j1)>1,
    error('Too many Xlim inputs!')
elseif length(j1)==1,
    Xlim = varargin{j1+1};
else %find plot limits by cutting out outliers
    Xlim = FindPlotLimits(X, OutlierFraction);
end;
set(gca,'Xlim',Xlim)

j1=find(strcmp(varargin, 'Ylim'));
if length(j1)>1,
    error('Too many Ylim inputs!')
elseif length(j1)==1,
    Ylim = varargin{j1+1};
    set(gca,'Ylim',Ylim)
elseif ~HistPlot %find plot limits by cutting out outliers
    Ylim = FindPlotLimits(Y, OutlierFraction);
    set(gca,'Ylim',Ylim)
end;



j1=find(strcmp(varargin, 'Nbins'));
if length(j1)>1,
    error('Too many Nbins inputs!')
elseif length(j1)==1,
    Nbins = varargin{j1+1};
end;

j1=find(strcmp(varargin, 'MaxPointsToPlot'));
if length(j1)>1,
    error('Too many MaxPointsToPlot inputs!')
elseif length(j1)==1,
    MaxPointsToPlot = varargin{j1+1};
end;

if length(X) > MaxPointsToPlot,
    JJ = round(1 + (length(X)-1)*rand(MaxPointsToPlot, 1));
    X = X(JJ);
    X = X(:);
    if ~HistPlot,
        Y = Y(JJ);
        Y = Y(:);
    end;
    clear JJ;
end;

% start plotting
% find data that are not outliers
Jx = (X>Xlim(1))&(X<Xlim(2));
if ~HistPlot,
    Jy = (Y>Ylim(1))&(Y<Ylim(2));
end;

if ~HistPlot,
    if PlotGate,
        if (size(Gate, 1) < size(Gate, 2)),
            Gate = Gate';
        end;
        if (size(Gate, 2)~=2)|(size(Gate, 1) < 3),
            error('Polygonal gate should be a two-column matrix of at least three rows');
        end;
        
        Xgate = [Gate(:, 1); Gate(1, 1)];
        Ygate = [Gate(:, 2); Gate(1, 2)];
    end
    
    if (LogXscale & LogYscale), % double asinh plot
        XX = asinh(X/LRx);
        YY = asinh(Y/LRy); 
        if PlotGate,
            Xgate = asinh(Xgate/LRx);
            Ygate = asinh(Ygate/LRy);
        end;
    elseif ((~LogXscale) & (~LogYscale)), % double linear plot
        XX = X;
        YY = Y;
    elseif (LogXscale & (~LogYscale)), % asinh X plot
        XX = asinh(X/LRx);
        YY = Y;
        if PlotGate,
            Xgate = asinh(Xgate/LRx);
        end;
    else % asinh Y plot
        XX = X;
        YY = asinh(Y/LRy);
        if PlotGate,
            Ygate = asinh(Ygate/LRy);
        end;
    end;  
        
    h = PlotColorCoded(XX, YY, Jx, Jy, Nbins);
    if PlotGate,
        hold on;
        hline = line(Xgate, Ygate);
        set(hline, 'Color', 'black','LineWidth',2)
        hold off;
    end
    figure(gcf)
else % if plotting Histogram
    if PlotGate, % assume range gate
        if (size(Gate, 1)>2)| (size(Gate, 2)>2),
            error('Range gate should contain no more than two numbers!!!'), 
        end;
        Xgate = Gate;
    end;
    if LogXscale, % asinh plot
        XX = asinh(X/LRx); 
        if PlotGate, %assume range gate
            Xgate = asinh(Xgate/LRx);
        end;
    else %linear plot
        XX = X;
    end;
    [HX, Xc] = hist(XX, Nbins);
    plot(Xc, HX,'LineWidth', 1.5);

    if PlotGate,
        Xgate = Xgate(:);
        Ygate = max(HX)/5*ones(2, 1); 
        hold on
        bartickle = 100;
        hline(1) = line(Xgate, Ygate);
        
        hline(2) = line([Xgate(1) Xgate(1)] , [Ygate(1)-bartickle Ygate(1)+bartickle]);
        hline(3) = line([Xgate(2) Xgate(2)] , [Ygate(2)-bartickle Ygate(2)+bartickle]);
         Width = 2;
        set(hline, 'Color', 'black','Linewidth', Width)

        hold off;
    end;
end;


h = gca;
if LogXscale,
    [XTicks, XTickLabel]= CalculateTicksOld(LRx);
    set(gca, 'XTick', XTicks, 'XTickLabel', XTickLabel);
end;
if (~HistPlot)&LogYscale,
    [YTicks, YTickLabel]= CalculateTicks(LRy);
    set(gca, 'YTick', YTicks, 'YTickLabel', YTickLabel);
end;



function Xlim = FindPlotLimits(Xdata, OF)
%find plot limits by cutting out outliers
len = length(Xdata);
Xdata = sort(Xdata, 'ascend');
Xlim = Xdata([ceil(len*OF) floor(len*(1-OF))]);

function h = PlotColorCoded(XX, YY, Jxx, Jyy, Nbins)
Jxy = Jxx&Jyy;
%prepare colorcoding
[DensityMatrix, Bins] = hist3([XX(Jxy) YY(Jxy)], [Nbins Nbins]);
Color = interp2(Bins{2}, Bins{1}, DensityMatrix, YY(Jxy), XX(Jxy), '*linear');
colormap('jet');
h = scatter(XX(Jxy), YY(Jxy), 1.5, Color,'fill');
%h = fastscatter(XX(Jxy), YY(Jxy),Color, 'markersize',4);%AOY
grid off;

axis tight
% hold on;
% % now plot outliers
% Xlim = get(gca, 'Xlim');
% Xtemp = XX(~Jxx);
% Ytemp = YY(~Jxx);
% Jtemp = (Xtemp < Xlim(1));
% plot(Xlim(1)*ones(1, sum(Jtemp)), Ytemp(Jtemp), 'x');
% Jtemp = (Xtemp > Xlim(2));
% plot(Xlim(2)*ones(1, sum(Jtemp)), Ytemp(Jtemp), 'x');
% 
% Ylim = get(gca, 'Ylim');
% Xtemp = XX(~Jyy);
% Ytemp = YY(~Jyy);
% Jtemp = (Ytemp < Ylim(1));
% plot(Xtemp(Jtemp), Ylim(1)*ones(sum(Jtemp), 1), 'x');
% Jtemp = (Ytemp > Ylim(2));
% plot(Xtemp(Jtemp), Ylim(2)*ones(sum(Jtemp), 1), 'x');
% hold off;


function [Ticks, TickLabel] = CalculateTicksOld(LinearRange) % not used: to remove
a1 = -(9:-1:1)'*10.^(4:-1:1);
%a1 = sort(- (1:9)'*10.^(1:2)); 
a2 =(1:9)'*10.^(1:5);
%a = [sort(a1(:)); 0; a2(:)];
a = [-100000; a1(:); 0; a2(:)];

emptyStr = {''; ''; '';''; ''; '';''; ''};
TickLabel = ['-100000'; emptyStr;'-10000'; emptyStr; '-1000'; emptyStr; '-100'; emptyStr; {''}; '0'; {''}; emptyStr; '100'; emptyStr; '1,000'; emptyStr; '10,000'; emptyStr; '100,000'; emptyStr];

Ticks = asinh(a/LinearRange);

function [Ticks, TickLabel] = CalculateTicks(LinearRange)
% find the power of 10 closest to LinearRange (in log sense)
maxpwr = 6;
pwr = round(log10(LinearRange));
a1 = -10.^(maxpwr:-1:pwr);
a3 =10.^(pwr:maxpwr);
a = [ a1(:); 0; a3(:)];


% calculate minor ticks
a1 = -(9:-1:1)'*10.^(maxpwr:-1:pwr);
a2 = (-9:9)*10^(pwr-1);
a3 =(1:9)'*10.^(pwr:maxpwr);
%a = [sort(a1(:)); 0; a2(:)];
ma = [a1(:); a2(:); a3(:)];

TickLabel = cell(size(ma));

TickLabel(:) = {''};

J = ismember(ma, a);
TickLabel(J) = cellstr(num2str(ma(J)));
Ticks = asinh(ma/LinearRange);