function h = spyx(x,y,S,opt1,opt2)
% SPYX Visualize sparsity pattern on defined axes.
%
%    Syntax:
%       hline = spyx(x,y,S,opt1,opt2);
%
%    Input:
%       x    - A vector of length N containing the x-axis values. 
%       y    - A vector of length M containing the y-axis values.
%       S    - A 2D matrix of size MxN to be check for sparsity.
%       opt  - Two optional input defining the 'LineSpec' and MarkerSize.
%              Can be none, one, or both, and in any order (see SPY for  
%              details).
%
%    Output
%       hline - Handle of the lines drew by the SPY function (optional).
%
%    Description:
%       This program is equal as SPY function (plots marks on the non zero
%       elements of S, even on Nan's) but uses the user defined axes
%       instead of the normal rows and columns indexes. It uses the SPY
%       function. Note: sets the axes limits to those of x,y and sets
%       'YDir' to normal (increasing upwards).
%
%    Example:
%       figure
%       x = linspace(-10,10,10);
%       y = linspace(20,60,10);
%       S = round(round(rand(10)).*rand(10)*100)/10; 
%       S([14 15 24 25 48 68 100]) = NaN; num2str(S,' %4.1f')
%
%       % NORMAL SPY:
%       subplot(221)
%        spy(S) 
%        ylabel('y-axis reverted.')  
%        title('SPY(S)')
%
%       % NEW SPYX:
%       subplot(222)
%        spyx(x,y,S)
%        title('SPYX(x,y,S)  (compare the axes ticks)')
%        ylabel('y-axis set to normal.')  
%
%       % MATRIX S WITH IMAGESC:
%       subplot(223)
%        imagesc(x,y,flipud(S)), set(gca,'YDir','normal')
%        title('Not seen NaN''s of S...')
%        xlabel x, ylabel y, lim = axis;
%
%       % MATRIX S WITH IMAGESC AND NAN'S MARKS:
%       subplot(224)
%        hold on
%         imagesc(x,y,flipud(S)), 
%         h = spyx(x,y,isnan(S),'*g',20);
%        hold off
%        title('Oh! There you are.')
%        xlabel x, ylabel y,  axis(lim)
%        set(gca,'Layer','top','Box','on')
%    
%        % To draw the marks at high 15 instead of 0, for example, use:
%        set(h,'ZData',repmat(15,1,7)) % 7 NaN's 
%
%
%    See also SPY.

%   Copyright 2008 Carlos Adrian Vargas Aguilera
%   $Revision: 1.0$  $Date: 2008/03/12 15:30:00 $

%   Written by
%   M.S. Carlos Adrian Vargas Aguilera
%   Physical Oceanography PhD candidate
%   CICESE 
%   Mexico, 2008
%   nubeobscura@hotmail.com
%
%   Download from:
%   http://www.mathworks.com/matlabcentral/fileexchange/loadAuthor.do?objec
%   tType=author&objectId=1093874

%% Errors checking and defaults
% Check number of inputs:
Nargin = nargin;
if ~ismember(Nargin,[3:5])
 error('Spyxy:InputNumber','Incorrect number of inputs.')
end
% Compare input sizes:
[m,n] = size(S);
if (~isempty(x) && length(x)~=n) || (~isempty(y) && length(y)~=m)
 error('Spyxy:InputSize','Incorrect size of inputs.')
end
% Ignores empty options:
if Nargin==4 
 if isempty(opt1), Nargin = 3; end 
end
if Nargin==5 
 if isempty(opt1), Nargin = 4; end
 if isempty(opt2), Nargin = Nargin-1; end 
end

%% Call to SPY
if Nargin<4
 spy(S)
elseif Nargin<5
 spy(S,opt1)
else
 spy(S,opt1,opt2)
end

%% Get line handle:
hline = findobj(get(gca,'Children'),'Type','line');

%% Check for correct handle if there are more lines
if length(hline)>1
 % Get current lines specifications:
 markersize = get(hline,'MarkerSize'); 
 marker     = get(hline,'Marker');
 linestyle  = get(hline,'LineStyle');
 color      = get(hline,'Color');
 
 % Get SPY line defaults:
  units = get(gca,'units');
  set(gca,'units','points');
  pos = get(gca,'position');
 dmarkersize = max(4,min(14,round(6*min(pos(3:4))/max(m+1,n+1))));
  set(gca,'units',units);
 dmarker     = '.';
 dlinestyle  = 'none'; 
 dcolor = get(gca,'colororder'); dcolor = dcolor(1,:); 
 % Change them if options:
 if Nargin>3 
  if ischar(opt1)
   [dlinestyle,dcolor,dmarker] = colstyle(opt1);
   if Nargin==5
    dmarkersize = opt2;
   end
  else
   dmarkersize = opt1;
   if Nargin==5
    [dlinestyle,dcolor,dmarker] = colstyle(opt2);
   end
  end
 end
 
 % Compare both:
 ih = 1:length(hline);
 ih = ih(ismember(cell2mat(markersize),dmarkersize));
 if length(ih)>1
  ih = ih(ismember(marker(ih),dmarker));
  if length(ih)>1
   ih = ih(ismember(linestyle(ih),dlinestyle));
   if length(ih)>1
    ih = ih(ismember(cell2mat(color(ih)),dcolor,'rows'));
    if length(ih)>1
     warning('Spyxy:IndistinguishLines',...
      'Couldn''t distinguish SPYXY lines from others. Axes ignored.')
     if nargout
      h = hline(ih);
     end
     return
    end
   end
  end
 end
 hline = hline(ih);
end

%% Sets the new axes and limits:
if isempty(x)
 if isempty(y)
  return
 else
  y = y(end:-1:1);
  j = get(hline,'YData');
  set(hline,'Ydata',y(j))
  set(gca,'YDir','normal')
  set(gca,'Ylim',[min(y) max(y)])
  axis normal
 end
else
 if isempty(y)
  i = get(hline,'XData');
  set(hline,'Xdata',x(i))
  set(gca,'Xlim',[min(x) max(x)])
  axis normal
 else
  y = y(end:-1:1);
  i = get(hline,'XData');
  j = get(hline,'YData');
  set(hline,'Xdata',x(i))
  set(hline,'Ydata',y(j))
  set(gca,'YDir','normal')
  set(gca,'Xlim',[min(x) max(x)],'Ylim',[min(y) max(y)])
  axis normal
 end
end

%% Returns the handle:
if nargout
 h = hline;
end