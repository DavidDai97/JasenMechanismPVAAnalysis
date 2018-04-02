function [] = formatplot( hfig,hax,htitle,hxlabel,hylabel)
%FORMATPLOT Format plot to user specified default parameters
% User must supply a figure handle (figh), axis handle (axh),title handle
% (htitle), xlabel handle (hxlabel), ylabel handle(hylabel)

set([htitle, hxlabel, hylabel], ...
    'FontName'   , 'Arial'); % sets font type for title, x/ylabels

set([hxlabel, hylabel]  , ...
    'FontSize'   , 10          ); % sets font size for x/y labels

set( htitle                    , ...
    'FontSize'   , 12          , ... %sets font size for title 
    'FontWeight' , 'bold'      ); % sets font weight for title

%Axis Parameters
set(hax, ...
  'FontSize'   , 9 ,... %number font size
  'FontName'    ,'Arial',...% number font name
  'Box'         ,'off', ... %turns on/off right and upper plot edges
  'TickDir'     ,'out', ... % tick direction in/out
  'TickLength'  ,[.02 .02], ... %tick length
  'XMinorTick'  ,'on', ... % x Minor ticks on/off
  'YMinorTick'  ,'on', ... % y Minor ticks on/off
  'YGrid'       ,'off', ... % y grid on/off
  'XGrid'       ,'off',... % x grid on/off
  'XColor'      ,[.3 .3 .3], ...% xaxis color 
  'YColor'      ,[.3 .3 .3], ... % yaxis color
  'LineWidth'   ,1); %line width of axes in pts

 set(hfig,...
     'Color',[1 1 1],...
     'PaperPositionMode', 'auto'); % makes exported figure look like 
                                   % it does on screen
                                   
hdata=get(hax,'Children'); %gets handles to all data in axis

for i=1:length(hdata) % Iterate over all data types in plot
    
    if strcmp(hdata(i).Type,'scatter') % Finds scatter type data
        set(hdata(i),...
            'SizeData', 50,... %Sets Marker size
            'LineWidth', 1); % Sets line width 
    end
  
    if strcmp(hdata(i).Type,'line') % Finds line type data
        set(hdata(i),...
            'LineWidth', 1); % Sets line width 
    end
    
end

