fhan = gcf;
axhan = get(fhan, 'Children');

fhan.Units = 'centimeters';
fhan.PaperUnits = 'centimeters';
wid = 8.5;
textsize = 8;
ticksize = 8;
MSZ      = 4; % Marker size

fhan.InnerPosition = [0,0,wid,wid*0.6321];
% fhan.PaperPosition = [0,0,wid,wid*0.6321];


theinterp = 'latex';

if length(axhan)>1
    for idx = 1 : length(axhan)
        try
            axhan(idx).XAxis.FontSize = textsize;
            axhan(idx).XAxis.Label.Interpreter = theinterp;
            axhan(idx).XAxis.Label.FontSize    = textsize;
            % axhan(idx).XAxis.Label.String      a= '';    
            axhan(idx).YAxis.Label.Interpreter = theinterp;
            axhan(idx).YAxis.Label.FontSize    = textsize;
            % axhan(idx).YAxis.Label.String      = '';        
            axhan(idx).XAxis.TickLabelInterpreter = theinterp;
            axhan(idx).YAxis.TickLabelInterpreter = theinterp;   
            for jdx = 1:length(axhan(idx).Children)
                axhan(idx).Children(jdx).MarkerSize = MSZ;
            end
        end
    end
else % There can be plots with yyaxis right...
    axhan.XAxis.FontSize = textsize;
    axhan.XAxis.Label.Interpreter = theinterp;
    axhan.XAxis.Label.FontSize    = textsize;
    axhan.XAxis.TickLabelInterpreter = theinterp;
        
    try
        axhan.YAxis.Label.Interpreter = theinterp;
        axhan.YAxis.Label.FontSize    = textsize;
        axhan.YAxis.TickLabelInterpreter = theinterp;
        for jdx = 1:length(axhan(idx).Children)
           axhan(idx).Children(jdx).MarkerSize = MSZ;
        end
    catch
        axhan.YAxis(1).Label.Interpreter = theinterp;
        axhan.YAxis(2).Label.Interpreter = theinterp;
        axhan.YAxis(1).Label.FontSize    = textsize;
        axhan.YAxis(2).Label.FontSize    = textsize;
        axhan.YAxis(1).TickLabelInterpreter = theinterp;
        axhan.YAxis(2).TickLabelInterpreter = theinterp;
        for jdx = 1:length(axhan(idx).Children)
            axhan(idx).Children(jdx).MarkerSize = MSZ;
        end
    end
end

