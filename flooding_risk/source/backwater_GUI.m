function backwater_GUI()

    % interactive backwater module written by Andrew J. Moodie
    % see module and website for more information
    % classroom module for this model can be found at 
    % http://www.coastalsustainability.rice.edu/outreach/
    % the model setup below is parameterized to the Lower Mississippi River
    % as established by Nittrouer et al., Spatial and temporal
    % trends...GSAB, 2012
    % S0 = 7e-5, B0 = 1100 m, Qbf = 35000
    
    L = 1600e3; % length of domain
    nx = 400; % number of nodes
    dx = L/nx; % width of cells
    x = 0:dx:L; % define x-coordinates
    
    start = 43; % pin-point to start eta from
    S0 = 7e-5; % bed slope
    eta = linspace(start, start - S0*(nx*dx), nx+1); % channel bed
    
    mou = 0.75; % fraction of x channelized (i.e. mouth position)
    thet = 2; % plume spreading angle
    
    B0 = 1100; % basic channel width
    [B] = set_B(B0, mou, thet, nx, dx); % channel width
    [S] = get_slope(eta, nx, dx); % bed slope at each node
    H0 = 0; % 
    Cf = 0.0047;
    Qwinit = 10000;
    Qw = Qwinit;
    Qwbf = 35000;
    Qwmax = 60000;
    Qwmin = 5000;
    
    [H] = get_backwater_dBdx(eta, S, B, H0, Cf, Qw, nx, dx);
    [Xs] = find_backwaterregion(H, dx);
    [zed] = 0.5 + get_backwater_dBdx(eta, S, B, H0, Cf, Qwbf, nx, dx);
    
    [nitt] = load_nittrouer2011();
    
    % create and then hide the GUI as it is being initialized
    f = figure('Visible','off','Position',[360,500,800,600], 'Units', 'normalized');
    % create the axes
    hand.ax = axes('Units','Pixels','Position',[60,300,700,250], 'Units', 'normalized');
                hold on
                ylim([-50 100])
    hand.Qw = uicontrol('Style', 'slider', ...
                'Min', 5000, 'Max', 60000, 'Value', Qw, 'Position', [60, 200, 320, 25], ...
                'Units', 'normalized', 'Callback', @updateQw);
    uicontrol('Parent', f, 'Style', 'text', 'Position', [20, 195, 40, 25], ...
                'String', format_number(Qwmin), 'BackgroundColor', f.Color, 'Units', 'normalized');
    uicontrol('Parent', f, 'Style', 'text', 'Position', [380, 195, 60, 25], ...
                'String', format_number(Qwmax), 'BackgroundColor', f.Color, 'Units', 'normalized');
    uicontrol('Parent', f, 'Style', 'text', 'Position', [60, 170, 320, 25], ...
                'String', 'discharge (m^3/s)', 'BackgroundColor', f.Color, 'Units', 'normalized');
    hand.QwValue = uicontrol('Parent', f, 'Style', 'text', 'Position', [540, 490, 200, 25], ...
                'String', ['Qw = ' format_number(Qw)], 'BackgroundColor', [1 1 1], 'FontSize', 14, 'Units',  'normalized');
    
    nittpanel = uipanel('Visible','off', ...
        'Title', 'Nittrouer et al., 2011 data', 'Position', [0.675 0.24 0.28 0.17]);
    nittthal = uicontrol(nittpanel, 'Style', 'checkbox', ...
                'String', 'show thalweg', ...
                'Value', 0, 'Position', [25 50 150 30], 'units', 'normalized', ...
                'Callback', @showThal);
    nittwater = uicontrol(nittpanel, 'Style', 'checkbox', ...
                'String', 'show water lines', ...
                'Value', 0, 'Position', [25 15 150 30], 'units', 'normalized', ...
                'Callback', @showWater);
    nittpanel.Visible = 'on';
            
    uicontrol('Parent', f, 'Style', 'pushbutton', 'Position', [585, 95, 100, 25], ...
                'String', 'Reset', 'units', 'normalized', 'Callback', @reset);
    uicontrol('Parent', f, 'Style', 'pushbutton', 'Position', [585, 55, 100, 25], ...
                'String', 'About', 'units', 'normalized', 'Callback', @moduleinformation);
            
    % create data table for information
    RKs = [0 165 368 425 505];
    RKidxs = (nx*mou) - round(RKs*1000/dx);
    tab = uitable(f);
    tab.Data = {'Head of Passes (RK 0)', '0', '0', false; 'New Orleans (RK 165)', '0', '0', false; 'Baton Rouge (RK 368)', '0', '0', false; ...
        'St. Francisville (RK 425)', '0', '0', false; 'Old River Diversion (RK 505)', '0', '0', false};  
    tab.ColumnName = {'location', 'flow depth (m)', 'stage (m)', 'over levee?'};
    tab.Position = [30 30 480 117]; % [450 20 300 210];
    tab.Units = 'normalized';
    tab.ColumnWidth = {185 100 80 80};
    tab.Data(1:end, 2) = (  sprintfc( '%10.1f', (H(RKidxs))' )  ); % insert proper depths
    tab.Data(1:end, 3) = (  sprintfc( '%10.1f', (eta(RKidxs)+H(RKidxs))' )  ); % insert proper stage
    
    % make the plot into the axes
    hand.etaLine = plot(x/1000, eta, 'k-', 'LineWidth', 1.2);
    hand.zedLine = plot(x(1:nx*mou)/1000, eta(1:nx*mou)+zed(1:nx*mou), 'k--', 'LineWidth', 1.2);
    hand.Hline = plot(x/1000, eta+H, 'b-', 'LineWidth', 1.2);
    hand.mouth = plot(repmat(L*mou/1000, 1, 2), ylim, 'k--');
    hand.bwbracket = plot([Xs(1) Xs(1) Xs(2) Xs(2)]/1000, [36 40 40 36], 'k-');
    hand.nittBed = plot(x/1000, nitt.bed, 'Marker', '.', 'LineStyle', 'none', ...
        'Color', [0.5 0.5 0.5], 'Visible', 'off');
    hand.nittWater = plot(repmat((x/1000)', 1, 3), nitt.surf', 'LineStyle', '--', ...
        'Color', [0.5 0.5 0.5], 'Visible', 'off');
    hand.bwvalue = text(((Xs(2)-Xs(1))/4 + Xs(1))/1000, 52, ...
        {'backwater from', ['RK ' num2str(L*mou/1000-round(Xs(1)/1000)), ' to ', num2str(L*mou/1000-round(Xs(2)/1000))]});
    plot(repmat(L/1000*mou-RKs(2:end), 2, 1), repmat(ylim', 1, length(RKs)-1), 'k:')
    xlim([L/1000*0.25 L/1000-(L/1000*0.125)])
    set(hand.ax, 'xTickLabels', cellfun(@num2str, num2cell(abs((cellfun(@str2num, (get(gca, 'XTickLabels')))) - (L/1000*mou))), 'UniformOutput', false))
    box on; set(hand.ax, 'FontSize', 10, 'LineWidth', 1.5);
    
    % make labels on plot
    text(L/1000*mou+10, 25, '< Head of Passes')
    text(L/1000*mou-RKs(2)-10, -40, '\^New Orleans', 'BackgroundColor', [1 1 1 0.85])
    text(L/1000*mou-RKs(3)-7, -34, '\^Baton Rouge', 'BackgroundColor', [1 1 1 0.85])
    text(L/1000*mou-570, -23, 'St. Francisville\^', 'BackgroundColor', [1 1 1 0.85]) % need to fix to calc for place
    text(L/1000*mou-690, -12, 'Old River Diversion\^', 'BackgroundColor', [1 1 1 0.85]) % need to fix to calc for place
    xlabel('distance from Head of Passes (km)')
    ylabel('elevation (m)')
    
    % Make the UI visible.
    f.Visible = 'on';

    
    % define callback functions
    function updateQw(source, ~)
        Qw = source.Value;
        [H] = get_backwater_dBdx(eta, S, B, H0, Cf, Qw, nx, dx);
        [Xs] = find_backwaterregion(H, dx);
        
        set(hand.Hline, 'YData', eta+H);
        set(hand.QwValue, 'String', ['Qw = ', format_number(round(Qw))]);
        set(hand.bwbracket, 'XData', [Xs(1) Xs(1) Xs(2) Xs(2)]/1000);
        hand.bwvalue.String(2) = {['RK ' num2str(L*mou/1000-round(Xs(1)/1000)), ' to ', num2str(L*mou/1000-round(Xs(2)/1000))]};
        hand.bwvalue.Position(1) = ((Xs(2)-Xs(1))/4 + Xs(1))/1000;
        updateTable()
    end
    
    function updateTable()
        tab.Data(1:end, 2) = (  sprintfc( '%10.1f', (H(RKidxs))' )  ); % new depths
        tab.Data(1:end, 3) = (  sprintfc( '%10.1f', (eta(RKidxs)+H(RKidxs))' )  ); % new stage
        over = ((eta+H) > eta+zed);
        tab.Data(1:end, 4) = arrayfun((@(x) {x}), over(RKidxs)');
    end

    function showThal(~, ~)
        state = get(nittthal, 'Value');
        if state == 1
            set(hand.nittBed,'Visible','on');
        elseif state == 0
            set(hand.nittBed,'Visible','off');
        end
    end

    function showWater(~, ~)
        state = get(nittwater, 'Value');
        if state == 1
            set(hand.nittWater,'Visible', 'on');
        elseif state == 0
            set(hand.nittWater,'Visible', 'off');
        end
    end

    function reset(~, ~)
        temp.Value = Qwinit;
        updateQw(temp);
        hand.Qw.Value = Qwinit;
        nittthal.Value = 0;
        set(hand.nittBed,'Visible','off')
        nittwater.Value = 0;
        set(hand.nittWater,'Visible','off')
    end

    function moduleinformation(~, ~)
        uiwait(msgbox({'Model parameterization follows Nittrouer et al., 2012 for Mississippi River. The parameters can be viewed in the source code at XXXXXXXgithub.', ...
            ' ', ...
            'The module was created by Andrew J. Moodie as part of an National Science Foundation funded research project assessing the sustainability of anthropogenically influenced deltas.', ...
            'The research is supported by Grant No. XXXXXXXX, and an NSF Graduate Research Fellowship to A.J.M. under Grant No. XXXXXXX.', ...
            ' ', ...
            'For more information, please visit http://www.coastalsustainability.rice.edu/outreach/'}, 'About this module', 'modal')); 
    end
end

function [H] = get_backwater_dBdx(eta, S, B, H0, Cf, Qw, nx, dx)
    % backwater formulated for changing width
    H = NaN(1,nx+1); % preallocate depth 
    H(nx+1) = abs(H0 - eta(nx+1)); % water depth at downstream boundary
    for i = nx:-1:1
        % predictor step: computation of a first estimation of the water depth Hp
        [Frsqp] = get_froude(Qw, H(i+1), B(i)); % calculate froude from conditions at i+1
        [dBdx] = (B(2:nx+1) - B(1:nx)) ./ dx;
        [dHdxp] = get_dHdx_dBdx(S(i+1), Cf, Frsqp, H(i+1), B(i), dBdx(i)); % get dHdx width changing form
        Hp = H(i+1) - dHdxp * dx; % solve for H prediction
        % corrector step: computation of H
        [Frsqc] = get_froude(Qw, Hp, B(i)); % calculate froude at i with prediction depth
        [dHdxc] = get_dHdx_dBdx(S(i), Cf, Frsqc, Hp, B(i), dBdx(i)); % doaa
        % convolution of prediction and correction, trapezoidal rule
        H(i) = H(i+1) - ( (0.5) * (dHdxp + dHdxc) * dx );
    end
    
    function [dHdx] = get_dHdx_dBdx(S_loc, Cf, Frsq, H_loc, B_loc, dBdx)
        % formulation to get dHdX for a changing width backwater formulation
        dHdx = ((S_loc-(Cf*Frsq))/(1-Frsq)) + (Frsq/(1-Frsq)) * (H_loc/B_loc) * (dBdx);
    end
  
    function [Frsq] = get_froude(Qw, H, B)
        g = 9.81; % gravitational acceleration constant
        Frsq = ( Qw^2 / (g * B^2* H^3) );
    end
end

function [Xs] = find_backwaterregion(H, dx)
    dHdx = (H(2:end)-H(1:end-1))/dx;
    dHdxdx = (dHdx(2:end)-dHdx(1:end-1))/dx;
    threshold = dHdxdx*1e9 < 0.080;
    changepts = abs(threshold(2:end) - threshold(1:end-1));
    Xs = find(changepts, 2) * dx;
end

function [B] = set_B(B0, mou, thet, nx, dx)
    B = zeros(1, nx+1);
    chanLen = mou*nx;
    B(1:chanLen) = B0;
    B(chanLen+1:end) = B0 + 2*((1:(nx+1-chanLen)) * dx * tand(thet));
end

function [S] = get_slope(eta, nx, dx)
    % return slope of input bed (eta)
    S = zeros(1,nx+1);
    S(1) = (eta(1) - eta(2)) / dx;
    S(2:nx) = (eta(1:nx-1) - eta(3:nx+1)) ./ (2*dx);
    S(nx + 1) = (eta(nx) - eta(nx + 1)) / dx;
end

function [string] = format_number(integer)
    string = sprintf(',%c%c%c',fliplr(num2str(integer)));
    string = fliplr(string(2:end));
end

function [nitt] = load_nittrouer2011()
    nitt.bed = linspace(43, 43 - 7e-5*(400*4000), 400+1)+rand(1, 401)*10-6;
    nitt.surf = [(linspace(43, 43 - 7e-5*(400*4000), 400+1)+rand(1, 401)*2+20); ...
        (linspace(43, 43 - 7e-5*(400*4000), 400+1)+rand(1, 401)*2+32); ...
        (linspace(43, 43 - 7e-5*(400*4000), 400+1)+rand(1, 401)*2+35)];
    
end
