function output = calc_profpars(data, graph, hrange, pars, setvalues)

    cm2ccm = 1e-6;
    
    isedens = 1;

    if isedens == 1        
        edens = data.edens * cm2ccm;
        err_edens = data.err_edens * cm2ccm;
    else        
        edens = dens2freq(data.edens, 2);
        err_edens = dens2freq(data.err_edens, 2);
    end
    
    height = data.height;
    time = data.time;
       
    global hvalues;
    global Nm hm;
    
    global hvalues_bottom;
    
    indh = find(height >=  hrange(1) & height <= hrange(2));
    if ~isempty(indh)
        edvalues = edens(indh)'; err_edvalues = err_edens(indh)';
        hvalues = height(indh)';
    else
        edvalues = edens'; err_edvalues = err_edens';
        hvalues = height';
    end

    hvals = hrange(1) : hrange(2);
        
    
    % Fitting to Chapman profile model
    %
    funcname = 'achapman_funct';
    jacname = 'achapman_jac';
    m0 = pars(1:3)';
    max_nloop = 400;
    ydata = edvalues;
    sigma = err_edvalues;
    
    [m, nloop] = lm_method(funcname, jacname, m0, max_nloop, ydata, sigma);
    
    Nm = m(1); hm = m(2);
    
    N_fit = achapman_feval(hvals, m);

    
    % Fitting to IRI F2-region bottomside profile
    %
    
    % Finding the minimum height of F2 bottomside region
    %
    fact_NmF2B = setvalues.FactNm; min_height = setvalues.MinHeight;
    
    min_NmF2B = fact_NmF2B * Nm;
    indh = find(height <= hm);
    indedb = find(edens(indh) >= min_NmF2B & height(indh) >= min_height);
    
    edvalues_bottom = edens(indh(indedb));
    err_edvalues_bottom = err_edens(indh(indedb));
    hvalues_bottom = height(indh(indedb));

    edvalues_bottom = cat(2, edvalues_bottom, Nm);
    err_edvalues_bottom = cat(2, err_edvalues_bottom, ...
        err_edvalues_bottom(numel(err_edvalues_bottom)));
    hvalues_bottom = cat(2, hvalues_bottom, hm)';
    
    myHVals = min(hvalues_bottom) : hm;        
   
    m0_bottom = pars(4:5)';
    ydata_bottom = edvalues_bottom';

    %options = optimset('Jacobian','on');

    [m_bottom, resnorm, residual, exitflag, output] = ....
        lsqcurvefit(@iri_f2bottom_functeval, m0_bottom, ...
        hvalues_bottom, ydata_bottom);

    N_fit_bottom = iri_f2bottom_functeval(m_bottom, myHVals);

    nloop_bottom = output.iterations;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Chi-square (bottomside)
    %
    % PENDING ...
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    output = struct('m', cat(1,m, m_bottom), 'nloop', ...
        [nloop nloop_bottom]);
    
    
    % Graphical output
    %
    
    if graph > 0
        
        npRows = 1; npCols = 2;
        
        close all;
        hf = figure(1); clf(hf, 'reset');
        hfPosition = [10 40 npCols*600 npRows*500];
        set(hf, 'Position', hfPosition);
        
        subplot(npRows, npCols, 1);
        
        errorbar(hvalues, edvalues, err_edvalues, 'Color', 'k', ...
            'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', setvalues.MarkerSize);
                
        line(hvals, N_fit, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        
        curr_YTick = get(gca, 'YTick');
        for j = 1 : numel(curr_YTick)
            tmpYTickLabel = sscanf(num2str(curr_YTick(j)), '%E');
            if j > 1
                myYTickLabel = cat(1, myYTickLabel, tmpYTickLabel);
            else
                myYTickLabel = tmpYTickLabel;
            end
        end
        
        set(gca, 'FontName', setvalues.FontName, 'FontSize', setvalues.FontSize, ...
            'XLim', hrange, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'YTickLabel', myYTickLabel);
           
        xlabel('Height (km)');         
        
        if isedens == 1            
            title(['F-region: ' datestr(data.time,0) ' ']); 
            ylabel('Ne(cm^-3)');            
        else            
            title(['JRO ' datestr(data.time,0) ' ']); 
            ylabel('MHz');            
        end
        
        curr_XLim = get(gca, 'XLim'); curr_YLim = get(gca, 'YLim');
        x0Text = curr_XLim(1) + 0.15 * abs(diff(curr_XLim));
        
        strText(1) = {'Initial guess: '};
        if isedens == 1
            strText(2) = {['Nm = ' num2str(m0(1), '%8.5E') ' cm^-^3']};
        else
            strText(2) = {['fm = ' num2str(m0(1), '%8.3f') ' MHz']};       
        end
        strText(3) = {['hm = ' num2str(m0(2), '%8.3f') ' km']};
        strText(4) = {['H = ' num2str(m0(3), '%8.3f') ' km']};

        strText(5) = {'Results: '};
        if isedens == 1
            strText(6) = {['Nm = ' num2str(m(1), '%8.5E') ' cm^-^3']};
        else
            strText(6) = {['fm = ' num2str(m(1), '%8.3f') ' MHz']};
        end
        strText(7) = {['hm = ' num2str(m(2), '%8.3f') ' km']};
        strText(8) = {['H = '  num2str(m(3), '%8.3f') ' km']};
        
        strText(9) = {['Iterations = ' num2str(nloop, '%i')]};

        yfactor = linspace(0.48, 0.18, numel(strText));
        

        for i = 1 : numel(strText)
            y0Text = curr_YLim(1) + yfactor(i) * abs(diff(curr_YLim));
            text(x0Text, y0Text, strText(i), 'FontName', setvalues.FontName, ...
                'FontSize', 0.65*setvalues.FontSize);
        end
        
        lgh = legend('ISR data','\alpha-Chapman profile','Location','Best');
        set(lgh, 'FontSize', 0.55*setvalues.FontSize)

        
        % Plot # 2
        %
        subplot(npRows, npCols, 2);

        errorbar(hvalues_bottom, edvalues_bottom, err_edvalues_bottom, ...
            'Color', 'k', 'LineStyle', 'none', 'Marker', 'o', ...
            'MarkerSize', setvalues.MarkerSize);
        
        line(myHVals, N_fit_bottom, 'Color', 'k', 'LineStyle', '-');

        curr_YTick = get(gca, 'YTick');
        for j = 1 : numel(curr_YTick)
            tmpYTickLabel = sscanf(num2str(curr_YTick(j)), '%E');
            if j > 1
                myYTickLabel = cat(1, myYTickLabel, tmpYTickLabel);
            else
                myYTickLabel = tmpYTickLabel;
            end
        end
        
        set(gca, 'FontName', setvalues.FontName, 'FontSize', setvalues.FontSize, ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', myYTickLabel);
        
        xlabel('Height (km)'); ylabel('Ne(cm^-3)');
        title('F-region bottomside'); 

        curr_XLim = get(gca, 'XLim'); curr_YLim = get(gca, 'YLim');        
        
        x0Text = curr_XLim(1) + 0.65 * abs(diff(curr_XLim));
        
        strTextB(1) = {'Initial guess: '};
        strTextB(2) = {['B0 = ' num2str(m0_bottom(1), '%8.3f') ' km']};
        strTextB(3) = {['B1 = ' num2str(m0_bottom(2), '%8.3f')]};       

        strTextB(4) = {'Results: '};
        strTextB(5) = {['B0 = ' num2str(m_bottom(1), '%8.3f') ' km']};
        strTextB(6) = {['B1 = ' num2str(m_bottom(2), '%8.3f')]};
        
        strTextB(7) = {['Iterations = ' num2str(nloop_bottom, '%i')]};

        yfactor = linspace(0.5, 0.2, numel(strText));                
        
        for i = 1 : numel(strTextB)
            y0Text = curr_YLim(1) + yfactor(i) * abs(diff(curr_YLim));
            text(x0Text, y0Text, strTextB(i), 'FontName', setvalues.FontName, ...
                'FontSize', 0.65*setvalues.FontSize);
        end
        
        lgb = legend('ISR data','IRI Bottomside profile','Location','Best');
        set(lgb, 'FontSize', 0.55*setvalues.FontSize)
        %

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Chi-square plot
        %
        % PENDING ...
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [yy, mm, dd, hh, mi, ss] = datevec(time);
        fname = ['fittprof_' num2str(yy,'%04i') num2str(mm,'%02i') ...
            num2str(dd,'%02i') num2str(hh,'%02i') num2str(mi,'%02i') ...
            num2str(round(ss),'%02i')];
        figname = [setvalues.GPath fname];
        save_figure(hf, graph, 300, 'portrait', figname)
        
    end
