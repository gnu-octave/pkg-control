## Copyright (C) 2012 - 2020 Torsten Lilge
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File}  {} rlocusx (@var{sys})
## @deftypefnx {Function File} {} rlocusx (@var{sys}, @var{increment}, @var{min_k}, @var{max_k})
## Interactive root locus plot of the specified SISO system @var{SYS}.
##
## This functions
## directly calls rlocus() from the control package
## which must be installed and loaded.
## In contrast to rlocus(), mouse clicks on the root locus display
## the related gain and all other poles resulting from this gain
## together with damping and frequency of conjugate complex pole pairs.@*
## All possible interaction by mouse clicks or keys are:
##
## @table @asis
## @item Left click: Gain, damping and Frequency
##                Displays related gain and all resulting
##                closed loop poles together with damping
##                and frequency
## @item @kbd{s}: Step response
##                Simulates the step response for the gain of
##                of the most recently selected pole locations
## @item @kbd{i}: Impulse response
##                Simulates the impulse response for the most
##                recently selected gain
## @item @kbd{b}: Bode plot
##                Provides the open loop bode plot for the most
##                recently selected gain
## @item @kbd{m}: Stability margins
##                Provides the open loop bode plot with stability
##                margins for the most recently selected gain
## @item @kbd{a}: All plots
##                Provide sall four aforementioned plots
## @item @kbd{c}: Clear
##                Removes all closed loop pole markers and annotations
## @item @kbd{d}: Delete
##                Removes all open figures with simulation and
##                bode plots
## @item @kbd{x}: Exit
##                Exits the interactive mode and re-activates
##                the octave prompt
## @end table
##
## There are no output parameters.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model.  Must be a single-input and single-output (SISO) system.
## @item increment
## The increment used in computing gain values.
## @item min_k
## Minimum value of @var{k}.
## @item max_k
## Maximum value of @var{k}.
## @end table
##
## @strong{Outputs}
##
## Plots the interactive root locus to the screen. @*
## Unlike rlocus(), this function does not have any output parameters.
## For output parameters please directly use rlocus().
##
## @seealso{rlocus}
##
## @end deftypefn

## Author: Torsten Lilge
## Date: April 2020
## Version: 0.1

function rlocusx(sys,varargin)

  global fig_h

  % Check the input parameters
  if ( nargin < 1 || nargin >4 )
    print_usage();        % wrong number, so print usage
  end;

  if (! isa (sys, "lti") || ! issiso (sys))
    error ('rlocusx: first argument must be a SISO LTI model');
  endif

  % Call rlocus depending on number of arguments
  switch (nargin)
  case 1
    [rldata, k_break, rlpol, gvec, real_ax_pts] = rlocus(sys);
  case 2                  % compatibility with matlab
    kk = varargin{1};     % get the desired gains
    len_kk = size (kk);
    if (max (len_kk) < 2) || (min (len_kk) > 1)
      error ('rlocusx: second parameter has to be a one dimensional list with at least two elements');
    endif
    [rldata, k_break, rlpol, gvec, real_ax_pts] = rlocus(sys,(kk(end)-kk(1))/(max(size(kk))-1),kk(1),kk(end));
  case 3
    print_usage();        % wrong number, so print usage
  case 4
    kinc=varargin{1};
    kstart=varargin{2};
    kend=varargin{3};
    [rldata, k_break, rlpol, gvec, real_ax_pts] = rlocus(sys,kinc,kstart,kend);
  endswitch;



  % Remove additional gain values (why do they exists sometimes?)
  while ( size(gvec,2) > size(rlpol,2) )
    gvec(end) = [];
  end;

  % Remove double entries (rlocus sometimes seems to "stuck" at the same gain for
  % several calculated points; duplicate entries breaks later use of interp1
  len_gvec = size(gvec,2);
  j = 2;
  while (j <= len_gvec)
    if (gvec(j) == gvec(j-1))   % if the same as the one before
      gvec(j) = [];             % delete this entry in the gains
      rlpol(:,j) = [];          % and delete the related closed loop poles
      len_gvec--;
    else
      j++;
    end
  end

  % Following calculations are based on open loop data
  % in "rlocus canonical form"
  [z,p,V,tsam]=zpkdata(sys,'v');
   if ( V < 0 )
    V = abs(V);           % take the absolute and move the minus into the feedback
    negfb = 0;            % thus, we get pos. feddback
  else
    negfb = 1;            % neg. feedback is assumed by default
  end
  n = length(p);          % n: number of open loop poles
  m = length(z);          % m: number of open loop zeros
  nm= n-m;                % we will need the difference more than once

  % Scaling of the plot: get some min/max values
  xlim('manual');               % manual scaling x-axis
  ylim('manual');               % manual scaling y-axis
  nrl = numel(rlpol);
  im = reshape(imag(rlpol),[1,nrl]);
  re = reshape(real(rlpol),[1,nrl]);
  xmin = min(re);
  xmax = max(re);
  ymin = min(im);
  ymax = max(im);

  % Intersection of asymptotes
  if ( nm > 0 )
    sigmaw = (sum(p)-sum(z))/nm;
    xmin = min(xmin,sigmaw);
    xmax = max(xmax,sigmaw);
  end

  % Scale while making sure that all relevant information is included
  for j=1:m
    xmin = min(xmin,real(z(j)));
    xmax = max(xmax,real(z(j)));
    ymin = min(ymin,imag(z(j)));
    ymax = max(ymax,imag(z(j)));
  end;
  dx = (xmax-xmin)/20;
  dy = (ymax-ymin)/20;
  if ( dx < 0.0001 )
    dx = dy;
  end;
  if ( dy < 0.0001 )
    dy = dx;
  end;
  xmin = xmin-dx;
  xmax = xmax+dx;
  ymin = ymin-dy;
  ymax = ymax+dy;

  % Some constants for colors and markers
  col_zp  = [0.85000   0.32500   0.09800];  % color of open loop zeros and poles
  col_clp = [0.10000   0.10000   0.10000];  % color of closed loop poles
  col_rl1 = [0.00000   0.44700   0.74100];  % color of first branch of rlocus
  col_rl2 = [0.30100   0.74500   0.93300];  % color of last branch of rlocus
  col_y   = [0.59400   0.18400   0.55600];  % color closed loop output y and simulated poles
  col_r   = [0.92900   0.69400   0.12500];  % color closed loop reference
  lw      = 2;              % general line width
  lw_mark = 2;              % line width of markers
  lw_asym = 1;              % line width of asymptotes
  lw_stab = 1;              % line width for region of stability

  % Marker size depending on engine
  if (strcmp(graphics_toolkit(),'gnuplot'))
    ms = 6;
  else
    ms = 12;
  end

  % define a cell array for storing the figures with the
  % Title and lables
  clf (gcf ());
  rlocus_fig = gcf ();                  % remember for resetting focus to this figure

  kmin = sprintf("%.3f",gvec(1));       % Interval of gain
  kmax = sprintf("%.3f",gvec(end));

  title({['Root locus (K = ',kmin,' ... ',kmax,')'];...
          'click: gain, s/i: step/imp., b/m: bode/margin, a: all, c: clear, d: del fig., x: end'});

  hold on;       	% from here on, do not delete what is drawn so far

  % Draw the stability region and put axes labels depending
  % on the domain (continuous ot discrete)
  if ( tsam > 0 ) % discrete time
    xlabel('Real axis');            % z-plane: just real and imag axis
    ylabel('Imag. axis');
    t = 0:0.05:6.3;
    plot(cos(t),sin(t),'-k');      % draw the unit circle (stability)
  else
    xlabel('Real axis: \sigma');    % s-plane: sigma + j omega
    ylabel('Imag. axis: \omega');
    plot([0,0],[ymin,ymax],'-k');  % draw imag axis (stability)
  end

  % Disable legend autoupdate, '' avoids a dummy entry 'data 1'
  legend ('','autoupdate','off');

  % Draw the open poles and zeros
  plot(real(p),imag(p),'x;open loop poles;','markersize',ms,'linewidth',lw,'color',col_zp);
  plot(real(z),imag(z),'o;open loop zeros;','markersize',ms,'linewidth',lw,'color',col_zp);

  % Draw root locii
  for j = 1:n
    if ( j == 1 )           % legend only once
      form = '-;locus;';
    else
      form = '-';
    end
    cj = (j-1)/(n-1); % each branch in a slightly different color
    col = col_rl2*cj + col_rl1*(1-cj);
    plot(real(rlpol(j,:)),imag(rlpol(j,:)),form,'linewidth',lw,'color',col);
  end
  grid on

  % Draw the asymptotes
  if ( nm > 0 )                     % only if there are asymptotes
    aslen  = max(abs(sigmaw - [ xmin xmax ymin ymax ]));  % the max visible length
    for j = 1:nm                    % loop for all asymptotes
      ang = (2*(j-1)+negfb)*pi/nm;  % the angle of the asymptotes
      endp = [sigmaw,0] + aslen*[cos(ang),sin(ang)]; % the endpoint
      if ( j == 1 )                 % legend only once
        form = '--;asymptotes;';
      else
        form = '--';
      end
      plot( [sigmaw, endp(1)], [0 endp(2)], form, 'color',[0.4,0.4,0.8] );  % finally plot it
    end
  end

  % Nice shaping and legend
  xlim([xmin,xmax]);
  ylim([ymin,ymax]);
  legend ('location','northwest');

  % Possible input values for ginput
  m_left =   1;     % left mouse button
  b_all  =  97;     % 'a'
  b_clr  =  99;     % 'c'
  b_del  = 100;     % 'd'
  b_exit = 120;     % 'x'
  b_step = 115;     % 's'
  b_imp  = 105;     % 'i'
  b_bode =  98;     % 'b'
  b_marg = 109;     % 'm'
  % Data of plots that can could be required for selected closed loop poles
  fig_h = [];             % Keep track of the figure handles that are used
                          % for simulations and bode plots for a specific K,
                          % rows: [K, fig step, fig imp, fog bode, fig margin]
  % Order of the plots, the first entry is a dummy for directly matching fig_h
  fig_b = [-1, b_step, b_imp, b_bode, b_marg];
  fig_n = {'','Step response (closed loop)','Imp. response (closed loop)',...
              'Bode plot (open loop)','Margin plot (open loop)'};

  % Some initializations
  button = 1;             % current button, not b_exit
  first = 1;              % counter for closed loop poles
  handles = [];           % list of handels for closed loop poles and text
  handle_sim_poles = 0;   % handles of closed loop poles with extra plots

  % As long as exit b_exit was not used
  while (button != b_exit)       % loop over all ginput values until 'x'

  [x,y,button] = ginput(1);      % wait for mouse/key event

  if length (button) == 0
    % No button -> figure was closed: Reset fig_h for not accessing
    % invalid handles in the close callbacks of remaining figures
    fig_h = [];
    return;
  endif

  if button == b_all             % plot all diagrams
    btns = fig_b;
  else
    btns = [button];
  endif

  for b = 1:length(btns)        % for all buttons (e.g., selected via "all")

    but = btns(b);              % the current button to handle here

    switch but                  % which mousen key or key?

    % left mouse key: markers and info for related root locus
    case ( m_left )
      s = x+y*i;                    % actual position in s-plane

      % rule 12 -> K_WOK, then divide by V (gain in zpk-form) gives real K
      K_tmp = prod(abs(s-p))/prod(abs(s-z))/V;

      % check if mouse click was really near to these locations
      distance = xmax - xmin;
      for j=1:n                     % for all poles
        poles(j) = interp1(gvec',rlpol'(:,j),K_tmp,'extrap');  % look up closed loop poles at our gain
        if abs (s-poles(j)) < distance
          distance = abs (s-poles(j));
        endif
      endfor

      % Plot all poles if we are "quite near" to the poles that belong to
      % the gain that was calculated from the click position
      if distance < (xmax - xmin)/25

        % Store gain
        K = K_tmp;

        % Plot closed loop together with related gain
        for j=1:n                     % for all poles
          x = real(poles(j));         % get real part
          y = imag(poles(j));         % and imaginary part
          if ( first )                % legend only once
            form = 'x;closed loop poles;';
            first = 0;
          else
            form = 'x';
          end

          % Plot the pole in the s plane
          h1 = plot (x,y,form, "markersize",ms*2/3,'linewidth',lw_mark,'color',col_clp);

          % Put the text with related parameters
          if ( y == 0*i )                 % if on real axis rotate text with gain
            str = sprintf("  K=%.2f",K);  % make a string from our gain
            h2 = text (x,y,str,'color',col_clp,'fontsize',10,'rotation',90);
          else                            % if not on real axis, text normal
            if ( tsam > 0 )               % calc related cont. pole if we
              s = log(x+y*i)/sys.tsam;    %   are in discrete time
              xs= real(s);
              ys= imag(s);
            else
             xs= x;
             ys= y;
            endif
            w0 = abs(xs+ys*i);            % get omega_0
            D  = -xs/w0;                  % and damping
            str = sprintf("  K=%.2f, D=%.2f, w_0=%.2f",K,D,w0);
            h2 = text(x,y,str,'color',col_clp,'fontsize',10); % put the text
          endif

          handles = horzcat(handles,[h2 h1]);  % store marker & text handles for
                                              % being able to remove them later
        endfor

      endif   % if near enough to the closed loop poles

    % 's', 'i', 'b', 'm': simulate closed loop or/and plot bode
    % with most recently selected poles
    case ({b_step, b_imp, b_bode, b_marg})

      if (! isempty (handles)) && (length (handles) >= 2*n)
        % handles for text and poles are not empty and contains 2n entries
        % (n poles and n texts), thus simulate and give related poles a
        % different color

        % Look for current K in the array of used figures
        if length (fig_h) == 0
          K_exists = 0;
          K_idx = 0;
        else
          [K_exists, K_idx] = ismember (K, fig_h(:,1));
        endif
        [b_exists, b_idx] = ismember (but, fig_b);

        if ! (K_exists && ishandle (fig_h(K_idx, b_idx)))
          % The desired figure does not yet exist

          if ! K_exists
            % Current K does not exist yet
            fig_h(end+1,1) = K; % new row
            K_idx = size (fig_h,1);
            % Initialize with invalid handles (fig_b is one element too long)
            fig_h(K_idx,2:2+length (fig_b)-2) = -1*ones (1,length (fig_b)-1);
            for j = 2:size(fig_h,2)

            endfor
            % update color of the related poles
            handle_sim_poles = handles(end-2*n+1:end);
            for j = 1:length (handle_sim_poles)
              set (handle_sim_poles(j), 'color', col_y);
            endfor
          endif

          % Create figure for desired plot
          fig_h(K_idx, b_idx) = ...
              figure ('DeleteFcn', {@__sim_fig_close_callback__, handle_sim_poles, col_clp},...
                      'name',['K = ',num2str(K),': ',fig_n{b_idx}],...
                      'numbertitle','off');

          % Do the desired plots
          switch (but)

            case {b_step, b_imp}
              % Simulate and plot the closed loop response

              grid on;
              hold on;

              if tsam > 0
                xlabel ('k');
              else
                xlabel ('t');
              endif

              closed_loop = feedback (K*sys, 1);
              if (but == b_step)
                [y,t] = step (closed_loop);
                plot ([t(1),t(end)],[1 1],'linewidth',lw,'color',col_r);
                plot (t,y,'linewidth',lw,'color',col_y);
                ylabel ('closed loop output and reference');
                title ({['Closed loop step response y for K = ',num2str(K)] });
                legend ('reference y_r','output y');
              elseif (but == b_imp)
                [y,t] = impulse (closed_loop);
                plot (t,y,'linewidth',lw,'color',col_y);
                ylabel ('impulse responce output and reference');
                title ({['Closed loop impulse response y for K = ',num2str(K)] });
                legend ('output y');
              endif

              figure (rlocus_fig)   % reset focus to root locus for ginput ()

            case {b_bode}
              % Bode plot of open loop
              bode (K*sys)
              figure (rlocus_fig)   % reset focus to root locus for ginput ()

            case {b_marg}
              % Bode plot of open loop
              margin (K*sys)
              figure (rlocus_fig)   % reset focus to root locus for ginput ()

          endswitch

        endif

      endif

    % Delete all closed loop markers
    case b_clr

      for j = 1:numel(handles)
        delete(handles(j));
      endfor
      handles=[];

    % Delete all existing figures related to closed loops
    case b_del

      % Collect extra figures that have to be closed from fig_h
      % and close them afterwards since closing will change fig_h
      % in the close callback
      figs_to_close = [];
      for j = 1:size (fig_h,1)
        for jj = 2:size (fig_h,2)
          if isfigure (fig_h(j,jj))
            figs_to_close(end+1) = fig_h(j,jj);
          endif
        endfor
     endfor
      % Close the figures and clean up
      for j = 1:length (figs_to_close)
        close (figs_to_close(j));
      endfor
      fig_h = [];

    endswitch

  endfor

  endwhile

  % 'x' -> Cleaning up: Reset fig_h for not accessing invalid handles
  % in the close callbacks of remaining figures. Reset the hold property,
  % following plots should delete this remaining one
  fig_h = [];
  hold off;

endfunction;



##
## Callback when closing figures with the closed loop step response.
## In this callback, the colors of the related closed loop poles together
## with their labels K, D, w0n are reeset to the 'normal' color
##
function __sim_fig_close_callback__ (h, e, handles, col_clp)

  global fig_h

  % Test whether fig_h is a valid array. If not, main figure may already
  % be closed
  if length (fig_h) == 0
    return; % Nothing to do here
  endif

  % Search for handle of closed figure in our array
  for j = 1:size (fig_h,1)
    [h_exists, h_idx] = ismember (h, fig_h(j,:));
    if h_exists
      fig_h(j,h_idx) = -1;  % Remove handle of closed figure
      K_row = j;
      break;
    endif
  endfor

  % If not found, return (should not happen)
  if ! exist ('K_row')
    return;
  endif

  % Is there another figure for this K?
  no_figure = 1;
  for j = 2:length (fig_h(K_row,:))
    if ishandle (fig_h(K_row,j))
      no_figure = 0;
      break;
    endif
  endfor

  if no_figure
    % No more figure, so reset color of closed loop poles and remove K
    for j = 1:length (handles)
      if ishandle (handles(j))
        set (handles(j), 'color', col_clp);
      endif
    endfor
    fig_h(K_row,:) = [];
  endif

endfunction
