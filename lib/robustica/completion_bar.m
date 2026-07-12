function [times, old_done, last_text]= completion_bar(done, total, times, old_done, last_text)
  
% Indicates % of work done and remaining processing time.
%  
% Performs linear interpolation to estimate remaining time.  
%
%
% SYNTAX: [new_times, new_done, current_text] = completion_bar(done, total, old_times, old_done, last_text);
%
%
% OUTPUTS:
%
%        new_times: new vector of time values (one per row)
%
%        new_done : new vector of values of 'done' (one per row)
%
%        current_text: currently typed text (to be erased in the following calls).
%
%
% INPUTS:
%
%        done     : value indicating amount of work done, from 0 (before starting the job)
%                   to 'total' (when finishing the loop)
%
%        total    : value representing the total amount of work to be done
%
%        old_times: vector of previous time values; first call: empty
%
%        old_done : previous values of 'done'; first call: empty
%  
%        last_text: previously typed text (to be erased); first call: empty.
%
%
%
% Before starting the loop, call the function as:
%
%         [t0, done_vector, last_text] = completion_bar(0, total);  
%
% Then, at the end of each loop, use:
%
%         [t0, done_vector, last_text] = completion_bar(done, total, t0, done_vector, last_text); 
%
%
% HISTORY:
%
%    <Please add modification date here>: - <please add modification details here>
%
% -- 2014/11/21: Version 3 release ----------------------------------------------------------------
%
% -- 2010/02/16: Version 2 release ----------------------------------------------------------------
%
% 2010/02/10: - created and added to RobustICA package by Vicente Zarzoso
%               (I3S Laboratory, University of Nice Sophia Antipolis, CNRS, France).



% adjust percentage  
  
t = 'S|';
ast = fix(20*done/total); %in 5% steps
for k = 1:ast-1
  t = [t, '-'];
end
if ast > 0 
  t = [t, '>'];
end
for k = length(t)-1:20
  t = [t, ' '];
end
t = [t, '|E: ', num2str(5*ast),'% done'];


if exist('last_text')
   del_text= ['', 8*ones(1, length(last_text)+1)];  % to erase previous text
                                                    % (8 is the ASCII code for backspace)
else
   del_text = '';
end   


% adjust timer

et = clock;
tt = '';

if done ~= 0
 
  tt = [tt, ', '];

  et = etime(et, times(1:6));   % first 6 elements of time keep starting clock
  times = [times, et];          % keep record of delay values
  old_done = [old_done, done];  % add new 'done' reference
  
  T = length(old_done);
   
  tim = times(7:T+6);
  don = old_done;
   
  a = (mean(don.*tim) - mean(don)*mean(tim))/( mean(don.^2) - mean(don)^2 );
  b = mean(tim) - a*mean(don);
   
  eft = a*total + b;                % final time LS estimate
  ert = max([round(eft - et), 0]);  % remaining time estimate; always >= 0

  [days, hours, mins, secs] = days_hours_mins_secs(ert);
  
  if days ~= 0, tt = [tt, num2str(days), 'd ']; end
  if hours ~= 0, tt = [tt, num2str(hours), 'h '];  end
  if mins ~= 0, tt = [tt, num2str(mins), 'm '];  end
   
  tt = [tt, num2str(secs), 's remaining               ']; 
 
else % if done  
  times = [et, 0];  % initialize time and 'done' counters 
  old_done = 0;
  disp(' ');
  disp(['Simulation started on ', datestr(now)]);
end

t = [t, tt];
disp([del_text, t]);
last_text = t;

if done == total                                             
  [d, h, m, s] = days_hours_mins_secs(et);
  disp(['Simulation ended on   ', datestr(now)]);
  disp(['Total duration        ', num2str(d), ' days ', num2str(h), 'h ', num2str(m), 'm ', num2str(s), 's']);
end
              



function [d, h, m, s] = days_hours_mins_secs(secs)

% Converts a time duration given in seconds into days, hours, minutes and seconds.
%
% 
% SYNTAX: [d, h, m, s] = days_hours_mins_secs(secs);
%
%          d   : days
%          h   : hours
%          m   : minutes
%          s   : seconds
%
%          secs: time duration in seconds to be converted.


sd = 3600*24;   % seconds in a day
sh = 3600;      % seconds in an hour
sm = 60;        % seconds in a minute
  
d = fix(secs/sd);
secs = secs - d*sd;
h = fix(secs/sh);
secs = secs - h*sh;
m = fix(secs/sm);
s = round(secs - m*sm);
