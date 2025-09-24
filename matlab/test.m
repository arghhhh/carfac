classdef test
    methods( Static )
        
function status = test_CAR_freq_response(do_plots)
% Test: Make sure that the CAR frequency response looks right.

status = 0;

CF = CARFAC_Design(1, 22050);  % 1 ear; use old test fs.
CF = CARFAC_Init(CF);
n_points = 2^14;
fft_len = n_points * 2;  % Lots of zero padding for good freq. resolution.
impulse_responses = zeros(fft_len, CF.n_ch);
ear = CF.ears(1);
ear.CAR_coeffs.linear = 1;
impulse = 1;
for i = 1:n_points
  [car_out, state] = CARFAC_CAR_Step( ...
    impulse, ear.CAR_coeffs, ear.CAR_state);
  ear.CAR_state = state;
  impulse = 0;
  impulse_responses(i, :) = car_out';
end
complex_spectra = fft(impulse_responses);
db_spectra = 20 * log10(abs(complex_spectra) + 1e-50);

if do_plots
  figure
  db_spectra(db_spectra < -20) = -20;
  imagesc(db_spectra(1:(n_points+1), :))
  set(gca,'YDir','normal')
  colorbar()
  xlabel('Channel Number')
  ylabel('Frequency bin')
  drawnow
end

expected = [ ...
  [10, 5604, 39.34, 705.6, 7.9];
  [20, 3245, 55.61, 429.8, 7.6];
  [30, 1809, 60.46, 248.1, 7.3];
  [40, 965, 59.18, 138.7, 7.0];
  [50, 477, 52.81, 74.8, 6.4];
  [60, 195, 38.98, 37.5, 5.2];
  [70, 32, 7.95, 14.6, 2.2];
  ];

% Test: check overall frequency response of a cascade of CAR filters.
% Match Figure 16.6 of Lyon's book
spectrum_freqs = (0:n_points)' * CF.fs / fft_len;
for i = 1:size(expected, 1)
  channel = expected(i, 1) + 1;  % Channel number was 0 based from Python.
  correct_cf = expected(i, 2);
  correct_gain = expected(i, 3);
  correct_bw = expected(i, 4);
  correct_q = expected(i, 5);
  cf_amp_bw = test.find_peak_response(spectrum_freqs, ...
    db_spectra(:, channel), 3);  % 3 dB width
  cf = round(cf_amp_bw(1));  % Zero decimals on this one.
  % Round to 1 or 2 decimal places, require exact match then.
  % 2 Decimal places for gain.
  gain = round(cf_amp_bw(2), 2);
  bw = round(cf_amp_bw(3), 1);
  q = round(cf_amp_bw(1) / cf_amp_bw(3), 1);
  fprintf(1, ...
    ['%d: cf is %.1f Hz, peak gain is %.1f dB,' ...
    ' 3 dB bandwidth is %.1f Hz (Q = %.1f).\n'], ...
    channel, cf, gain, bw, q);
  if cf ~= correct_cf
    status = 1;
    fprintf(1, 'Mismatch cf %f should be %f.\n', cf, correct_cf)
  end
  if gain ~= correct_gain
    status = 1;
    fprintf(1, 'Mismatch gain %f should be %f.\n', gain, correct_gain)
  end
  if bw ~= correct_bw
    status = 1;
    fprintf(1, 'Mismatch bw %f should be %f.\n', bw, correct_bw)
  end
  if q ~= correct_q
    status = 1;
    fprintf(1, 'Mismatch q %f should be %f.\n', q, correct_q)
  end
end
test.report_status(status, 'test_CAR_freq_response')
return

end


function cf_amp_bw = find_peak_response(freqs, db_gains, bw_level)
% Returns center frequency, amplitude at this point, and the 3dB width."""
[~, peak_bin] = max(db_gains);
[peak_frac, amplitude] = test.quadratic_peak_interpolation( ...
  db_gains(peak_bin - 1), db_gains(peak_bin), db_gains(peak_bin + 1));
peak_loc = peak_bin + peak_frac;
if peak_frac < 0
  cf = (1 + peak_frac)*freqs(peak_bin) - ...
    peak_frac*freqs(peak_bin - 1);
else
  cf = (1 - peak_frac)*freqs(peak_bin) + ...
    peak_frac*freqs(peak_bin + 1);
end
% cf = linear_interp(freqs, peak_bin + peak_frac)
freqs_3db = test.find_zero_crossings(freqs, db_gains - amplitude + bw_level);
if length(freqs_3db) >= 2
  bw = freqs_3db(2) - freqs_3db(1);
else
  bw = 0;
end
cf_amp_bw = [cf, amplitude, bw];
return
end


% Formula from:
% https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html
function [location, amplitude] = quadratic_peak_interpolation(...
  alpha, beta, gamma)
location = 1 / 2 * (alpha - gamma) / (alpha - 2 * beta + gamma);
amplitude = beta - 1.0 / 4 * (alpha - gamma) * location;
return
end


function index = find_closest_channel(values, desired)
[~, index] = min((values - desired).^2);
return
end

function zclist = find_zero_crossings(x, y)
locs = find(y(2:end) .* y(1:end-1) < 0, 2);
a = y(locs);
b = y(locs+1);
frac = -a ./ (b - a);
zclist = x(locs) .* (1 - frac) + x(locs + 1) .* frac;
return
end


function status = test_spike_rates(do_plots)
% Test: Assure the 3 class rates versus level look good.

status = 0;
fs = 22050;
fp = 1000;  % Probe tone frequency
duration = 0.25;
dbstep = 10;   % 10 is good
dbfs = -104:dbstep:6;  % 0 to 110 dB SPL

t = (0:(1/fs):(duration - 1/fs))';  % Sample times for short duration
sinusoid = sin(2 * pi * t * fp);
signal = [];
time = [];
t_start = 0;
for db = dbfs  % Levels spanning a huge range
  amplitude = sqrt(2) * 10.^(db/20);
  signal = [signal; amplitude*sinusoid];
  time = [time; t + t_start];
  t_start = t_start + duration;
end

CF = CARFAC_Design(1, fs, 'do_syn');  % v3 3-class synapse model
CF = CARFAC_Init(CF);
[nap, CF, bm, ohc, agc, firings] = CARFAC_Run_Segment(CF, signal);  % nap has 3 columns of firings

if do_plots
  %%
  chan = find(CF.pole_freqs * 1.06 < fp, 1); % probably best channel
  chan_firings = squeeze(firings(:, chan, :, 1));  % Just one channel, 3 class columns.
  healthy_n_fibers = CF.SYN_params.healthy_n_fibers;
  rates = chan_firings ./ healthy_n_fibers;

  figure();
  plot(time, chan_firings);
  title('Instantaneous rates of 3 fiber-group classes')
  xlabel('time in seconds, with 10 dB steps from -100 to 0 dB FS')
  ylabel('firings per sample')
  for db = dbfs + 104
    text(duration * (db/dbstep + 0.4), 12, num2str(db))
  end

  figure();
  plot(time, fs*smooth1d(rates, fs*0.005)) % Per fiber
  title('Mean rates of 3 fiber classes')
  xlabel('time in seconds, with 10 dB steps from 0 to 110 dB SPL rms')
  ylabel('firings per second per fiber')
  for db = dbfs + 104
    text(duration * (db/dbstep + 0.4), 100, num2str(db))
  end
  octave_basal_chan = find(CF.pole_freqs * 1.06 < fp*2, 1);
  half_octave_basal_chan = find(CF.pole_freqs * 1.06 < fp*sqrt(2), 1);
  best_chan = find(CF.pole_freqs * 1.06 < fp, 1);
  half_octave_apical_chan = find(CF.pole_freqs * 1.06 < fp/sqrt(2), 1);
  channels = [octave_basal_chan, half_octave_basal_chan, best_chan, ...
    half_octave_apical_chan];
  figure()
  plot(time(1:8:end), agc(1:8:end, channels))
  text(2.55, 0.15, [num2str(channels(4)), ': apical 0.5'])
  text(2.55, 0.5, [num2str(channels(3)), ': best'])
  text(2.58, 0.74, [num2str(channels(2)), ': basal 0.5'])
  text(2.45, 0.93, [num2str(channels(1)), ': basal 1'])
  for db = dbfs + 104
    text(duration * (db/dbstep + 0.4), 0.4, num2str(db))
  end
end

report_status(status, 'test_spike_rates')
return
end


function report_status(status, name, extra)
if nargin < 3, extra = 0; end
if extra
  if status
    disp(['FAIL ' name '; at least one test failed.'])
  else
    disp(['PASS ' name '; all tests passed.'])
  end
else
  if status
    disp(['FAIL ' name])
    if status > 1
      disp('(status > 1 => error in test or expected results size)')
    end
  else
    disp(['PASS ' name])
  end
end
return
end




    end % methods( Static )
end % classdef test

