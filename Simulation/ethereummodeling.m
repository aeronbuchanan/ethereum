
# An Octave (free matlab clone) simulation of ethereum (and bitcoin for that matter)
# https://www.gnu.org/software/octave/
#
# Copyright 2014 Aeron Buchanan
# Released under the Creative Commons Attribution-ShareAlike 3.0 License
# http://creativecommons.org/licenses/by-sa/3.0/

# Instructions:
# Download this script file somewhere
# Download and run [QT]Octave
# > pwd
# > cd path/to/script/directory
# > global identified
# > ethereummodeling()

# opmode can be
# 'simulation' [default]
# 'example'

function ethereummodeling(opmode)
	global identified;
	identified = false;

	if nargin < 1
		opmode = 'simulation'
	end

	# number of samples in time
	time_res = 1000;

	# number of nodes
	N = 100000;

	# number of different computer speeds
	Nsp = 100;

	# distribution of processor speeds
	# range of desktop speeds over last 3 years fastest/slowest = 10
	# assuming this accounts for 95% (2-sigma) then sigma = sqrt(sqrt(10)) = 1.78
	sigma = sqrt(sqrt(10));
	sp_mults = logspace( 6*log10(1/sigma), 6*log10(sigma), Nsp ); # 6-sigma range to get more than 99.999999% of speeds 
	freqs = stdnormal_pdf(linspace(-6, 6, Nsp));
	freqs = N * freqs / sum(freqs);
	# dither
	r = 0;
	for i = 1:length(freqs)
		v = floor(freqs(i));
		r = r + freqs(i) - v;
		if r > 1
			r = r - 1;
			v = v + 1;
		end
		freqs(i) = v;
	end
	# add a single fast node
	if sum(freqs) < N
		freqs(floor(Nsp*(9 + rand())/10)) = 1;
	end

	# target time = 1 min
	target_time = 1;

	# difficulty is chance of success per attempt = p 
	# chance of failure q = 1 - p

	# the simulation is based on large powers of the chance of failure q = 1 - p
	# and it is important that this is not rounded to 0 prematurely
	# but at the same time, must not be rounded to 1 before we even start:
	# smallest p = 1.11e-16 (2^-53) => (1 - p)^6.72e18 == 0
	# larger values of p (values of q further from 1) get rounded to 0 faster
	# smaller values of p are beyond double precision on my (64bit) machine

	# therefore the number of attempts for the fastest machine must be much
	# fewer than 6.7e18 = 2^62 to avoid rounding errors, which is fine

	# average speed of a processor (ignoring ASICS) is 
	# ~ 60 * 200 * 10^6 hash operations per minute
	# https://en.bitcoin.it/wiki/Mining_hardware_comparison
	# so to be able to work with the precision of my (64bit) machine
	# set minimum time interval to be an this "avg" time for a hash calculation 
	# (thus keeping calculations well within the 52 binary digits of precision)
	tbase = (2^-33.5); # of target_time 
	speeds = sp_mults; # hash operations per tbase

	# hash times now range from 2^-28.5 (slowest) to 2^-38.5 (fastest)
	# allowing for between 3.8e8 (fewest = 6 MH/s) and 3.8e11 (most = 6 GH/s) attempts within target time
	# which means success probs of 2^-17 (most difficult) to 2^-30 (least difficult)
	# will fail to be modeled accuractely due to premature rounding mentioned above
	# i.e difficulty must be between 2^-53 and 2^-32
	# so we might have to adjust the number of hashes of an atomic Eth PoW attempt
	diff_mult = 500;
	tbase = tbase * diff_mult;

	# difficulty to be such that expected time is t = 1
	# ... some complicated (unknown) function of N and speeds

	p = 3.5799711927710392416478381480891520546178119488445e-13; # 2^-41.3

	disp(['Starting difficulty (d = 1/p = 1/' num2str(p) ') = 2^' num2str(-log2(1 - (1 - p)^(1/diff_mult))) '.']);

	# chance of failure after n attempts = q^n
	# chance of failure after time t = q^(speed*t/tbase)
	# chance of any successes in N_n nodes of speed n after time t = 1 - q^(freqs(n)*speeds(n)*t/tbase)

	# calculate chance any node has succeeded after time t
	# (minimum time interval is target / (speeds(end) * tbase) )

	max_time = 1000 * target_time;
	ns = calc_ns(tbase, max_time, time_res);

	# toggle the commenting of the next two lines for different graphs 
	switch opmode
		case 'example'
			example(p, ns, tbase, target_time, speeds, freqs)
		case 'simulation'
			moving_diff_simulation(p, ns, speeds, freqs, target_time, tbase, 100000, 1000, 1000);
	end

	#keyboard
end

function moving_diff_simulation(p, ns, speeds, freqs, target_time, tbase, num_blocks, avg_win, update_period)

	# simulation
	difficulty_update_strategy = 'bitcoin';

	overhead = 2e3; # something negligible to start with
	proportion_comp = 0.0001; # this proportion of target_time on non-PoW comp
	# set proportion_comp to zero for faster sim times - it doesn't really affect anything

	# get a good starting value for p
	[pdf_t, cdf_t, pdf_sp_given_t, cdf_sp_given_t, pdf_sp, cdf_sp] = calc_dists(p, ns, target_time, tbase, speeds, freqs, overhead, proportion_comp);
	
	sp_delays = target_time * proportion_comp * speeds;
	broadcast_delay_for_average_node = sum(pdf_sp(:) .* sp_delays(:))
	delay = overhead + broadcast_delay_for_average_node / tbase;

	p = optimize_p(p, ns, target_time, tbase, speeds, freqs, delay, proportion_comp)

	figure(1); clf; hold on
	h_t = plot(0, 0, 'm');
	h_t_title = title('block times'); xlabel('number of weeks'); ylabel(['time (target = ' num2str(target_time) ')']);

	figure(2); clf
	#h_pdf_sp = tloglog(1, 1, 'c');
	h_pdf_sp = plot(1, 1, 'r');
	h_pdf_sp_title = title('pdf over speeds'); xlabel("speed as a multiple of the average node"); ylabel('probability density');

	figure(3); clf
	h_diffs = plot(0,0,'r');
	h_diffs_title = title('difficulty'); xlabel('number of weeks'); ylabel('probability of success');

	diffs = [p];
	t_times = [target_time];
	t_total = t_times;
	t_avgs = [0];

	[pdf_t, cdf_t, pdf_sp_given_t, cdf_sp_given_t, pdf_sp, cdf_sp] = calc_dists(p, ns, target_time, tbase, speeds, freqs, overhead, proportion_comp);

	while length(t_times) < num_blocks
		update_needed = false;

		diffs(end + 1) = p;

		# how long did this block take?
		qv = rand();
		[t_times(end + 1), index] = quartile(ns * tbase, cdf_t, qv);
		t_total(end + 1) = t_total(end) + t_times(end);

		# what speed processor got it?
		qv = rand();
		winning_speeds(end + 1) = quartile(speeds, cdf_sp_given_t(index,:), qv);

		# assume they did some non PoW ops
		if proportion_comp > 0
			old_delay = delay;
			delay = overhead + ceil(target_time * proportion_comp * winning_speeds(end) / tbase);
			if old_delay != delay
				update_needed = true;
			end
		end

		t_sample = t_times(max(1, length(t_times)-avg_win):end);
		t_avgs(end + 1) = sum(t_sample) / length(t_sample); 

		if length(t_avgs) > avg_win && mod(length(t_times), update_period) == 0
			# update diff
			p = update_(p, t_avgs(avg_win:end), target_time, difficulty_update_strategy);

			update_needed = true;
		end

		if update_needed
			[pdf_t, cdf_t, pdf_sp_given_t, cdf_sp_given_t, pdf_sp, cdf_sp] = calc_dists(p, ns, target_time, tbase, speeds, freqs, delay, proportion_comp);
		end

		if mod(length(t_times), num_blocks/200) == 0
			#update_loglog(h_pdf_sp, speeds, pdf_sp);
			set(h_pdf_sp, 'xdata', speeds, 'ydata', pdf_sp);
			set(h_pdf_sp_title, 'string', ['pdf over speeds (p = ' num2str(p) ', delay = ' num2str(delay*tbase) ')']);
			
			nn = 1000;
			if length(t_times) > nn
				bns = floor(linspace(1,length(t_avgs),nn));
				t_avgs_data = t_avgs(bns);
				diffs_data = diffs(bns);
			else
				bns = 1:length(t_avgs);
				t_avgs_data = t_avgs;
				diffs_data = diffs;
			end

			set(h_t_title, 'string', ['block times (latest avg = ' num2str(t_avgs(end)) ')']);
			set(h_t, 'ydata', t_avgs_data, 'xdata', t_total(bns)/(60*24*7));

			set(h_diffs_title, 'string', ['difficulty (latest p = ' num2str(p) ')']);
			set(h_diffs, 'ydata', diffs_data, 'xdata', t_total(bns)/(60*24*7));
			pause(0)
		end

	end

	mean_t_times = mean(t_times)
	median_t_times = median(t_times)
	mode_t_times = mode(t_times)
	std_t_times = std(t_times)

	mean_t_avgs = mean(t_avgs(avg_win:end))
	median_t_avgs = median(t_avgs(avg_win:end))
	mode_t_avgs = mode(t_avgs(avg_win:end))
	std_t_avgs = std(t_avgs(avg_win:end))

	figure(1)
	lims = axis;
	text(0.4*lims(2), 0.4, ['mean avg-' num2str(avg_win) ' upd-' num2str(update_period) ' = ' num2str(mean_t_avgs)]);
	text(0.4*lims(2), 0.3, ['std avg-' num2str(avg_win) ' upd-' num2str(update_period) ' = ' num2str(std_t_avgs)]);
	text(0.4*lims(2), 0.2, ['median avg-' num2str(avg_win) ' upd-' num2str(update_period) ' = ' num2str(median_t_avgs)]);
	text(0.4*lims(2), 0.1, ['mode avg-' num2str(avg_win) ' upd-' num2str(update_period) ' = ' num2str(mode_t_avgs)]);

	text(0.1*lims(2), 0.4, ['mean = ' num2str(mean_t_times)]);
	text(0.1*lims(2), 0.3, ['std = ' num2str(std_t_times)]);
	text(0.1*lims(2), 0.2, ['median = ' num2str(median_t_times)]);
	text(0.1*lims(2), 0.1, ['mode = ' num2str(mode_t_times)]);

	if proportion_comp > 0
		text(0.7*lims(2), 0.1, ['constant relative comp of target time / ' num2str(proportion_comp)]);
	end

	#keyboard
end

function new_p = update_(p, avgs, target, strategy)
	global identified;

	switch strategy

		case "bitcoin"
			if ! identified, printf('using bitcoin update strategy'), identified = true; end
			# the bitcoin update
			update_factor = avgs(end) / target;
			update_factor = max(0.25, min(4, update_factor));
			new_p = p * update_factor;

		case "PID"
			if ! identified, printf('using PID update strategy'), identified = true; end
			# PID controller
			error = target - avgs(end);
			integral = sum(target - avgs);
			derivative = avgs(end) - avgs(end - 1);
			# parameters estimated using wikipedia PID article
			K_p = 3e7; # instability by 1e8
			K_i = 1.2 * K_p / 2200; # periodic instability
			K_d = 1e-3;
			# because p is bounded
			# work on a notion of difficulty d = 1/p
			update = (K_p * error + K_i * integral + K_d * derivative);
			new_p = p / (1 + p * update);

		case "proportional"
			if ! identified, printf('using proportional update strategy'), identified = true; end
			# proportional-style controller 
			update = (avg - target) / 100;
			new_p = p * (1 + update);

		case "whitepaperV2"
			if ! identified, printf('using V2 update strategy'), identified = true; end
			# ethereum whitepaper version 2
			# note: p = 1/difficulty
			# note: for small p, diff_mult does not have a significant effect
			# note: this is effectively the same as the default below
			if avgs(end) > target
				new_p = p * (1000 / 999); # diff -= floor(diff/1000)
			else
				new_p = p * (1000 / 1001); # diff += floor(diff/1000)
			end

		otherwise
			if ! identified, disp('using default update strategy'), identified = true; end
			update = p / 1000;

			if update < 2^-53
				# simulation is insensitive to any part of p < 2^-53
				disp('Warning! round update up to 2^-53');
				update = 2^-53; # could be doubling the difficulty!
			end

			if avg > target
				new_p = p + update; # easier
			else
				new_p = p - update; # harder
			end
	end

	# p is a probability
	new_p = max(0, min(1, new_p));
end


function p = optimize_p(p, ns, target_time, tbase, speeds, freqs, delay, percentage_non_pow)
	# optimize p (binary search)
	ts = ns * tbase;

	margin = 0.005;

	pdf_t = calc_dists(p, ns, target_time, tbase, speeds, freqs, delay, percentage_non_pow);
	mean_t = find_mean(ts, pdf_t);

	old_p = p;
	while mean_t < target_time * (1 - margin) || mean_t > target_time * (1 + margin)
		new_p = p;
		if mean_t > target_time * (1 - margin)
			if new_p > old_p
				update = (new_p - old_p) * 1.5;
			else
				update = (old_p - new_p) / 2;
			end
			new_p = new_p + max(update, 2^-53); # easier 
		elseif mean_t < target_time * (1 + margin)
			if new_p > old_p
				update = (new_p - old_p) / 2;
			else
				update = (old_p - new_p) * 1.5;
			end
			new_p = new_p - max(update, 2^-53); # harder
		end
		old_p = p;
		p = new_p;

		pdf_t = calc_dists(p, ns, target_time, tbase, speeds, freqs, delay, percentage_non_pow);
		mean_t = find_mean(ts, pdf_t);

		if p > 1
			return
		end

		printf('.');
	end
	printf(' ')
end

function example(p, ns, tbase, target_time, speeds, freqs)

	overhead = 2000; # small number 

	for proportion_comp = linspace(0,1,10)
		[pdf_t, cdf_t, pdf_sp_given_t, cdf_sp_given_t, pdf_sp, cdf_sp] = calc_dists(p, ns, target_time, tbase, speeds, freqs, overhead, proportion_comp);
		
		pdf_sps(:, end+1) = pdf_sp;
		sp_delays = target_time * proportion_comp * speeds;
		this_delay_for_average_node = proportion_comp
		next_delay_for_average_node = sum(pdf_sp(:) .* sp_delays(:))

		mean_t = find_mean(ns, pdf_t) * tbase;
		mode_sp = find_mode(speeds, pdf_sp);
		disp([num2str(proportion_comp*100) '% comp: mean = ' num2str(mean_t) ', modal speed = ' num2str(mode_sp) ])
	end

	proportion_comp = 1e-5; 
	[pdf_t, cdf_t, pdf_sp_given_t, cdf_sp_given_t, pdf_sp, cdf_sp] = calc_dists(p, ns, target_time, tbase, speeds, freqs, overhead, proportion_comp);

	ts = ns * tbase;
	
	figure(1)
	clf
	hold on
	tloglog(ts, pdf_t, 'b');
	title('probability that PoW will be solved in a particular time slot')
	xlabel("time as a fraction of target time")
	ylabel('probability density')
	# add line for target_time
	lims = axis;
	plot([target_time, target_time], lims([3 4]), 'y')

	figure(2)
	clf
	hold on
	tloglog(ts, cdf_t, 'b');
	title('cumulative probability that PoW will be solved in a particular time slot')
	xlabel("time as a fraction of target time")
	ylabel('cumulative probability')
	# add line for target_time
	lims = axis;
	plot([target_time, target_time], lims([3 4]), 'y')

	figure(3)
	clf
	hold on
	tloglog(speeds, freqs / sum(freqs), 'c');
	tloglog(speeds, pdf_sp, 'k');
	title("probability that PoW will be solved by a particular speed node\n(note shift from distribution of speeds, in cyan)")
	xlabel('speed multiple')
	ylabel('probability density')
	# add line for average speed
	lims = axis;
	plot([1 1], lims([3 4]), 'y')

	figure(4)
	clf
	hold on
	tloglog(speeds, cdf_sp);
	title('cumulative probability that PoW will be solved by a particular speed node')
	xlabel('speed multiple')
	ylabel('cumulative probability')
	# add line for average speed
	lims = axis;
	plot([1 1], lims([3 4]), 'y')

	figure(5)
	clf
	hold on
	tloglog(speeds(sum(pdf_sp_given_t,1)!=0), pdf_sp_given_t(:, sum(pdf_sp_given_t,1)!=0));
	title('probability that PoW will be solved in a particular time slot by a particular speed node')
	xlabel('speed multiple')
	ylabel('probability density')
	# add line for average speed
	lims = axis;
	plot([1 1], lims([3 4]), 'y')

	figure(6)
	clf
	hold on
	plot(speeds, cdf_sp_given_t(sum(cdf_sp_given_t,2)!=0,:));
	title('cumulative probability that PoW will be solved for various times')
	xlabel('speed multiple')
	ylabel('cumulative probability')
	# add line for average speed
	lims = axis;
	plot([1 1], lims([3 4]), 'y')
	n = 1;
	used_spis = [];
	while n < length(ns)
		if sum(cdf_sp_given_t(n,:)) > 0
			[spi, i] = quartile(speeds, cdf_sp_given_t(n,:), 0.5);
			if length(find(used_spis==spi)) == 0
				text(spi, min(0.5, cdf_sp_given_t(n,i)), ['t = ' num2str(n)]);
				used_spis(end+1) = spi;
			end
			n = floor(n * 1.4);
		else
			n = n + 1;
		end
	end

	A = trapezium_integration(ts, pdf_t)
	mean_t = find_mean(ns, pdf_t) * tbase
	median_t = find_median(ns, cdf_t) * tbase
	mode_t = find_mode(ns, pdf_t) * tbase
	skew_t = skewness(ts, cdf_t)
	
	disp('test')
	xs = logspace(-10,2,500);
	xs = [-1e5 -fliplr(xs) 1 xs 1e5];
	cdf = stdnormal_cdf(xs);
	pdf = stdnormal_pdf(xs);
	xs = xs + 2;

	A = trapezium_integration(xs, pdf)
	mean = find_mean(xs, pdf)
	median = find_median(xs(2:end), cdf(1:end-1))
	mode = find_mode(xs, pdf)
	skew = skewness(xs, cdf)

	disp(["Type 'dbquit' to exit debug mode"])
	keyboard
end


function h = tloglog(xs, ys, opts)
	if nargin < 3
		opts = '-';
	end

	if ( length(xs) == size(ys,1) && size(ys,2) > 1)
		check = sum(ys != 0, 1);
		ys = ys(:, check != 0);
	elseif  ( length(xs) == size(ys,2) && size(ys,1) > 1)
		check = sum(ys != 0, 2);
		ys = ys(check != 0, :);
	end

	ys( ys == 0 ) = NaN;
	xs( xs == 0 ) = NaN;

	h = loglog(xs, ys, opts);
end

function update_loglog(handle, xs, ys)
	is = (xs != 0) & (ys != 0);
	xxs = xs( is );
	yys = ys( is );

	set(handle, 'xdata', xxs, 'ydata', yys);
end

function c = random_category(probs) 
	ri = rand();
	cps = [0 cumsum(probs)];
	cps(end) = 10; # must not be chosen
	is = find(cps < ri);
	c = is(end);
end

function ns = calc_ns(tbase, max_time, time_res)
	ns = floor(logspace(0, log10(max_time / tbase), time_res));
end

function [pdf_t, cdf_t, pdf_sp_given_t, cdf_sp_given_t, pdf_sp, cdf_sp] = calc_dists(p, ns, target_time, tbase, speeds, freqs, delay, proportion_comp)

	# delay = number of Eth PoW ops used on non PoW work (verify + compute)

	N = length(ns);
	Nsp = length(speeds);

	ts = ns * tbase;

	# calculate number of attempts for an average node
	delay = delay + target_time * proportion_comp / tbase;
	ns = max(0, ns - delay);	

	n_attempts = floor(ns(:) * speeds(:)'); # per node in each speed category
	n_attempts_sp = n_attempts .* (ones(N,1) * freqs(:)'); # per speed category
	n_attempts_total = sum(n_attempts_sp,2); # by a particular time
	p_all_fail = (1 - p) .^ n_attempts_total; # prob that all nodes in a given speed category fail

	# cumulative prob that at least one node got it
	cdf_t = 1 - p_all_fail; 

	# prob that it will be solved in a particular time slot
	pdf_t = [0; diff(cdf_t(:))] ./ [1; diff(ts(:))];
	pdf_t( isnan(pdf_t) ) = 0;

	# prob that a node of a given speed got it in a particular time slot
	pdf_sp_given_t = (pdf_t * ones(1,Nsp)) .* n_attempts_sp;
	pdf_sp_given_t = pdf_sp_given_t ./ (sum(pdf_sp_given_t,2) * ones(1,Nsp));
	pdf_sp_given_t( isnan(pdf_sp_given_t) ) = 0;

	cdf_sp_given_t = cumsum(pdf_sp_given_t,2);

	pdf_sp_and_t = pdf_sp_given_t .* (pdf_t * ones(1,Nsp));
	pdf_sp = sum(pdf_sp_and_t .* ([0 diff(ts)]' * ones(1,Nsp)),1);

	cdf_sp = cumsum(pdf_sp);

end

function [f, f1] = chance_of_failure(p, speeds, freq, ns)
	# discrete version (slow and unnecessary)
	
	attempts = floor(ns * speed);
	[as, is, js] = unique(attempts);

	if freq * (as(end) + 1) > 6.72e18
		disp(['Warning! ' num2str(ns(end)) '!'])
	end

	q = 1 - p;
	uf = q .^ (as * freq); # cdf after n attempts
	uf1 = q .^ ((as + 1) * freq); # cdf after n+1 attempts

	f = uf(js); # "full" cdf
	f1 = uf(js); # "full" cdf
	f1(is) = uf1; # fill in transition indices with the next cdf value
end

function A = trapezium_integration(xs, ys)
	# area under the graph using linear interpolation
	A = sum( trapezia_areas(xs, ys) );
end

function As = trapezia_areas(xs, ys)
	xs = reshape(xs, 1, []);
	ys = reshape(ys, 1, []);
	yys = 0.5 * (ys(1:end-1) + ys(2:end));
	hs = diff(xs);
	As = yys .* hs;
end

function [x, q] = quartile(xs, cdf, p)
	qis = find(cdf > p);
	if ( length(qis) > 0 )
		q = qis(1);
		x = xs(q);
	else
		x = xs(end);
	end
end

function m = find_median(xs, cdf)
	m = quartile(xs, cdf, 0.5);
end

function m = find_mean(xs, pdf)
	probs = trapezia_areas(xs, pdf);
	S = sum(probs);
	if S > 0
		probs = probs / S;
		mid_xs = 0.5 * (xs(1:end-1) + xs(2:end));
		m = sum( mid_xs .* probs );
	else
		m = NaN;
	end
end

function m = find_mode(xs, pdf)
	[p, i] = max(pdf);
	m = xs(i);	
end

function s = skewness(xs, cdf)
	u = 0.75;
	s = (quartile(xs, cdf, u) + quartile(xs, cdf, 1 - u) - 2*quartile(xs, cdf, 0.5) ) / (quartile(xs, cdf, u) + quartile(xs, cdf, 1 - u));
end


