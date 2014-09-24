function pass = check_direct_fit_params(params)

pass = 0;

props.NBAND.desc = 'The # of frequency sample points';


propsToCheck = cell(4, 1);
propsToCheck{1}.name = 'input';
propsToCheck{1}.msg = 'Must specify an input for spectrographic representation.';
propsToCheck{2}.name = 'samp_rate';
propsToCheck{2}.msg = 'Must specify a sample rate (Hz) for spectrographic representation.';
propsToCheck{3}.name = 'increment';
propsToCheck{3}.msg = 'Must specify a time increment (# samples) for spectrographic representation.';
propsToCheck{4}.name = 'window_length';
propsToCheck{4}.msg = 'Must specify a window length for spectrographic representation.';

%if (~check_props(params, propsToCheck))

