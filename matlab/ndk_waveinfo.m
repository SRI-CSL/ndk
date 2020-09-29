%% Returns a 4-tuple of peak-to-trough amplitude difference,
%% peak-to-trough time difference, min amplitude, and max amplitude
%% for a given spike ID (index) on a given channel.
function tuple = ndk_waveinfo(esobj, spikeID, channel)
    tuple = py.ndk.features.waveinfo(esobj, spikeID, channel)
end
