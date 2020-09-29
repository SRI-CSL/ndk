%% Returns the min and max value for all waveforms on the givel channel.
function pair = ndk_waveform_minmax(esobj, channel)
    pair = esobj.waveform_min_max(channel)
end
