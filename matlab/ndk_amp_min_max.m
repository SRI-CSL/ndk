%% Returns the min and max amplitude on a given channel.
function result = ndk_amp_min_max(esobj, channel)
    result = py.ndk.features.amp_min_max(esobj, channel);
end

