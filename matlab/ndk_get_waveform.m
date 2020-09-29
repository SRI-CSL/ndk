%% Returns an array that contains the waveform for the specified
%% spike (spike_id) on the specified channel:
function y = ndk_get_waveform(esobj, spike_id, channel)
    arr = esobj.get_spike_waveform(spike_id, channel)
    y = double(py.array.array('d', py.numpy.nditer(arr, pyargs('order','F'))))
end
