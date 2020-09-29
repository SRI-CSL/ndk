%% Returns an array of interspike intervals on the specified channel:
function r = ndk_isi(esobj, channel)
    tmp = py.ndk.features.isi(esobj,channel)
    r = double(py.array.array('d', py.numpy.nditer(tmp, pyargs('order','F'))))
end
