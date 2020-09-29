%% Just returns the version number by querying the python package:

function version = ndk()
    version = py.ndk.version();
end

