%% Returns an event store object.  Can be closed with esobj.close():

function esobj = ndk_open(filename)
    esobj = py.ndk.es.db3.db3(filename);
end

