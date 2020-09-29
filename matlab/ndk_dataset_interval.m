%% Returns the min and max timestamps for events in the dataset.
function tuple = ndk_dataset_interval(esobj, channel)
    tuple = esobj.get_dataset_interval();
end
