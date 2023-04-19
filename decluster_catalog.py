def decluster_catalog(catalog, space_threshold, days_threshold):
    """
            decluster_catalog is a function that is trying to find cluster of earthquake in a catalog based on closeness
            in space and time. Finally it removes all the events in the cluster except the maxium magnitude

            Parameters:
                catalog: catalog of earthquakes

                space_threshold: the threshold for space in km.

                days_threshold: the threshold for time in days.


            Returns:
                catalog: the declustered catalog
        """

    import numpy as np
    import geopy.distance


    # Get specific parts of the data in numpy array
    dates = catalog['date'].values
    lat = catalog['lat'].values
    lon = catalog['lon'].values
    mag = catalog['mw'].values

    # List of index of cluster
    cluster_EQ_index = list()

    # For each EQ in catalog
    for i in range(len(catalog)):
        cluster_temp = list()

        # Look if there is a similar EQ in close location in space and time
        for j in range(len(catalog)):

            # Calculate distance
            dist = geopy.distance.geodesic((lat[i], lon[i]), (lat[j], lon[j])).km
            # Calculate separation date
            time_shift_day = abs((dates[i] - dates[j]).astype('timedelta64[D]').astype('int'))

            # If two events are close in time and space
            if i != j and dist <= space_threshold and time_shift_day <= days_threshold:
                # Create cluster temp
                cluster_temp.append(j)

        # Cluster EQ index
        cluster_EQ_index.append(cluster_temp)

    # Create a list of indices to remove
    index_to_remove = list()

    # For all the element with possible cluster
    for i in range(0, len(cluster_EQ_index)):
        # If there is a cluster
        if cluster_EQ_index[i]:
            # Add the element to the cluster
            cluster_EQ_index[i].append(i)

            # Create an array with all the magnitudes
            mags = np.array(cluster_EQ_index[i])

            # Find the maximum argument (maximum magnitude)
            max_mag_index = np.argmax(mag[mags])

            # Remove the maximum argument from mags
            mags = np.delete(mags, max_mag_index, None)

            # Save the indices to delete
            for j in range(mags.shape[0]):
                index_to_remove.append(mags[j])

    # Create a set of indices to remove
    index_to_remove = list(set(index_to_remove))

    # Reset index numbering from 0
    catalog = catalog.reset_index()

    # Drop index from cluster
    catalog.drop(index_to_remove, axis=0, inplace=True)

    return catalog


