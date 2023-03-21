def find_conjunction(date,celestial_bodies,threshold_angle,model_ephemeris,other_info=''):
    """
        Find_conjunction is a function that is returning a list of conjunction for a given date.
        Parameters:
            date: a date in astropy.time.Time format.

            celestial_bodies: list of celestial objects among
            ['moon', 'mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune', 'sun'].

            threshold_angle: the threshold angle for defining a conjunction, in degree.

            model_ephemeris: the name of a model from NASA JPL Horizons ephemeris ('de430' or 'de432s')

        Returns:
            conjunction_list: a list of conjunctions
    """
    from astropy.coordinates import get_body_barycentric
    import itertools
    import numpy as np

    # Conjunction list
    conjunction_list = list()

    # Create solar system
    solar_system = {}

    # Get position of every body
    for body in celestial_bodies:
        solar_system[body] = get_body_barycentric(body, date, ephemeris=model_ephemeris).get_xyz().value

    # Find all combinaison of 3 planets among all the planets
    for triplet in itertools.combinations(solar_system.keys(), 3):
        # Create list of vector
        vec = list()

        # Calculate vector 1
        vec.append(solar_system[triplet[1]] - solar_system[triplet[0]])

        # Calculate vector 2
        vec.append(solar_system[triplet[2]] - solar_system[triplet[1]])

        # Calculate vector 3
        vec.append(solar_system[triplet[0]] - solar_system[triplet[2]])

        # Find maximum norm
        max_norm_i = np.argmax(np.array([np.linalg.norm(vec[0]), np.linalg.norm(vec[1]), np.linalg.norm(vec[2])]))

        # Make scalar product with the two smallest norm
        del vec[max_norm_i]
        angle = 180 / np.pi * np.arccos(np.dot(vec[0] / np.linalg.norm(vec[0]), vec[1] / np.linalg.norm(vec[1])))

        # If there is a conjonction
        if 180 - angle < threshold_angle or angle < threshold_angle:
            if max_norm_i == 0:
                conjunction_list.append('{4}    {0}-{1}-{2}    Angle:{3}'.format(triplet[0], triplet[2], triplet[1], angle, str(other_info)))
            elif max_norm_i == 1:
                conjunction_list.append('{4}    {0}-{1}-{2}    Angle:{3}'.format(triplet[1], triplet[0], triplet[2], angle, str(other_info)))
            elif max_norm_i == 2:
                conjunction_list.append('{4}    {0}-{1}-{2}    Angle:{3}'.format(triplet[0], triplet[1], triplet[2], angle, str(other_info)))

    return conjunction_list

