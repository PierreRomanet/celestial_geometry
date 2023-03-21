def find_moon_phase(date, threshold_angle, model_ephemeris, other_info=''):
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
    from astropy.coordinates import get_body_barycentric, get_body_barycentric_posvel
    import numpy as np

    # Conjunction list
    moon_phase = list()

    # Bodies
    celestial_bodies = ['Moon', 'Earth', 'Sun']

    # Create solar system
    solar_system = {}

    # Get position of every body
    for body in celestial_bodies:
        solar_system[body] = get_body_barycentric(body, date, ephemeris=model_ephemeris).get_xyz().value

    # Get velocity of the Earth
    vec_vel_earth = get_body_barycentric_posvel('Earth', date, ephemeris=model_ephemeris)[1].get_xyz().value
    vec_vel_earth = vec_vel_earth/np.linalg.norm(vec_vel_earth)

    # Get vector sun to moon, and sun to earth
    vec_earth_sun = solar_system['Sun'] - solar_system['Earth']
    vec_earth_sun = vec_earth_sun/np.linalg.norm(vec_earth_sun)

    # Project the moon in the Ecliptic plan
    plan_vec = np.cross(vec_vel_earth, vec_earth_sun)

    moon_proj = solar_system['Moon']-solar_system['Earth'] - plan_vec * (
                (np.dot(plan_vec, solar_system['Moon']) - np.dot(plan_vec, solar_system['Earth'])) / np.dot(plan_vec, plan_vec))

    # print(np.dot(moon_proj,plan_vec)-np.dot(plan_vec,solar_system['Earth']))
    # Calculate angle
    angle = 180 / np.pi * np.arccos(np.dot(moon_proj / np.linalg.norm(moon_proj), vec_earth_sun ))

    # angle = 180 / np.pi * np.arccos(
    #         np.dot(vec_earth_moon / np.linalg.norm(vec_earth_moon), vec_earth_sun / np.linalg.norm(vec_earth_sun)))

    # If there is a conjonction
    if 180 - angle < threshold_angle:
        moon_phase.append('Moon-Earth-Sun     {}'.format(str(other_info)))
    elif angle < threshold_angle:
        moon_phase.append('Earth-Moon-Sun     {}'.format(str(other_info)))

    return moon_phase
