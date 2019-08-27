'''This file is for all of the satellite functions to be in one place to not have to keep staying
in the same notebook'''
'''Quick snapshot of all functions in this file
load_satellites
get_satellites
from_localtime
separation_matrix # this is from another file but is needed for other functions
get_pos
get_visible
radec
'''
#importing all necessary functions
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pylab, mlab
import scipy.interpolate
import skyfield.api
import skyfield.toposlib
from skyfield.api import Time
import datetime

#Globals
ts=skyfield.api.load.timescale()


def load_satellites(reload = True, groups = None, debris = True):
	"""reload defaulted to true, when running in jupyter set to False.
	When groups is set to None, all of the groups from celestrak will be included in the dictionary. To specify groups, find the filename from celestrak and set 
	groups = ('group1','group2', 'groupx'). 
	Debris will be included if groups = None OR if debris = True. Debris will be default to True, and will only be excluded if debris=False. 
	If specifying groups, do not set debris =True or False."""
	satellites = {}
	added = set()
	if groups is not None:
		for group in groups:
			loaded=skyfield.api.load.tle(f'http://celestrak.com/NORAD/elements/{group}.txt', reload=reload)
			unique=set(loaded.values())
			toadd=unique-added
			nadded=0
			for satellite in toadd:
				if satellite.target_name in satellites:
					if satellite.epoch.tai <= satellites[satellite.target_name].epoch.tai:
						continue
				satellites[satellite.target_name] = satellite
				added.add(satellite)
				nadded+=1
			print(f'Added {nadded} satellites from {group}.')
		print(f'Found {len(satellites)} satellites')
		return satellites
	if groups is None:
		if debris is True:
			for group in (
			# Weather & earth resources
			'weather', 'noaa', 'goes', 'resource', 'sarsat', 'dmc', 'tdrss', 'argos', 'planet', 'spire',
			# Communications
			'geo', 'gpz', 'gpz-plus', 'intelsat', 'ses', 'iridium', 'iridium-NEXT', 'starlink', 'orbcomm', 'globalstar',
			'amateur', 'x-comm', 'other-comm', 'satnogs', 'gorizont', 'raduga', 'molniya',
			# Navigation
			'gps-ops', 'glo-ops', 'galileo', 'beidou', 'sbas', 'nnss', 'musson',
			# Scientific
			'science', 'geodetic', 'engineering', 'education',
			# Miscellaneous
			'military', 'radar', 'cubesat', 'other',
			# Debris
			'2019-006', '1999-025', 'iridium-33-debris', 'cosmos-2251-debris', '2012-044',
			# Special interest (should already be included above)
			'tle-new','stations', 'visual', 'active', 'analyst'):
				loaded=skyfield.api.load.tle(f'http://celestrak.com/NORAD/elements/{group}.txt', reload=reload)
				unique=set(loaded.values())
				toadd=unique-added
				nadded=0
				for satellite in toadd:
					if satellite.target_name in satellites:
						if satellite.epoch.tai <= satellites[satellite.target_name].epoch.tai:
							continue
					satellites[satellite.target_name] = satellite
					added.add(satellite)
					nadded+=1
				print(f'Added {nadded} satellites from {group}.')
			print(f'Found {len(satellites)} satellites')
			return satellites
		if debris is False:
			for group in (
			# Weather & earth resources
			'weather', 'noaa', 'goes', 'resource', 'sarsat', 'dmc', 'tdrss', 'argos', 'planet', 'spire',
			# Communications
			'geo', 'gpz', 'gpz-plus', 'intelsat', 'ses', 'iridium', 'iridium-NEXT', 'starlink', 'orbcomm', 'globalstar',
			'amateur', 'x-comm', 'other-comm', 'satnogs', 'gorizont', 'raduga', 'molniya',
			# Navigation
			'gps-ops', 'glo-ops', 'galileo', 'beidou', 'sbas', 'nnss', 'musson',
			# Scientific
			'science', 'geodetic', 'engineering', 'education',
			# Miscellaneous
			'military', 'radar', 'cubesat', 'other',
			# Special interest (should already be included above)
			'tle-new','stations', 'visual', 'active', 'analyst'):
				loaded=skyfield.api.load.tle(f'http://celestrak.com/NORAD/elements/{group}.txt', reload=reload)
				unique=set(loaded.values())
				toadd=unique-added
				nadded=0
				for satellite in toadd:
					if satellite.target_name in satellites:
						if satellite.epoch.tai <= satellites[satellite.target_name].epoch.tai:
							continue
					satellites[satellite.target_name] = satellite
					added.add(satellite)
					nadded+=1
				print(f'Added {nadded} satellites from {group}.')
			print(f'Found {len(satellites)} satellites')
			return satellites


def get_satellites(satellites, string):
	''' Returns the satellite data(name, type?, target_name, number, and epoch, for any
	satellite with the string in its name.
	Parameters
	----------
	string: string, example: 'ISS'
	Returns
	-------
	List of satellites data (name, type?, target_name, number, and epoch)
	with 'string' in its name.
	Notes 
	-----
	Must have defined satellites from load_satellites first
	any satellite with the string in its name will be returned. For example,
	both ISS DEB and ISS (ZARYA) will be returned with input as 'ISS'.'''
	found = [name for name in satellites if string in name]
	return {name: satellites[name] for name in found}


def from_localtime(*args, tzoffset = -7):
	ts=skyfield.api.load.timescale()
	return ts.tai_jd(ts.utc(*args).tai - tzoffset / 24)


def get_pos(where, when, string, satellites):
    which = get_satellites(satellites, string)
    for satellite in which.values():
        # Calculate satellite position relative to the observer at each time step.
        topocentric = (satellite - where).at(when)
        # Calculate (alt, az)
        alt, az, _ = topocentric.altaz()
        # Calculate astrometric (ra, dec)
        ra, dec, _ = topocentric.radec()
        print(satellite.name, 'ALT', alt.degrees, 'AZ', az.degrees, 'RA', ra._degrees, 'DEC', dec._degrees)


def separation_matrix(ra1, dec1, ra2, dec2, max_separation=None):
    """Build a matrix of pair-wise separation between (ra,dec) pointings.
    The ra1 and dec1 arrays must have the same shape. The ra2 and dec2 arrays
    must also have the same shape, but it can be different from the (ra1,dec1)
    shape, resulting in a non-square return matrix.
    Uses the Haversine formula for better accuracy at low separations. See
    https://en.wikipedia.org/wiki/Haversine_formula for details.
    Equivalent to using the separations() method of astropy.coordinates.ICRS,
    but faster since it bypasses any units.
    Parameters
    ----------
    ra1 : array
        1D array of n1 RA coordinates in degrees (without units attached).
    dec1 : array
        1D array of n1 DEC coordinates in degrees (without units attached).
    ra2 : array
        1D array of n2 RA coordinates in degrees (without units attached).
    dec2 : array
        1D array of n2 DEC coordinates in degrees (without units attached).
    max_separation : float or None
        When present, the matrix elements are replaced with booleans given
        by (value <= max_separation), which saves some computation.
    Returns
    -------
    array
        Array with shape (n1,n2) with element [i1,i2] giving the 3D separation
        angle between (ra1[i1],dec1[i1]) and (ra2[i2],dec2[i2]) in degrees
        or, if max_separation is not None, booleans (value <= max_separation).
    """
    ra1, ra2 = np.deg2rad(ra1), np.deg2rad(ra2)
    dec1, dec2 = np.deg2rad(dec1), np.deg2rad(dec2)
    if ra1.shape != dec1.shape:
        raise ValueError('Arrays ra1, dec1 must have the same shape.')
    if len(ra1.shape) != 1:
        raise ValueError('Arrays ra1, dec1 must be 1D.')
    if ra2.shape != dec2.shape:
        raise ValueError('Arrays ra2, dec2 must have the same shape.')
    if len(ra2.shape) != 1:
        raise ValueError('Arrays ra2, dec2 must be 1D.')
    havRA12 = 0.5 * (1 - np.cos(ra2 - ra1[:, np.newaxis]))
    havDEC12 = 0.5 * (1 - np.cos(dec2 - dec1[:, np.newaxis]))
    havPHI = havDEC12 + np.cos(dec1)[:, np.newaxis] * np.cos(dec2) * havRA12
    if max_separation is not None:
        # Replace n1 x n2 arccos calls with a single sin call.
        threshold = np.sin(0.5 * np.deg2rad(max_separation)) ** 2
        return havPHI <= threshold
    else:
        return np.rad2deg(np.arccos(np.clip(1 - 2 * havPHI, -1, +1)))


def get_visible(where, tbegin, exptime, ra0, dec0, fov, satellites, title, groups=None, nsteps=10, oversampling=5, location = None):
	"""
	"""
	if groups is not None:
		satellitess={}
		for name in groups:
			satellitess.update(get_satellites(satellites, name))
		satellites=satellitess
	cos_dec0 = np.cos(np.deg2rad(dec0))
	# Calculate the end time of this exposure.
	tend = ts.utc(tbegin.utc_datetime() + datetime.timedelta(seconds=exptime))
	# Calculate the exposure midpoint.
	tmid = ts.tai_jd(0.5 * (tbegin.tai + tend.tai))
	# Calculate equal coarse time steps covering [tbegin, tend].
	tcoarse = np.linspace(tbegin.tai, tend.tai, nsteps)
	rcoarse = np.empty((3, nsteps))
	tsteps = ts.tai_jd(tcoarse)
	# Build a fine grid of time steps during the exposure for interpolation.
	tfine = np.linspace(tbegin.tai, tend.tai, oversampling * nsteps)
	# Loop over satellites.
	names, ravec, decvec = [], [], []
	plt.figure(figsize=(10,10))
	for i, satellite in enumerate(satellites.values()):
		# Calculate satellite position relative to the observer at each time step.
		topocentric = (satellite - where).at(tsteps)
		# Calculate the corresponding astrometric (ra, dec)
		ra, dec, _ = topocentric.radec()
		# Check for propagation errors.
		if np.any(np.isnan(ra.radians) | np.isnan(dec.radians)):
			#print(f'Unable to propagate {satellite.name}')
			continue
		# Tabulate the satellite position on the finer grid using cubic interpolation.
		rcoarse[0] = np.sin(dec.radians)
		rcoarse[1] = np.cos(ra.radians)
		rcoarse[2] = np.sin(ra.radians)
			#'''Why is cos(dec) not used?'''
		interpolator = scipy.interpolate.interp1d(
			tcoarse, rcoarse, axis=1, kind='cubic', copy=False,
			bounds_error=False, assume_sorted=True)
		sin_dec, cos_ra, sin_ra = interpolator(tfine)
		dec_deg = np.rad2deg(np.arcsin(sin_dec))
		ra_deg = np.fmod(360 + np.rad2deg(np.arctan2(sin_ra, cos_ra)), 360)
		# Does this satellite ever enter the fov?
		#sep = separation_matrix([ra0], [dec0], ra_deg, dec_deg, max_separation=0.6 * fov).reshape(-1)
		sep = separation_matrix([ra0], [dec0], ra_deg, dec_deg).reshape(-1)
		if sep.min() < 2 * fov:
			plt.plot(ra_deg, dec_deg, '.-', label=satellite.name)
			#names.append(satellite.name)
			ravec.append(ra_deg)
			decvec.append(dec_deg)
		if groups is not None:
			names.append(satellite.name)
	plt.xlim(ra0 - 2 * fov / cos_dec0, ra0 + 2 * fov / cos_dec0)
	plt.ylim(dec0 - 2 * fov, dec0 + 2 * fov)
	plt.xlabel('RA [deg]')
	plt.ylabel('DEC [deg]')
	plt.title(title)
	ax = plt.gca()
	ax.invert_xaxis()
	ax.set_aspect(1 / cos_dec0)
	#plt.xlim(ra0+15, ra0-15, nsteps)
	#plt.ylim(dec0-15, dec0+15, nsteps)
	# Draw circle showing fov.
	theta = np.linspace(0, 2 * np.pi, 50)
	x = ra0 + 0.5 * fov * np.cos(theta) / cos_dec0
	y = dec0 + 0.5 * fov * np.sin(theta)
	ax.plot(x, y, ls=':', c='gray', label='FOV')
	if location is not None:
		plt.legend(loc=location)
	else:
		plt.legend(loc='best')
	plt.grid()
	plt.savefig(title)
	return names, ravec, decvec

def radec(satellite, where, tbegin, exptime):
	#satellite is an EarthSatellite object from skyfield
	#where is the skyfield.Topos object with lat, lon, and alt
	#tbegin is the beginning time to calculate position
	#exptime is how long to calculate for in seconds
	ra, dec = [], []
	time_delta = np.linspace(0, exptime, exptime)
	for second in time_delta:
		time = ts.utc(tbegin.utc_datetime() + datetime.timedelta(seconds=second))
		difference = satellite - where
		topocentric = difference.at(time)
		ra2, dec2, d2 = topocentric.radec() #ICRF COORDS SINCE EPOCH NOT SPECIFIED
		ra.append(ra2._degrees)
		dec.append(dec2._degrees)
	return ra, dec
