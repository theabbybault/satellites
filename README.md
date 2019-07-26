# satellites

the purpose of this is to track satellites using tle files from celestrak. 

The functions and their descriptions are below.

load_satellites(reload, groups, debris): 
Loads satellites from the tle files at celestrak.
Reload is defaulted to 'True'. Groups is defaulted to 'None' which will pull every satellite from the tle files at celestrak. To specify groups, enter the names of the files of the particular groups wanted as a tuple of strings (see examples below). There is no need for file extensions here, and any number of groups can be specified. Debris is defaulted to 'True' which includes the debris tle files from celestrak. If debris is set to 'False', debris will be excluded. If specifying groups, there is no need to set debris = True or False. Examples of different loading options below:
All satellites including debris: load_satellites(reload = False) There is no need to specify groups or debris when you want all satellites.
All satellites excluding debris: load_satellites(reload = False, debris = False)
Specifying specific satellites: load_satellites(reload = False, groups = ('group1', 'group2', 'group3'))

get_satellites(satellites, string):
Gets the satellite data for the given satellites. Data includes satellite name, type?, target_name, number, and epoch for any satellite with 'string' in its name. Satellites is the dictionary of satellites that are loaded with load_satellites.

from_localtime(^args, tzoffset = -7):
Gets the utc time given your local time. For Irvine, CA, the offset is 7 hours.

get_pos(where, when, string):
Gets the next position of satellites with 'string' in their name. Prints the satellite name, altitude, azimuth, RA, and DEC in degrees.

separation_matrix(ra1, dec1, ra2, dec2, max_separation = None):
Builds a matrix of pair-wise separation between (ra,dec) pointings. This was used in the next function.

get_visible(where, tbegin, exptime, ra0, dec0, fov, satellites, title, string, nsteps, oversampling, location):
Returns the names, RA, and DEC of satellites visible during a specified viewing time over a telescope, as well as a plot of their position. 
where is the location of the telescope. tbegin is the time object of the start of the exposure (use from_localtime since the time needs to be in utc time). exptime is the length of the exposure in seconds. ra0 and dec0 are the ra and dec the telescope is pointing to. fov is the field of view of the telescope. satellites are the dictionary of satellites defined using load_satellites. title is the title of the plot. string is to specify what satellites you would like to see (if they are over the telescope at that time, for example, if you would like to see the ISS then you would put string = ['ISS']. nsteps is the number of time steps you would like to have, default is 10. oversampling is default to 5. location is the location of the legend in the plot, default is None which will give the legend location as 'best'. 
