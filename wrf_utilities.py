#--------  Useful python functions for WRF data analysis -------
#
# Louis Marelle, 2022/05/12
#

#--------
def get_wrf_proj(WRF_FILENAME):
  """Returns WRF grid projection info wrf_proj for WRF file WRF_FILENAME"""
  # Imports
  from netCDF4 import Dataset
  # Define custom wrf_projection class
  class wrf_projection():
    def __init__(self, map_proj, imax, kmax, dx, dy, moad_cen_lat, truelat1,
        truelat2, stdlon, ref_lat, ref_lon, ref_i, ref_j, pole_lat, pole_lon,
        hemi, mminlu):
      self.map_proj = map_proj
      self.imax = imax
      self.jmax = jmax
      self.kmax = kmax
      self.dx = dx
      self.dy = dy
      self.moad_cen_lat = moad_cen_lat
      self.truelat1 = truelat1
      self.truelat2 = truelat2
      self.stdlon = stdlon
      self.ref_lat = ref_lat
      self.ref_lon = ref_lon
      self.ref_i = ref_i
      self.ref_j = ref_j
      self.pole_lat = pole_lat
      self.pole_lon = pole_lon
      self.hemi = hemi
      self.mminlu = mminlu
  # Open NetCDF and get WRF projection attribute values
  ncfile = Dataset(WRF_FILENAME)
  map_proj = ncfile.__getattribute__('MAP_PROJ')
  imax = ncfile.__getattribute__('WEST-EAST_GRID_DIMENSION') - 1
  jmax = ncfile.__getattribute__('SOUTH-NORTH_GRID_DIMENSION') - 1
  kmax = ncfile.__getattribute__('BOTTOM-TOP_GRID_DIMENSION') - 1
  dx = ncfile.__getattribute__('DX')/1000.
  dy = ncfile.__getattribute__('DY')/1000.
  moad_cen_lat = ncfile.__getattribute__('MOAD_CEN_LAT')
  truelat1 = ncfile.__getattribute__('TRUELAT1')
  truelat2 = ncfile.__getattribute__('TRUELAT2')
  stdlon = ncfile.__getattribute__('STAND_LON')
  ref_lat = ncfile.__getattribute__('CEN_LAT')
  ref_lon = ncfile.__getattribute__('CEN_LON')
  pole_lat = ncfile.__getattribute__('POLE_LAT')
  pole_lon = ncfile.__getattribute__('POLE_LON')
  mminlu = ncfile.__getattribute__('MMINLU')
  ncfile.close()
  ref_i = float((imax + 1.0) / 2.0)
  ref_j = float((jmax + 1.0) / 2.0)
  if(truelat1 > 0):
    hemi = 1.0
  else:
    hemi = -1.0
  # Create wrf_projection instance wrf_proj
  wrf_proj = wrf_projection(map_proj, imax, kmax, dx, dy, moad_cen_lat, 
    truelat1, truelat2, stdlon, ref_lat, ref_lon, ref_i, ref_j, pole_lat,
    pole_lon, hemi, mminlu)
  # Return wrf_proj
  return wrf_proj


def wrf_llij(lat, lon, wrf_proj):
  """Convert lat and lon into WRF i,j for the WRF grid defined in wrf_proj"""
  import numpy as np
  import math
  wrfi = []
  wrfj = []
  rad_per_deg =  math.pi/180.
  # Earth radius in kilometers divided by dx
  rebydx = 6370./wrf_proj.dx

  # Convert lists to numpy arrays
  lat = np.asarray(lat)
  lon = np.asarray(lon)

  if(wrf_proj.map_proj == 1):
  #-------- Lambert conformal projection --------
  # Subroutines set_lc, lc_cone and llij_lc from WPS: geogrid/src/module_map_utils.f90

    #---- Subroutine to compute the cone factor of a Lambert Conformal projection

    # Input Args
    #REAL, INTENT(IN)             :: wrf_proj.truelat1  # (-90 -> 90 degrees N)
    #REAL, INTENT(IN)             :: wrf_proj.truelat2  #   "   "  "   "     "
    # Output Args
    #REAL, INTENT(OUT)            :: cone

    # First, see if this is a secant or tangent projection.  For tangent
    # projections, wrf_proj.truelat1 = wrf_proj.truelat2 and the cone is tangent to the 
    # Earth's surface at this latitude.  For secant projections, the cone
    # intersects the Earth's surface at each of the distinctly different
    # latitudes
    if(np.abs(wrf_proj.truelat1 - wrf_proj.truelat2) > 0.1):
       cone = np.log10(np.cos(wrf_proj.truelat1 * rad_per_deg)) - \
              np.log10(np.cos(wrf_proj.truelat2 * rad_per_deg))
       cone = cone / (np.log10(np.tan((45.0 - np.abs(wrf_proj.truelat1) / 2.0) * rad_per_deg)) - \
              np.log10(np.tan((45.0 - np.abs(wrf_proj.truelat2) / 2.0) * rad_per_deg)))
    else:
       cone = np.sin(np.abs(wrf_proj.truelat1)*rad_per_deg )

    #---- Initialize the remaining items in the proj structure for a
    #---- lambert conformal grid.

    # Compute longitude differences and ensure we stay out of the
    # forbidden "cut zone"
    deltalon1 = wrf_proj.ref_lon - wrf_proj.stdlon
    if(deltalon1 > +180.):
      deltalon1 = deltalon1 - 360.
    if(deltalon1 < -180.):
      deltalon1 = deltalon1 + 360.

    # Convert wrf_proj.truelat1 to radian and compute COS for later use
    tl1r = wrf_proj.truelat1 * rad_per_deg
    ctl1r = np.cos(tl1r)

    # Compute the radius to our known lower-left (sw) corner
    rsw = rebydx * ctl1r/cone * \
           (np.tan((90.0*wrf_proj.hemi-wrf_proj.ref_lat)*rad_per_deg/2.) / \
            np.tan((90.0*wrf_proj.hemi-wrf_proj.truelat1)*rad_per_deg/2.))** cone

    # Find pole point
    arg = cone*(deltalon1*rad_per_deg)
    polei = wrf_proj.hemi*wrf_proj.ref_i - wrf_proj.hemi * rsw * np.sin(arg)
    polej = wrf_proj.hemi*wrf_proj.ref_j + rsw * np.cos(arg)


    #---- Subroutine to compute the geographical latitude and longitude values
    #---- to the cartesian x/y on a Lambert Conformal projection.

    # Input Args
    #REAL, INTENT(IN)              :: lat      # Latitude (-90->90 deg N)
    #REAL, INTENT(IN)              :: lon      # Longitude (-180->180 E)

    # Output Args                 
    #REAL, INTENT(OUT)             ::wrfi       # Cartesian X coordinate
    #REAL, INTENT(OUT)             :: wrfj        # Cartesian Y coordinate

    #---- Compute deltalon between known longitude and standard lon and ensure
    # it is not in the cut zone
    deltalon = lon - wrf_proj.stdlon
    if np.any(deltalon > 180.):
      deltalon[deltalon > 180.] = deltalon[deltalon > 180.] - 360.
    if np.any(deltalon < -180.):
      deltalon[deltalon < -180.] = deltalon[deltalon < -180.] + 360.

    # Convert wrf_proj.truelat1 to radian and compute COS for later use
    tl1r = wrf_proj.truelat1 * rad_per_deg
    ctl1r = np.cos(tl1r)

    # Radius to desired point
    rm = rebydx * ctl1r/cone * \
         (np.tan((90.0*wrf_proj.hemi-lat)*rad_per_deg/2.) / \
          np.tan((90.0*wrf_proj.hemi-wrf_proj.truelat1)*rad_per_deg/2.)) ** cone

    arg = cone*(deltalon*rad_per_deg)
    wrfi = polei + wrf_proj.hemi * rm * np.sin(arg)
    wrfj = polej - rm * np.cos(arg)

    # Finally, if we are in the southern hemisphere, flip the i/j
    # values to a coordinate system where (1,1) is the SW corner
    # (what we assume) which is different than the original NCEP
    # algorithms which used the NE corner as the origin in the 
    # southern hemisphere (left-hand vs. right-hand coordinate?)
    wrfi = wrf_proj.hemi * wrfi
    wrfj = wrf_proj.hemi * wrfj

  elif(wrf_proj.map_proj == 2):
    #-------- Polar stereographic projection --------
    # Compute the reference longitude by rotating 90 degrees to the east
    # to find the longitude line parallel to the positive x-axis.
    reflon = wrf_proj.stdlon + 90.0
    # Compute numerator term of map scale factor
    scale_top = 1. + wrf_proj.hemi * np.sin(wrf_proj.truelat1 * rad_per_deg)
    # Compute radius to lower-left (SW) corner
    ala1 = wrf_proj.ref_lat * rad_per_deg
    rsw = rebydx * np.cos(ala1) * scale_top/(1.+wrf_proj.hemi*np.sin(ala1))
    # Find the pole point
    alo1 = (wrf_proj.ref_lon - reflon) * rad_per_deg
    polei = wrf_proj.ref_i - rsw * np.cos(alo1)
    polej = wrf_proj.ref_j - wrf_proj.hemi * rsw * np.sin(alo1)
    # Find radius to desired point
    ala = lat * rad_per_deg
    rm = rebydx * np.cos(ala) * scale_top/(1. + wrf_proj.hemi*np.sin(ala))
    alo = (lon - reflon) * rad_per_deg
    wrfi = polei + rm * np.cos(alo)
    wrfj = polej + wrf_proj.hemi * rm * np.sin(alo)

  else:
    #TODO throw error if projection wrf_proj.map_proj is not included
    print('wrf_llij: wrf_proj.map_proj '+str(wrf_proj.map_proj)+' value not handled')

  return wrfi, wrfj


def wrf_ijll(wrfi, wrfj, wrf_proj):
  """Convert WRF i,j into lat and lon for the WRF grid defined in wrf_proj"""
  #TODO convert my matlab routine to python, starting with Lambert and Polar projections
  import numpy as np
  import math

  if(np.shape(wrfi) != np.shape(wrfj)):
    print('Error, np.shape(wrfi) != np.shape(wrfj)')
    exit

  # Convert lists to numpy arrays
  if np.shape(wrfi == ()):
    wrfi = [wrfi]
    wrfj = [wrfj]
  wrfi = np.asarray(wrfi)
  wrfj = np.asarray(wrfj)
  wrflat = np.empty((np.shape(wrfi)))
  wrflon = np.empty((np.shape(wrfj)))
  if np.shape(wrflat == ()):
    wrflat = [wrflat]
    wrflon = [wrflon]
    wrflat = np.asarray(wrflat)
    wrflon = np.asarray(wrflon)
  wrflat[:] = np.nan
  wrflon[:] = np.nan

  # Earth radius in kilometers divided by dx
  rebydx = 6370./wrf_proj.dx
  deg_per_rad =  180.0/math.pi
  rad_per_deg =  math.pi/180.

  # To know what each map_proj corresponds to, check in WPS code
  # geogrid/src/misc_definitions_module.f90
  if(wrf_proj.map_proj == 1):
    #-------- Lambert conformal projection --------
    # Subroutines from WPS: geogrid/src/module_map_utils.f90
    chi1 = (90.0 - wrf_proj.hemi * wrf_proj.truelat1) * rad_per_deg
    chi2 = (90.0 - wrf_proj.hemi * wrf_proj.truelat2) * rad_per_deg
    
    # See if we are in the southern hemispere and flip the indices if we are. 
    inew = wrf_proj.hemi * wrfi
    jnew = wrf_proj.hemi * wrfj
        
    # Compute radius**2 to i/j location 
    reflon = wrf_proj.stdlon + 90.0
    ala1 = wrf_proj.ref_lat * rad_per_deg
    alo1 = (wrf_proj.ref_lon - reflon) * rad_per_deg
    scale_top = 1.0 + wrf_proj.hemi * np.sin(wrf_proj.truelat1 * rad_per_deg)

    # Compute longitude differences and ensure we stay out of the
    # forbidden "cut zone"
    deltalon1 = wrf_proj.ref_lon - wrf_proj.stdlon
    if(deltalon1 > +180.):
      deltalon1 = deltalon1 - 360.
    if(deltalon1 < -180.):
      deltalon1 = deltalon1 + 360.

    # First, see if this is a secant or tangent projection.  For tangent
    # projections, wrf_proj.truelat1 = wrf_proj.truelat2 and the cone is tangent to the 
    # Earth's surface at this latitude.  For secant projections, the cone
    # intersects the Earth's surface at each of the distinctly different
    # latitudes
    if(np.abs(wrf_proj.truelat1-wrf_proj.truelat2) > 0.1):
      cone = np.log10(np.cos(wrf_proj.truelat1*rad_per_deg)) \
             - np.log10(np.cos(wrf_proj.truelat2*rad_per_deg))
      cone = cone / (np.log10(np.tan((45.0 - np.abs(wrf_proj.truelat1) / 2.0) * rad_per_deg)) \
             - np.log10(np.tan((45.0 - np.abs(wrf_proj.truelat2) / 2.0) * rad_per_deg)))
    else:
      cone = np.sin(np.abs(wrf_proj.truelat1) * rad_per_deg )

    # Convert wrf_proj.truelat1 to radian and compute COS for later use
    tl1r = wrf_proj.truelat1 * rad_per_deg
    ctl1r = np.cos(tl1r)

    # Compute the radius to our known lower-left (sw) corner
    rsw = rebydx * ctl1r/cone \
          * (np.tan((90.0*wrf_proj.hemi-wrf_proj.ref_lat)*rad_per_deg/2.0) \
          / (np.tan((90.0*wrf_proj.hemi-wrf_proj.truelat1)*rad_per_deg/2.0)))**cone

    # Find pole point
    arg = cone*(deltalon1*rad_per_deg)
    polei = wrf_proj.hemi*wrf_proj.ref_i - wrf_proj.hemi*rsw*np.sin(arg)
    polej = wrf_proj.hemi*wrf_proj.ref_j + rsw*np.cos(arg)

    xx = inew - polei
    yy = polej - jnew
    r2 = (xx*xx + yy*yy)
    r = np.sqrt(r2)/rebydx

    # Convert to lat/lon
    wrflat[r2==0.0] = wrf_proj.hemi * 90.0
    wrflon[r2==0.0] = wrf_proj.stdlon

    wrflon[r2!=0.0] = wrf_proj.stdlon + deg_per_rad \
                      * np.arctan(wrf_proj.hemi * xx[r2!=0.0]/yy[r2!=0.0]) / cone
    wrflon[r2!=0.0] = np.mod(wrflon[r2!=0.0] + 360.0, 360.0)
    
    chi = np.empty((np.shape(wrfi)))
    chi[:] = np.nan
    if (chi1==chi2) :
      chi[r2!=0.0] = 2.0 * np.arctan( (r[r2!=0.0]/np.tan(chi1))**(1.0/cone) \
                     * np.tan(chi1*0.5) )
    else:
      chi[r2!=0.0] = 2.0 * np.arctan( (r[r2!=0.0]*cone/np.sin(chi1))**(1.0/cone) \
                     * np.tan(chi1*0.5))

    wrflat[r2!=0.0] = (90.0-chi[r2!=0.0]*deg_per_rad)*wrf_proj.hemi


  elif(wrf_proj.map_proj == 2):
    #-------- Polar stereographic projection --------
    # Compute the reference longitude by rotating 90 degrees to the east
    # to find the longitude line parallel to the positive x-axis.
    reflon = wrf_proj.stdlon + 90.0
    # Compute numerator term of map scale factor
    scale_top = 1.0 + wrf_proj.hemi * np.sin(wrf_proj.truelat1 * rad_per_deg)
    # Calculate pole i and j
    # Compute radius to lower-left (SW) corner
    ala1 = wrf_proj.ref_lat * rad_per_deg
    rsw = rebydx*np.cos(ala1)*scale_top / (1.0 + wrf_proj.hemi*np.sin(ala1))
    # Find the pole point
    alo1 = (wrf_proj.ref_lon - reflon) * rad_per_deg
    polei = wrf_proj.ref_i - rsw*np.cos(alo1)
    polej = wrf_proj.ref_j - wrf_proj.hemi * rsw * np.sin(alo1)
    # Compute radius to point of interest
    xx = wrfi - polei
    yy = (wrfj - polej) * wrf_proj.hemi
    r2 = xx**2.0 + yy**2.0
    # Now the magic code
    wrflat[r2==0] = wrf_proj.hemi * 90.0
    wrflon[r2==0] = reflon
    gi2 = (rebydx * scale_top) ** 2.0
    wrflat[r2!=0] = deg_per_rad * wrf_proj.hemi \
                    * np.arcsin( (gi2-r2[r2!=0]) / (gi2+r2[r2!=0]) )
    arccos = np.arccos(xx / np.sqrt(r2))
    if(np.any(np.logical_and(r2!=0.0, yy>0.0))):
        wrflon[np.logical_and(r2!=0.0, yy>0.0)] = reflon + deg_per_rad \
                                        * arccos[np.logical_and(r2!=0.0, yy>0.0)];
    if(np.any(np.logical_and(r2!=0.0, yy<=0.0))):
        wrflon[np.logical_and(r2!=0.0, yy<=0.0)] = reflon - deg_per_rad \
                                        * arccos[np.logical_and(r2!=0.0, yy<=0.0)];

  elif(wrf_proj.map_proj == 3):
    #-------- Mercator projection --------
    # From module_map_utils.f90
    clain = np.cos(rad_per_deg * wrf_proj.truelat1)
    dlon = 1.0 / (rebydx * clain)
    #! Compute distance from equator to origin, and store in the 
    #! proj%rsw tag.
    rsw = 0.
    if(wrf_proj.ref_lat != 0.):
       rsw = (np.log10(np.tan(0.5*((wrf_proj.ref_lat+90.)*rad_per_deg))))/dlon
    wrflat = 2.0 * np.arctan(np.exp(dlon * (rsw + wrfj - wrf_proj.ref_j))) \
             * deg_per_rad - 90.
    wrflon = (wrfi - wrf_proj.ref_i) * dlon * deg_per_rad + wrf_proj.ref_lon

  else:
    #TODO throw error if projection wrf_proj.map_proj is not included
    print('wrf_proj.map_proj '+str(wrf_proj.map_proj)+' value not handled')

  # Convert to a -180 -> 180 East convention
  if np.any(wrflon > 180.):
    wrflon[wrflon > 180.0] = wrflon[wrflon > 180.0] - 360.0
  if np.any(wrflon > 180.):
    wrflon[wrflon < -180.0] = wrflon[wrflon < -180.0] + 360.0

  return wrflat, wrflon


def wrf_interp3(wrf_x, wrf_y, wrf_z, wrf_var, x, y, z):
  """Interpolate WRF variable WRF var at coordinates x,y,z"""
  #TODO try not to reimport this each time, try putting at the beginning of the
  #     module and time, or find a way to import globally
  #TODO try to do everything in vectorial calculations to speed it up a bit?
  #TODO maybe remove wrf_x, wrf_y, since they are implicit in wrf_var
  #TODO or recode as wrf_xlat, wrf_xlong, wrf_z
  #TODO check that x, y, z have the same shapes, check that wrf_var has the
  #     shape of wrf_x, wrf_y, wrf_z
  import numpy as np
  # Initialize
  wrf_x = np.asarray(wrf_x)
  wrf_y = np.asarray(wrf_y)
  wrf_z = np.asarray(wrf_z)
  x = np.asarray(x).astype(float)
  y = np.asarray(y).astype(float)
  z = np.asarray(z).astype(float)
  z[z < min(wrf_z)] = min(wrf_z)
  wrf_var_interp = np.empty(x.shape)
  wrf_x_weights_below = np.empty(x.shape)
  wrf_x_weights_above = np.empty(x.shape)
  wrf_y_weights_below = np.empty(y.shape)
  wrf_y_weights_above = np.empty(y.shape)
  wrf_z_weights_below = np.empty(z.shape)
  wrf_z_weights_above = np.empty(z.shape)
  wrf_x_below = np.empty(x.shape)
  wrf_x_above = np.empty(x.shape)
  wrf_y_below = np.empty(y.shape)
  wrf_y_above = np.empty(y.shape)
  wrf_z_below = np.empty(z.shape)
  wrf_z_above = np.empty(z.shape)
  # wrf_var_interp[:] = np.nan
  # wrf_x_weights_below[:] = np.nan
  # wrf_x_weights_above[:] = np.nan
  # wrf_y_weights_below[:] = np.nan
  # wrf_y_weights_above[:] = np.nan
  # wrf_z_weights_below[:] = np.nan
  # wrf_z_weights_above[:] = np.nan
  # wrf_x_below[:] = np.nan
  # wrf_x_above[:] = np.nan
  # wrf_y_below[:] = np.nan
  # wrf_y_above[:] = np.nan
  # wrf_z_below[:] = np.nan
  # wrf_z_above[:] = np.nan
  wrf_i_below = np.empty(x.shape)
  wrf_i_above = np.empty(x.shape)
  wrf_j_below = np.empty(y.shape)
  wrf_j_above = np.empty(y.shape)
  wrf_k_below = np.empty(z.shape)
  wrf_k_above = np.empty(z.shape)
  # wrf_i_below[:] = np.nan
  # wrf_i_above[:] = np.nan
  # wrf_j_below[:] = np.nan
  # wrf_j_above[:] = np.nan
  # wrf_k_below[:] = np.nan
  # wrf_k_above[:] = np.nan
  # Ignore out of bounds values
  outofbounds_mask = (x>max(wrf_x)) | (x<min(wrf_x)) | \
                     (y>max(wrf_y)) | (y<min(wrf_y)) | \
                     (z>max(wrf_z)) | (z<min(wrf_z))
  x[outofbounds_mask] = np.nan
  y[outofbounds_mask] = np.nan
  z[outofbounds_mask] = np.nan
  inbounds_mask = np.logical_not(outofbounds_mask)
  #----Trilinear interpolation
  # Find the closest i,j,k values in wrf_x wrf_y, wrf_z
  wrf_i_above[inbounds_mask] = np.searchsorted(wrf_x, x[inbounds_mask], side="left")
  wrf_j_above[inbounds_mask] = np.searchsorted(wrf_y, y[inbounds_mask], side="left")
  wrf_k_above[inbounds_mask] = np.searchsorted(wrf_z, z[inbounds_mask], side="left")
  wrf_i_above = wrf_i_above.astype(int)
  wrf_j_above = wrf_j_above.astype(int)
  wrf_k_above = wrf_k_above.astype(int)
  wrf_i_below = wrf_i_above - 1 
  wrf_j_below = wrf_j_above - 1
  wrf_k_below = wrf_k_above - 1
  #-- Calculate weights and below/above values
  # General case: calculate weights in x,y,z
  wrf_x_below = np.asarray(wrf_x[wrf_i_below])
  wrf_x_above = np.asarray(wrf_x[wrf_i_above])
  wrf_y_below = np.asarray(wrf_y[wrf_j_below])
  wrf_y_above = np.asarray(wrf_y[wrf_j_above])
  wrf_z_below = np.asarray(wrf_z[wrf_k_below])
  wrf_z_above = np.asarray(wrf_z[wrf_k_above])
  # Use mask to avoid division by 0
  mask_x = (np.around(x) != x)
  mask_y = (np.around(y) != y)
  mask_z = (np.around(z) != z)
  wrf_x_weights_below[mask_x] = np.asarray((wrf_x_above[mask_x] - x[mask_x])/(wrf_x_above[mask_x] - wrf_x_below[mask_x]))
  wrf_x_weights_above = np.asarray(1.0 - wrf_x_weights_below)
  wrf_y_weights_below[mask_y] = np.asarray((wrf_y_above[mask_y] - y[mask_y])/(wrf_y_above[mask_y] - wrf_y_below[mask_y]))
  wrf_y_weights_above = np.asarray(1.0 - wrf_y_weights_below)
  wrf_z_weights_below[mask_z] = np.asarray((wrf_z_above[mask_z] - z[mask_z])/(wrf_z_above[mask_z] - wrf_z_below[mask_z]))
  wrf_z_weights_above = np.asarray(1.0 - wrf_z_weights_below)
  # Special cases for exact values: use below=above and weights = [0, 1]
  wrf_x_weights_below[np.around(x) == x] = 0.0
  wrf_x_weights_above[np.around(x) == x] = 1.0
  wrf_y_weights_below[np.around(y) == y] = 0.0
  wrf_y_weights_above[np.around(y) == y] = 1.0
  wrf_z_weights_below[np.around(z) == z] = 0.0
  wrf_z_weights_above[np.around(z) == z] = 1.0
  wrf_x_below[np.around(x) == x] = x[np.around(x) == x]
  wrf_x_above[np.around(x) == x] = x[np.around(x) == x]
  wrf_y_below[np.around(y) == y] = y[np.around(y) == y]
  wrf_y_above[np.around(y) == y] = y[np.around(y) == y]
  wrf_z_below[np.around(z) == z] = z[np.around(z) == z]
  wrf_z_above[np.around(z) == z] = z[np.around(z) == z]
  #-- Trilinear interpolation
  # Get values for interpolation
  wrf_var_interp_xb_yb_zb = wrf_var[wrf_k_below.astype(int), wrf_j_below.astype(int), wrf_i_below.astype(int)]
  wrf_var_interp_xa_yb_zb = wrf_var[wrf_k_below.astype(int), wrf_j_below.astype(int), wrf_i_above.astype(int)]
  wrf_var_interp_xb_ya_zb = wrf_var[wrf_k_below.astype(int), wrf_j_above.astype(int), wrf_i_below.astype(int)]
  wrf_var_interp_xa_ya_zb = wrf_var[wrf_k_below.astype(int), wrf_j_above.astype(int), wrf_i_above.astype(int)]
  wrf_var_interp_xb_yb_za = wrf_var[wrf_k_above.astype(int), wrf_j_below.astype(int), wrf_i_below.astype(int)]
  wrf_var_interp_xa_yb_za = wrf_var[wrf_k_above.astype(int), wrf_j_below.astype(int), wrf_i_above.astype(int)]
  wrf_var_interp_xb_ya_za = wrf_var[wrf_k_above.astype(int), wrf_j_above.astype(int), wrf_i_below.astype(int)]
  wrf_var_interp_xa_ya_za = wrf_var[wrf_k_above.astype(int), wrf_j_above.astype(int), wrf_i_above.astype(int)]
  # Interpolate in x
  wrf_var_interp_x_yb_zb = wrf_x_weights_below * wrf_var_interp_xb_yb_zb + wrf_x_weights_above * wrf_var_interp_xa_yb_zb
  wrf_var_interp_x_ya_zb = wrf_x_weights_below * wrf_var_interp_xb_ya_zb + wrf_x_weights_above * wrf_var_interp_xa_ya_zb
  wrf_var_interp_x_yb_za = wrf_x_weights_below * wrf_var_interp_xb_yb_za + wrf_x_weights_above * wrf_var_interp_xa_yb_za
  wrf_var_interp_x_ya_za = wrf_x_weights_below * wrf_var_interp_xb_ya_za + wrf_x_weights_above * wrf_var_interp_xa_ya_za
  # Interpolate in y
  wrf_var_interp_x_y_zb = wrf_y_weights_below * wrf_var_interp_x_yb_zb + wrf_y_weights_above * wrf_var_interp_x_ya_zb
  wrf_var_interp_x_y_za = wrf_y_weights_below * wrf_var_interp_x_yb_za + wrf_y_weights_above * wrf_var_interp_x_ya_za
  # Interpolate in z
  wrf_var_interp = wrf_z_weights_below * wrf_var_interp_x_y_zb + wrf_z_weights_above * wrf_var_interp_x_y_za

  return wrf_var_interp

def wrf_interp2(model_x, model_y, model_var, x, y):
  """Interpolate model variable model_var at coordinates x,y"""
  #TODO try not to reimport this each time, try putting at the beginning of the
  #     module and time, or find a way to import globally
  #TODO check that x, y, have the same shapes, check that model_var has the
  #     shape of model_x, model_y
  import numpy as np
  # Initialize
  model_x = np.asarray(model_x)
  model_y = np.asarray(model_y)
  x = np.asarray(x).astype(float)
  y = np.asarray(y).astype(float)
  model_var_interp = np.empty(x.shape)
  model_x_weights_below = np.empty(x.shape)
  model_x_weights_above = np.empty(x.shape)
  model_y_weights_below = np.empty(y.shape)
  model_y_weights_above = np.empty(y.shape)
  model_x_below = np.empty(x.shape)
  model_x_above = np.empty(x.shape)
  model_y_below = np.empty(y.shape)
  model_y_above = np.empty(y.shape)
  model_i_below = np.empty(x.shape)
  model_i_above = np.empty(x.shape)
  model_j_below = np.empty(y.shape)
  model_j_above = np.empty(y.shape)
  # Ignore out of bounds values
  outofbounds_mask = (x>max(model_x)) | (x<min(model_x)) | \
                     (y>max(model_y)) | (y<min(model_y))
  x[outofbounds_mask] = np.nan
  y[outofbounds_mask] = np.nan
  inbounds_mask = np.logical_not(outofbounds_mask)
  #----Trilinear interpolation
  # Find the closest i,j,k values in model_x model_y
  model_i_above[inbounds_mask] = np.searchsorted(model_x, x[inbounds_mask], side="left")
  model_j_above[inbounds_mask] = np.searchsorted(model_y, y[inbounds_mask], side="left")
  model_i_above = model_i_above.astype(int)
  model_j_above = model_j_above.astype(int)
  model_i_below = model_i_above - 1 
  model_j_below = model_j_above - 1
  #-- Calculate weights and below/above values
  # General case: calculate weights in x,y
  model_x_below = np.asarray(model_x[model_i_below])
  model_x_above = np.asarray(model_x[model_i_above])
  model_y_below = np.asarray(model_y[model_j_below])
  model_y_above = np.asarray(model_y[model_j_above])
  # Use mask to avoid division by 0
  mask_x = (np.around(x) != x)
  mask_y = (np.around(y) != y)
  model_x_weights_below[mask_x] = np.asarray((model_x_above[mask_x] - x[mask_x])/(model_x_above[mask_x] - model_x_below[mask_x]))
  model_x_weights_above = np.asarray(1.0 - model_x_weights_below)
  model_y_weights_below[mask_y] = np.asarray((model_y_above[mask_y] - y[mask_y])/(model_y_above[mask_y] - model_y_below[mask_y]))
  model_y_weights_above = np.asarray(1.0 - model_y_weights_below)
  # Special cases for exact values: use below=above and weights = [0, 1]
  model_x_weights_below[np.around(x) == x] = 0.0
  model_x_weights_above[np.around(x) == x] = 1.0
  model_y_weights_below[np.around(y) == y] = 0.0
  model_y_weights_above[np.around(y) == y] = 1.0
  model_x_below[np.around(x) == x] = x[np.around(x) == x]
  model_x_above[np.around(x) == x] = x[np.around(x) == x]
  model_y_below[np.around(y) == y] = y[np.around(y) == y]
  model_y_above[np.around(y) == y] = y[np.around(y) == y]
  #-- Bilinear interpolation
  # Get values for interpolation
  model_var_interp_xb_yb = model_var[model_j_below.astype(int), model_i_below.astype(int)]
  model_var_interp_xa_yb = model_var[model_j_below.astype(int), model_i_above.astype(int)]
  model_var_interp_xb_ya = model_var[model_j_above.astype(int), model_i_below.astype(int)]
  model_var_interp_xa_ya = model_var[model_j_above.astype(int), model_i_above.astype(int)]
  # Interpolate in x
  model_var_interp_x_yb = model_x_weights_below * model_var_interp_xb_yb + model_x_weights_above * model_var_interp_xa_yb
  model_var_interp_x_ya = model_x_weights_below * model_var_interp_xb_ya + model_x_weights_above * model_var_interp_xa_ya
  # Interpolate in y
  model_var_interp = model_y_weights_below * model_var_interp_x_yb + model_y_weights_above * model_var_interp_x_ya

  return model_var_interp

def wrf_interp2_field(wrf_x, wrf_y, wrf_var, x, y):
  """Interpolate WRF variable WRF var at coordinates x,y,z"""
  #TODO try not to reimport this each time, try putting at the beginning of the
  #     module and time, or find a way to import globally
  #TODO try to do everything in vectorial calculations to speed it up a bit?
  #TODO maybe remove wrf_x, wrf_y, since they are implicit in wrf_var
  #TODO or recode as wrf_xlat, wrf_xlong,
  #TODO check that x, y, have the same shapes, check that wrf_var has the
  #     shape of wrf_x, wrf_y
  import numpy as np
  # Initialize
  wrf_x = np.asarray(wrf_x)
  wrf_y = np.asarray(wrf_y)
  x = np.asarray(x).astype(float)
  y = np.asarray(y).astype(float)
  wrf_var_interp = np.empty(x.shape)
  wrf_x_weights_below = np.empty(x.shape)
  wrf_x_weights_above = np.empty(x.shape)
  wrf_y_weights_below = np.empty(y.shape)
  wrf_y_weights_above = np.empty(y.shape)
  wrf_x_below = np.empty(x.shape)
  wrf_x_above = np.empty(x.shape)
  wrf_y_below = np.empty(y.shape)
  wrf_y_above = np.empty(y.shape)
  wrf_i_below = np.empty(x.shape)
  wrf_i_above = np.empty(x.shape)
  wrf_j_below = np.empty(y.shape)
  wrf_j_above = np.empty(y.shape)
  # Ignore out of bounds values
  outofbounds_mask = (x>max(wrf_x)) | (x<min(wrf_x)) | \
                     (y>max(wrf_y)) | (y<min(wrf_y))
  x[outofbounds_mask] = np.nan
  y[outofbounds_mask] = np.nan
  inbounds_mask = np.logical_not(outofbounds_mask)
  #----Trilinear interpolation
  # Find the closest i,j,k values in wrf_x wrf_y
  wrf_i_above[inbounds_mask] = np.searchsorted(wrf_x, x[inbounds_mask], side="left")
  wrf_j_above[inbounds_mask] = np.searchsorted(wrf_y, y[inbounds_mask], side="left")
  wrf_i_above = wrf_i_above.astype(int)
  wrf_j_above = wrf_j_above.astype(int)
  wrf_i_below = wrf_i_above - 1 
  wrf_j_below = wrf_j_above - 1
  #-- Calculate weights and below/above values
  # General case: calculate weights in x,y
  wrf_x_below = np.asarray(wrf_x[wrf_i_below])
  wrf_x_above = np.asarray(wrf_x[wrf_i_above])
  wrf_y_below = np.asarray(wrf_y[wrf_j_below])
  wrf_y_above = np.asarray(wrf_y[wrf_j_above])
  # Use mask to avoid division by 0
  mask_x = (np.around(x) != x)
  mask_y = (np.around(y) != y)
  wrf_x_weights_below[mask_x] = np.asarray((wrf_x_above[mask_x] - x[mask_x])/(wrf_x_above[mask_x] - wrf_x_below[mask_x]))
  wrf_x_weights_above = np.asarray(1.0 - wrf_x_weights_below)
  wrf_y_weights_below[mask_y] = np.asarray((wrf_y_above[mask_y] - y[mask_y])/(wrf_y_above[mask_y] - wrf_y_below[mask_y]))
  wrf_y_weights_above = np.asarray(1.0 - wrf_y_weights_below)
  # Special cases for exact values: use below=above and weights = [0, 1]
  wrf_x_weights_below[np.around(x) == x] = 0.0
  wrf_x_weights_above[np.around(x) == x] = 1.0
  wrf_y_weights_below[np.around(y) == y] = 0.0
  wrf_y_weights_above[np.around(y) == y] = 1.0
  wrf_x_below[np.around(x) == x] = x[np.around(x) == x]
  wrf_x_above[np.around(x) == x] = x[np.around(x) == x]
  wrf_y_below[np.around(y) == y] = y[np.around(y) == y]
  wrf_y_above[np.around(y) == y] = y[np.around(y) == y]
  #-- Bilinear interpolation
  # Get values for interpolation
  wrf_var_interp_xb_yb = wrf_var[wrf_j_below.astype(int), wrf_i_below.astype(int)]
  wrf_var_interp_xa_yb = wrf_var[wrf_j_below.astype(int), wrf_i_above.astype(int)]
  wrf_var_interp_xb_ya = wrf_var[wrf_j_above.astype(int), wrf_i_below.astype(int)]
  wrf_var_interp_xa_ya = wrf_var[wrf_j_above.astype(int), wrf_i_above.astype(int)]
  # Interpolate in x
  wrf_var_interp_x_yb = wrf_x_weights_below * wrf_var_interp_xb_yb + wrf_x_weights_above * wrf_var_interp_xa_yb
  wrf_var_interp_x_ya = wrf_x_weights_below * wrf_var_interp_xb_ya + wrf_x_weights_above * wrf_var_interp_xa_ya
  # Interpolate in y
  wrf_var_interp = wrf_y_weights_below * wrf_var_interp_x_yb + wrf_y_weights_above * wrf_var_interp_x_ya

  return wrf_var_interp

def calc_wrf_grid_edges(wrf_proj):
  """Return the latitude and longitude of the grid cell edges for the WRF grid defined in wrf_proj"""
  import numpy as np
  wrf_i_edge = np.linspace(0.5, wrf_proj.imax+0.5, wrf_proj.imax+1)
  wrf_j_edge = np.linspace(0.5, wrf_proj.jmax+0.5, wrf_proj.jmax+1)
  wrf_lat_edge = np.empty((wrf_proj.jmax+1, wrf_proj.imax+1))
  wrf_lon_edge = np.empty((wrf_proj.jmax+1, wrf_proj.imax+1))
  for ii in range(0,wrf_proj.imax+1):
    wrf_lat_edge[:,ii], wrf_lon_edge[:,ii] = wrf_ijll(wrf_i_edge[ii]*np.ones((np.shape(wrf_j_edge))),
                                                      wrf_j_edge, wrf_proj)
  return wrf_lat_edge, wrf_lon_edge


