"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os
import numpy as np


#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata)

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

#####
    # Lower and upper edge of computational domain:
    #Long
    clawdata.lower[0] = 94.6
    clawdata.upper[0] = 95.5
    #Lat
    clawdata.lower[1] = 28.125             
    clawdata.upper[1] = 30.32

    # Number of grid cells: Coarsest grid
    clawdata.num_cells[0] = 89
    clawdata.num_cells[1] = 254

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 3

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk43688'  # File to use for restart data

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = 6
        clawdata.tfinal = 172800.0
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [30000,31200,32400,36000,34800,36000]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 1
        clawdata.output_t0 = True
        

    clawdata.output_format == 'ascii'      # 'ascii' or 'netcdf' 

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'none'  # could be list
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0



    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 1



    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.016

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.7    
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 50000000
   



    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'

    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at tfinal.
        pass

    elif clawdata.checkpt_style == 2:
        # Specify a list of checkpoint times.  
        clawdata.checkpt_times = [0.1,0.15]

    elif clawdata.checkpt_style == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5


#####
    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 5

    # List of refinement ratios at each level (length at least mxnest-1)
    amrdata.refinement_ratios_x = [4,2,2,2]
    amrdata.refinement_ratios_y = [4,2,2,2]
    amrdata.refinement_ratios_t = [4,2,2,2]


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center','capacity','yleft','center']


    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 3

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0


    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    # == setregions.data values ==
    regions = rundata.regiondata.regions
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    
#AMR zones
    #1 Lake
    rundata.regiondata.regions.append([4, 4, 0, 1.e10, 94.72, 94.925, 30.17, 30.32])

    #Lake through gorge
    #rundata.regiondata.regions.append([1, 2, 0, 1.e10, 94.6, 95.5, 29.48, 30.32])
 
#Gauges
    rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    
    #Tongmai
    rundata.gaugedata.gauges.append([60, 95.064181, 30.100367, 0, 1.e10])
    rundata.gaugedata.gauges.append([61, 95.064773, 30.100117, 0, 1.e10])
    rundata.gaugedata.gauges.append([62, 95.065364, 30.099866 , 0, 1.e10])
    rundata.gaugedata.gauges.append([63, 95.065956, 30.099615, 0, 1.e10])
    rundata.gaugedata.gauges.append([64, 95.066548, 30.099364, 0, 1.e10])
    rundata.gaugedata.gauges.append([65, 95.067140, 30.099114 , 0, 1.e10])
    rundata.gaugedata.gauges.append([66, 95.067731, 30.098863, 0, 1.e10])
    
    #Landslides    
    rundata.gaugedata.gauges.append([70, 94.9571, 30.1663, 0, 1.e10])
    rundata.gaugedata.gauges.append([71, 95.0137, 30.0354, 0, 1.e10])
    rundata.gaugedata.gauges.append([72, 95.0171, 30.0283, 0, 1.e10])
    rundata.gaugedata.gauges.append([73, 95.0313, 30.0071, 0, 1.e10])
    rundata.gaugedata.gauges.append([74, 95.0371, 30.0046, 0, 1.e10])
    rundata.gaugedata.gauges.append([75, 95.0437, 30.0037, 0, 1.e10])
    rundata.gaugedata.gauges.append([76, 95.0521, 29.9988, 0, 1.e10])
    rundata.gaugedata.gauges.append([77, 95.0579, 29.9946, 0, 1.e10])
    rundata.gaugedata.gauges.append([78, 95.0604, 29.9904, 0, 1.e10])
    rundata.gaugedata.gauges.append([79, 95.0587, 29.9788, 0, 1.e10])
    rundata.gaugedata.gauges.append([80, 95.0687, 29.9729, 0, 1.e10])
    rundata.gaugedata.gauges.append([81, 95.0775, 29.9717, 0, 1.e10])
    rundata.gaugedata.gauges.append([82, 95.0796, 29.9696, 0, 1.e10])
    rundata.gaugedata.gauges.append([83, 95.0825, 29.9566, 0, 1.e10])
    rundata.gaugedata.gauges.append([84, 95.0804, 29.9529, 0, 1.e10])
    rundata.gaugedata.gauges.append([85, 95.0871, 29.9487, 0, 1.e10])
    rundata.gaugedata.gauges.append([86, 95.0892, 29.9483, 0, 1.e10])
    rundata.gaugedata.gauges.append([87, 95.0938, 29.9429, 0, 1.e10])
    rundata.gaugedata.gauges.append([88, 95.0971, 29.9379, 0, 1.e10])
    rundata.gaugedata.gauges.append([89, 95.1021, 29.9312, 0, 1.e10])
    rundata.gaugedata.gauges.append([90, 95.1325, 29.9108, 0, 1.e10])
    rundata.gaugedata.gauges.append([91, 95.1296, 29.8921, 0, 1.e10])
    rundata.gaugedata.gauges.append([92, 95.1454, 29.8754, 0, 1.e10])
    rundata.gaugedata.gauges.append([93, 95.1471, 29.8888, 0, 1.e10])
    rundata.gaugedata.gauges.append([94, 95.1787, 29.9013, 0, 1.e10])
    rundata.gaugedata.gauges.append([95, 95.1829, 29.8988, 0, 1.e10])
    rundata.gaugedata.gauges.append([96, 95.1921, 29.8904, 0, 1.e10])
    rundata.gaugedata.gauges.append([97, 95.1888, 29.8937, 0, 1.e10])
    rundata.gaugedata.gauges.append([98, 95.2179, 29.8796, 0, 1.e10])
    rundata.gaugedata.gauges.append([99, 95.2087, 29.8596, 0, 1.e10])
    rundata.gaugedata.gauges.append([100, 95.2887, 29.8152, 0, 1.e10])
    rundata.gaugedata.gauges.append([101, 95.3596, 29.7087, 0, 1.e10])
    rundata.gaugedata.gauges.append([102, 95.3754, 29.5804, 0, 1.e10])
    rundata.gaugedata.gauges.append([103, 95.0346, 30.1296, 0, 1.e10])
    rundata.gaugedata.gauges.append([104, 95.0512, 30.1254, 0, 1.e10])
    rundata.gaugedata.gauges.append([105, 95.1062, 29.9179, 0, 1.e10])
    rundata.gaugedata.gauges.append([106, 95.2829, 29.8312, 0, 1.e10])

    #Deposits
    rundata.gaugedata.gauges.append([1, 94.910801, 29.048672, 0, 1.e10])
    rundata.gaugedata.gauges.append([3, 94.90578199, 29.050469, 0, 1.e10])
    rundata.gaugedata.gauges.append([5, 94.90157822, 28.99481, 0, 1.e10])
    rundata.gaugedata.gauges.append([7, 94.89917192, 28.984449, 0, 1.e10])
    rundata.gaugedata.gauges.append([9, 94.86564025, 28.961909, 0, 1.e10])
    rundata.gaugedata.gauges.append([11, 94.77323329, 28.90818692, 0, 1.e10])
    rundata.gaugedata.gauges.append([13, 95.06934523, 28.576593, 0, 1.e10]) 
    rundata.gaugedata.gauges.append([15, 95.07955999, 28.40177, 0, 1.e10])
    rundata.gaugedata.gauges.append([17, 95.05756851, 28.354163, 0, 1.e10])
    rundata.gaugedata.gauges.append([19, 94.99651629, 28.23502186, 0, 1.e10])
    rundata.gaugedata.gauges.append([21, 95.03097289, 28.173463, 0, 1.e10])
    rundata.gaugedata.gauges.append([23, 95.14357644, 28.13399179, 0, 1.e10])


    #End
    rundata.gaugedata.gauges.append([107, 95.293702, 28.1296, 0, 1.e10])

    #Lake
    rundata.gaugedata.gauges.append([108, 94.912645, 30.1947, 0, 1.e10])


    return rundata
    # end of function setrun
    # ----------------------


#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    try:
        geo_data = rundata.geo_data
    except:
        print "*** Error, this rundata has no geo_data attribute"
        raise AttributeError("Missing geo_data attribute")

##### coordinate_system = 1 for Cartesian x-y in meters, coordinate_system = 2 for latitude-longitude on the sphere

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.04
    geo_data.friction_depth = 1000.0

    # Refinement data
    refinement_data = rundata.refinement_data
    refinement_data.wave_tolerance = 1.e-2
    refinement_data.deep_depth = 1e2
    refinement_data.max_level_deep = 4
    refinement_data.variable_dt_refinement_ratios = True

    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]
    
    #
    topo_data.topofiles.append([3, 1, 5, 0., 1.e10, 'srtm010616fill.txt'])

    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [topotype, minlevel,maxlevel,fname]

##### 0 = No perturbation specified, 1 = Perturbation to depth h, 4 = Perturbation to surface level

    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]

    # == setfixedgrids.data values ==
    fixedgrids = rundata.fixed_grid_data.fixedgrids
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]

    # == fgmax.data values ==
    fgmax_files = rundata.fgmax_data.fgmax_files
    # for fixed grids append to this list names of any fgmax input files
    rundata.fgmax_data.num_fgmax_val = 5
    #fgmax_files.append('fgmax2000_dep.txt')
    #fgmax_files.append('fgmax2000_inun.txt')
    #fgmax_files.append('fgmax2000_landslide.txt')

    return rundata
    # end of function setgeo
    # ----------------------



if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()

