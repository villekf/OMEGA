#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This example file can be used for seamless GATE --> OMEGA reconstruction
The GATE simulation is performed first where the output data is stored as
ROOT files. The offline coincidence sorter is then used to create coincidences.
A sinogram is then created by OMEGA from the simulated coincidences and then 
reconstructed.
This example uses the Philips Vereos PET scanner, but any cylindrical PET
scanner can be used here.
"""
import opengate as gate
import opengate.contrib.pet.philipsvereos as pet_vereos
from opengate.userhooks import check_production_cuts


if __name__ == "__main__":
    
    """
    This block contains only GATE 10 related code.
    
    GATE CODE BEGINS
    """
    
    sim = gate.Simulation()
    
    sim.g4_verbose = False
    sim.visu = False
    sim.check_volumes_overlap = False
    sim.random_seed = 123456789
    sim.output_dir = ''
    sim.progress_bar = True
    sim.number_of_threads = 16

    # units
    m = gate.g4_units.m
    mm = gate.g4_units.mm
    cm = gate.g4_units.cm
    Bq = gate.g4_units.Bq
    MBq = Bq * 1e6
    sec = gate.g4_units.second
    keV = gate.g4_units.keV

    #  change world size
    world = sim.world
    world.size = [3 * m, 3 * m, 3 * m]
    world.material = "G4_AIR"

    # add a PET VEREOS
    sim.volume_manager.add_material_database("GateMaterials_pet.db")
    pet = pet_vereos.add_pet(sim, "pet", create_housing=False, create_mat=False)
            
    crystal = sim.volume_manager.volumes[f"{pet.name}_crystal"]
    
    # Get the number of blocks, stacks, die and crystals
    nBlocks = int(sim.volume_manager.volumes[f"{pet.name}_module"].translation.shape[0])
    nStacksTr = int(sim.volume_manager.volumes[f"{pet.name}_module"].size[1] // sim.volume_manager.volumes[f"{pet.name}_stack"].size[1])
    nStacksAx = int(sim.volume_manager.volumes[f"{pet.name}_module"].size[2] // sim.volume_manager.volumes[f"{pet.name}_stack"].size[2])
    nDieTr = int(sim.volume_manager.volumes[f"{pet.name}_stack"].size[1] // sim.volume_manager.volumes[f"{pet.name}_die"].size[1])
    nDieAx = int(sim.volume_manager.volumes[f"{pet.name}_stack"].size[2] // sim.volume_manager.volumes[f"{pet.name}_die"].size[2])
    nCrystalTr = int(sim.volume_manager.volumes[f"{pet.name}_die"].size[1] // sim.volume_manager.volumes[f"{pet.name}_crystal"].size[1])
    nCrystalAx = int(sim.volume_manager.volumes[f"{pet.name}_die"].size[2] // sim.volume_manager.volumes[f"{pet.name}_crystal"].size[2])
    
    # physics
    sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option4"
    sim.physics_manager.set_production_cut("world", "all", 1 * m)

    # cuts
    reg2 = sim.physics_manager.add_region("reg2")
    reg2.production_cuts.all = 0.1 * mm
    reg2.associate_volume(f"{pet.name}_crystal")

    # Add to rod sources, no phantom
    source = sim.add_source('GenericSource', 'mysource')
    source.particle = "back_to_back"
    source.energy.type == "mono"
    source.energy.mono = 511 * keV
    source.activity = 5 * MBq / sim.number_of_threads
    source.position.type = "cylinder"
    source.position.radius = 15 * mm
    source.position.dz = 164 * mm / 2.0
    source.position.translation = [-3 * cm, -7 * cm, 0 * cm]
    source.direction.type = 'iso'

    source2 = sim.add_source('GenericSource', 'mysource2')
    source2.particle = "back_to_back"
    source2.energy.type == "mono"
    source2.energy.mono = 511 * keV
    source2.activity = 10 * MBq / sim.number_of_threads
    source2.position.type = "cylinder"
    source2.position.radius = 50 * mm
    source2.position.dz = 120 * mm / 2.0
    source2.position.translation = [7 * cm, 0 * cm, 0 * cm]
    source2.direction.type = 'iso'

    # add stat actor
    s = sim.add_actor("SimulationStatisticsActor", "Stats")
    s.track_types_flag = True

    # set user hook function
    sim.user_hook_after_init = check_production_cuts
    
    # Necessary hits collection
    hc = sim.add_actor("DigitizerHitsCollectionActor", "Hits")
    hc.attached_to = crystal.name
    hc.authorize_repeated_volumes = True
    hc.attributes = ["EventID", "PostPosition", "TotalEnergyDeposit", "PreStepUniqueVolumeID", "GlobalTime"]
    hc.root_output.write_to_disk = False

    # singles collection
    sc = sim.add_actor("DigitizerAdderActor", "Singles")
    sc.input_digi_collection = "Hits"
    # sc.policy = "EnergyWinnerPosition"
    sc.policy = "EnergyWeightedCentroidPosition"
    sc.root_output.write_to_disk = False

    # Energy windowing
    ew = sim.add_actor("DigitizerEnergyWindowsActor", "EnergyWindows")
    ew.attached_to = hc.attached_to
    ew.input_digi_collection = "Singles"
    ew.authorize_repeated_volumes = True
    ew.channels = [{"name": ew.name, "min": 450 * keV, "max": 613 * keV}]
    ew.attributes = ["EventID", "TotalEnergyDeposit", "PreStepUniqueVolumeID", "GlobalTime"]
    ew.output_filename = sim.output_dir + 'PETVereosOMEGA.root'
    ew.root_output.write_to_disk = True
    
    sim.run_timing_intervals = [[0 * sec, 30 * sec]]
    
    sim.run()
    stats = sim.get_actor("Stats")
    print(stats)
    
    import os, uproot
    import numpy as np
    from opengate.actors.coincidences import coincidences_sorter
    # Remove old coincidences file if it exists
    coincidenceFile = sim.output_dir + "coincidences.root"
    if os.path.exists(coincidenceFile):
        os.remove(coincidenceFile)
    
    # Create coincidences using GATE's offline sorter
    root_file = uproot.open(sim.output_dir + "PETVereonOMEGA.root")
    singles_tree = root_file["EnergyWindows"]
    ns = gate.g4_units.nanosecond
    time_window = 2 * ns
    # policy = "takeAllGoods"
    policy = "takeWinnerOfGoods"
    # policy = "removeMultiples"
    mm = gate.g4_units.mm
    transaxial_plane = "xy"
    min_transaxial_distance = 128 * mm
    max_axial_distance = 164 * mm
    
    # Apply coincidence sorter
    coincidences = coincidences_sorter(
        singles_tree,
        time_window,
        policy,
        min_transaxial_distance,
        transaxial_plane,
        max_axial_distance,
        chunk_size=5000000,
        return_type='dict',
        output_file_path= sim.output_dir + 'coincidences.root',
        # output_file_path=None,
        output_file_format='root'
    )
    
    """
    GATE CODE ENDS
    """
    
    """
    This block contains mainly OMEGA related code
    
    OMEGA CODE BEGINS
    """
    # Open the coincidence-file
    root_file = uproot.open(sim.output_dir + "coincidences.root")
    coincidences = root_file["Coincidences"]
    
    # Load necessary submodules
    from omegatomo.projector import proj
    from omegatomo.reconstruction import reconstructions_main
    from omegatomo.util.checkCUDA import checkCUDA
    
    # Create the projectorClass class object
    options = proj.projectorClass()

    # Needed to compute the diameter of the scanner
    dist = sim.volume_manager.volumes[f"{pet.name}_module"].translation
    # Scanner parameters neede for reconstruction
    # modules in axial direction
    options.blocks_per_ring = nBlocks
    # Number of crystals in axial direction
    options.cryst_per_block_z = nCrystalAx * nDieAx
    # Number of stacks in axial direction
    options.linear_multip = nStacksAx
    # number of crystals in transaxial direction
    options.cryst_per_block = nCrystalTr * nStacksTr * nDieTr
    # Number of detectors per ring
    detPerRing = int(options.blocks_per_ring * options.cryst_per_block)
    # Span factor/axial compression
    options.span = 1
    # Number of radial positions (views) in sinogram
    options.Ndist = 520
    # Number of angles (tangential positions) in sinogram
    options.Nang = detPerRing // 2
    # Number of crystal rings
    options.rings = options.cryst_per_block_z * options.linear_multip
    # Maximum ring difference
    options.ring_difference = options.rings - 1
    
    # Initialize the sinogram
    options.SinM = np.zeros(0, dtype=np.uint16)
    from omegatomo.util.sinogram import saveSinogram
    
    # Load the coincidences and create the sinogram in chunks
    for chunk in coincidences.iterate(step_size=70000000, library="pd"):
    
        num_entries = chunk.shape[0]
        
        joku1 = chunk["PreStepUniqueVolumeID1"].astype(str)
        
        # Get the module, stack, die and crystal numbers
        pattern = r'-([0-9]{1,2})_([0-9]{1,2})_([0-9]{1,2})_([0-9]{1,2})_([0-9]{1,2})_([0-9]{1,2})$'
        
        # Extract using regex
        extracted = joku1.str.extract(pattern)
        
        # Convert extracted columns to NumPy arrays with dtype uint16
        XX1 = extracted[2].astype(np.uint16).to_numpy()
        YY1 = extracted[3].astype(np.uint16).to_numpy()
        ZZ1 = extracted[4].astype(np.uint16).to_numpy()
        TT1 = extracted[5].astype(np.uint16).to_numpy()
        
        del joku1
        
        joku2 = chunk["PreStepUniqueVolumeID2"].astype(str)
        
        extracted = joku2.str.extract(pattern)
        
        XX2 = extracted[2].astype(np.uint16).to_numpy()
        YY2 = extracted[3].astype(np.uint16).to_numpy()
        ZZ2 = extracted[4].astype(np.uint16).to_numpy()
        TT2 = extracted[5].astype(np.uint16).to_numpy()
        
        del extracted
        del joku2
            
        koko = XX1.size
        
        # Ring number of the first photon
        ring_number1 = np.uint16(np.mod(YY1, nStacksAx) * nCrystalAx * nDieAx + 
                              ((np.mod(ZZ1, nDieAx)) * nCrystalAx + (np.mod(TT1, nCrystalAx))))
        
        # ring number of the second photon
        ring_number2 = np.uint16(np.mod(YY2, nStacksAx) * nCrystalAx * nDieAx + 
                              ((np.mod(ZZ2, nDieAx)) * nCrystalAx + (np.mod(TT2, nCrystalAx))))
        
        # Transaxial ring position of the first photon
        ring_pos1 = np.uint16(np.mod(XX1, nBlocks) * nCrystalTr * nStacksTr * nDieTr + 
                          np.floor(TT1 / nCrystalAx) + np.floor(YY1 / nStacksAx) * nCrystalTr * nDieTr 
                          + np.floor(ZZ1 / nDieAx) * nCrystalTr)
        
        # Transaxial ring position of the second photon
        ring_pos2 = np.uint16(np.mod(XX2, nBlocks) * nCrystalTr * nStacksTr * nDieTr + 
                          np.floor(TT2 / nCrystalAx) + np.floor(YY2 / nStacksAx) * nCrystalTr * nDieTr 
                          + np.floor(ZZ2 / nDieAx) * nCrystalTr)
        
        # Create the sinogram
        # The data is appended to the previous sinogram, if non-empty
        options.SinM = saveSinogram(ring_pos1, ring_pos2, ring_number1, ring_number2, options.Nang, options.Ndist, options.ring_difference, 
                           span = options.span, rings = options.rings, Sino = options.SinM)
    
    
    ### Crystal pitch/size in x- and y-directions (transaxial) (mm)
    options.cr_p = sim.volume_manager.volumes[f"{pet.name}_crystal"].size[1]
    
    ### Crystal pitch/size in z-direction (axial) (mm)
    options.cr_pz = sim.volume_manager.volumes[f"{pet.name}_crystal"].size[2]
    
    ### Ring diameter (distance between perpendicular detectors) (mm)
    options.diameter = (np.abs(dist[np.argmin(np.abs(dist[:,0])),1]) - sim.volume_manager.volumes[f"{pet.name}_module"].size[0] / 2.) * 2
    
    ### The size of the gap(s) between stacks (mm)
    gap = np.abs(sim.volume_manager.volumes[f"{pet.name}_stack"].translation[0][2]-sim.volume_manager.volumes[f"{pet.name}_stack"].translation[1][2])
    options.ringGaps = np.tile(gap - options.cryst_per_block_z * options.cr_pz, options.linear_multip - 1)
    
    ### Transaxial FOV size (mm), this is the length of the x (row) side
    # of the FOV (same is used for the column direction)
    options.FOVa_x = 576
    
    ### Axial FOV (mm)
    options.axial_fov = 164
    
    ### Scanner name
    # Used for naming purposes (measurement data)
    options.machine_name = 'Philips_Vereos_PET_example'
    
    ### Reconstructed image pixel count (X/row-direction)
    options.Nx = 128
    
    ### Y/column-direction
    options.Ny = 128
    
    ### Z-direction (number of slices) (axial)
    options.Nz = 79
    
    ### Flip the image (in column direction)?
    options.flip_image = False
    
    ### How much is the image rotated?
    # You need to run the precompute phase again if you modify this
    # NOTE: The rotation is done in the detector space (before reconstruction).
    # This current setting is for systems whose detector blocks start from the
    # right hand side when viewing the device from front.
    # Positive values perform the rotation in clockwise direction
    options.offangle = 0
    
    ### Show status messages
    # These are e.g. time elapsed on various functions and what steps have been
    # completed. It is recommended to keep this True.
    options.verbose = 2
    
    ### Device used 
    # Uncomment the below lines and run them to determine the available device
    # numbers:
    # from omegatomo.util.devinfo import deviceInfo
    # deviceInfo(True)
    options.deviceNum = 0
    
    ### Use CUDA
    # Selecting this to True will use CUDA kernels/code instead of OpenCL. This
    # only works if the CUDA code was successfully built. Recommended only for
    # Siddon as the orthogonal/volume-based ray tracer are slower in CUDA.
    options.useCUDA = checkCUDA(options.deviceNum)
    
    ### Use CPU
    # Selecting this to True will use CPU-based code instead of OpenCL or CUDA.
    options.useCPU = False
     
    ############################### PROJECTOR #################################
    ### Type of projector to use for the geometric matrix
    # 0 = Regular Siddon's algorithm (only available with implementation 1 and
    # when precomputed_lor = False) NOT RECOMMENDED.
    # 1 = Improved/accelerated Siddon's algorithm
    # 2 = Orthogonal distance based ray tracer
    # 3 = Volume of intersection based ray tracer
    # See the wiki for more information:
    # https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
    options.projector_type = 1
    
    ### Use point spread function (PSF) blurring
    # Applies PSF blurring through convolution to the image space. This is the
    # same as multiplying the geometric matrix with an image blurring matrix.
    options.use_psf = True
    
    # FWHM of the Gaussian used in PSF blurring in all three dimensions
    # options.FWHM = [options.cr_p options.cr_p options.cr_pz]
    options.FWHM = np.array([options.cr_p, options.cr_p, options.cr_pz])
     
    ######################### RECONSTRUCTION SETTINGS #########################
    ### Number of iterations (all reconstruction methods)
    options.Niter = 4
    
    ### Number of subsets
    options.subsets = 8
    
    ### Subset type (n = subsets)
    # 1 = Every nth (column) measurement is taken
    # 2 = Every nth (row) measurement is taken (e.g. if subsets = 3, then
    # first subset has measurements 1, 4, 7, etc., second 2, 5, 8, etc.) 
    # 3 = Measurements are selected randomly
    # 4 = (Sinogram only) Take every nth column in the sinogram
    # 5 = (Sinogram only) Take every nth row in the sinogram
    # 8 = Use every nth sinogram
    # 9 = Randomly select the full sinograms
    # 11 = Use prime factor sampling to select the full sinograms
    #Most of the time subsetType 1 or 4 is sufficient.
    options.subsetType = 1
    
    ### Initial value for the reconstruction
    options.x0 = np.ones((options.Nx, options.Ny, options.Nz), dtype=np.float32)
     
    ############################### ML-METHODS ################################
    ### Ordered Subsets Expectation Maximization (OSEM) OR Maximum-Likelihood
    ### Expectation Maximization (MLEM) (if subsets = 1)
    options.OSEM = False
    
    ### Accelerated COSEM (ACOSEM)
    options.ACOSEM = False
     
     
    ############################### MAP-METHODS ###############################
    # Any algorithm selected here will utilize any of the priors selected below
    # this. Note that only one algorithm and prior combination is allowed! You
    # can also use most of these algorithms without priors (such as PKMA or
    # PDHG).
    ### Modified BSREM (MBSREM)
    options.MBSREM = False
    
    ### Preconditioned Krasnoselskii-Mann algorithm (PKMA)
    options.PKMA = True
    
    ### Primal-dual hybrid gradient (PDHG)
    options.PDHG = False
    
    ### Primal-dual hybrid gradient (PDHG) with Kullback-Leibler minimization
    options.PDHGKL = False
     
     
    ################################# PRIORS ##################################
    ### Median Root Prior (MRP)
    options.MRP = False
    
    ### Non-local Means (NLM) prior
    options.NLM = False
    
    ### Relative difference prior
    options.RDP = False
    
    
    ############################ ENFORCE POSITIVITY ###########################
    ### Applies to PDHG, PDHGL1, PDDY, FISTA, FISTAL1, MBSREM, MRAMLA, PKMA
    # Enforces positivity in the estimate after each iteration
    options.enforcePositivity = True
     
     
    ############################ ACOSEM PROPERTIES ############################
    ### Acceleration parameter for ACOSEM (1 equals COSEM)
    options.h = 2
    
    
    ########################## RELAXATION PARAMETER ###########################
    ### Relaxation parameter for MRAMLA, RAMLA, ROSEM, BSREM, MBSREM and PKMA
    # Use scalar if you want it to decrease as
    # lambda / ((current_iteration - 1)/20 + 1). Use vector (length = Niter) if
    # you want your own relaxation parameters. Use empty array or zero if you
    # want to OMEGA to compute the relaxation parameter using the above formula
    # with lamda = 1. Note that current_iteration is one-based, i.e. it starts
    # at 1.
    options.lambdaN = np.zeros(0, dtype=np.float32)
     
    
    ######################## MRAMLA & MBSREM PROPERTIES #######################
    ### Upper bound for MRAMLA/MBSREM (use 0 for default (computed) value)
    options.U = 0
     
    
    ############################# PKMA PROPERTIES #############################
    ### Step size (alpha) parameter for PKMA
    # If a scalar (or an empty) value is used, then the alpha parameter is
    # computed automatically as alpha_PKMA(oo) = 1 + (options.rho_PKMA *((i -
    # 1) * options.subsets + ll)) / ((i - 1) * options.subsets + ll +
    # options.delta_PKMA), where i is the iteration number and l the subset
    # number. The input number thus has no effect. options.rho_PKMA and
    # options.delta_PKMA are defined below.
    # If, on the other hand, a vector is input then the input alpha values are
    # used as is without any modifications (the length has to be at least the
    # number of iterations * number of subsets).
    options.alpha_PKMA = 0
    
    ### rho_PKMA
    # This value is ignored if a vector input is used with alpha_PKMA
    options.rho_PKMA = 0.95
    
    ### delta_PKMA
    # This value is ignored if a vector input is used with alpha_PKMA
    options.delta_PKMA = 1
    
    ############################# PDHG PROPERTIES #############################
    # Primal value
    # If left zero, or empty, it will be automatically computed
    options.tauCP = 0
    # Primal value for filtered iterations, applicable only if
    # options.precondTypeMeas[2] = True. As with above, automatically computed
    # if left zero or empty.
    options.tauCPFilt = 0
    # Dual value. Recommended to set at 1.
    options.sigmaCP = 1
    # Next estimate update variable
    options.thetaCP = 1
    
    # Use adaptive update of the primal and dual variables
    # Currently only one method available
    # Setting this to 1 uses an adaptive update for both the primal and dual
    # variables.
    # Can lead to unstable behavior with using multi-resolution
    # Minimal to none use with filtering-based preconditioner
    options.PDAdaptiveType = 0
    
    ############################# PRECONDITIONERS #############################
    ### Applies to PDHG, PDHGL1, PDHGKL, PKMA, MBSREM, MRAMLA, PDDY, FISTA and
    ### FISTAL1
    # Measurement-based preconditioners
    # precondTypeMeas(0) = Diagonal normalization preconditioner (1 / (A1))
    # precondTypeMeas(1) = Filtering-based preconditioner
    options.precondTypeMeas[0] = False
    if options.PDHG or options.PDHGKL:
        options.precondTypeMeas[1] = True
    
    # Image-based preconditioners
    # Setting options.precondTypeImage(1) = true when using PKMA, MRAMLA or
    # MBSREM is recommended
    # precondTypeImage(0) = Diagonal normalization preconditioner (division with
    # the sensitivity image 1 / (A^T1), A is the system matrix) 
    # precondTypeImage(1) = EM preconditioner (f / (A^T1), where f is the current
    # estimate) 
    # precondTypeImage(2) = IEM preconditioner (max(n, fhat, f)/ (A^T1), where
    # fhat is an estimate of the final image and n is a small positive number) 
    # precondTypeImage(3) = Momentum-like preconditioner (basically a step size
    # inclusion) 
    # precondTypeImage(4) = Gradient-based preconditioner (Uses the normalized
    # divergence (sum of the gradient) of the current estimate) 
    # precondTypeImage(5) = Filtering-based preconditioner
    # precondTypeImage(6) = Curvature-based preconditioner
    options.precondTypeImage[0] = False
    if options.PKMA or options.MBSREM:
        options.precondTypeImage[1] = True
    options.precondTypeImage[2] = False
    options.precondTypeImage[3] = False
    options.precondTypeImage[4] = False
    options.precondTypeImage[5] = False
    options.precondTypeImage[6] = False
    
    # Reference image for precondTypeImage(3). Can be either a mat-file or a
    # variable
    options.referenceImage = ''
    
    # Momentum parameter for precondTypeImage(4)
    # Set the desired momentum parameters to the following variable (note that
    # the length should be options.Niter * options.subsets): 
    # options.alphaPrecond = np.empty(0, dtype=np.float32)
    # Otherwise set the following parameters:
    options.rhoPrecond = options.rho_PKMA
    options.delta1Precond = options.delta_PKMA
    
    # Parameters for precondTypeImage(5)
    # See the article for details
    options.gradV1 = 1.5
    options.gradV2 = 2
    # Note that these include subiterations (options.Niter * options.subsets)
    options.gradInitIter = 8
    options.gradLastIter = 24
    
    # Number of filtering iterations
    # Applies to both precondTypeMeas(2) and precondTypeImage(6)
    options.filteringIterations = 16
    
    
    ######################### REGULARIZATION PARAMETER ########################
    ### The regularization parameter for ALL regularization methods (priors)
    options.beta = 1
     
     
    ######################### NEIGHBORHOOD PROPERTIES #########################
    ### How many neighboring pixels are considered 
    # With MRP, QP, L, FMH, NLM, GGMRF and weighted mean
    # E.g. if Ndx = 1, Ndy = 1, Ndz = 0, then you have 3x3 square area where
    # the pixels are taken into account (I.e. (Ndx*2+1)x(Ndy*2+1)x(Ndz*2+1)
    # area).
    # NOTE: Currently Ndx and Ndy must be identical.
    # For NLM this is often called the "search window".
    options.Ndx = 2
    options.Ndy = 2
    options.Ndz = 1
     
     
    ############################## NLM PROPERTIES #############################
    ### Filter parameter
    options.NLMsigma = 900
    
    ### Patch radius
    options.Nlx = 1
    options.Nly = 1
    options.Nlz = 1
    
    ### Standard deviation of the Gaussian filter
    options.NLM_gauss = 2
    
    # Search window radius is controlled by Ndx, Ndy and Ndz parameters
    # Use anatomical reference image for the patches
    options.NLM_use_anatomical = False
    
    ### Specify filename for the reference image here (same rules apply as with
    # attenuation correction above)
    options.NLM_reference_image = 'reference_image.mat'
    
    # Note that only one of the below options for NLM can be selected!
    ### Use Non-local total variation (NLTV)
    # If selected, will overwrite regular NLM regularization as well as the
    # below MRP version.
    options.NLTV = False
    
    ### Use MRP algorithm (without normalization)
    # I.e. gradient = im - NLM_filtered(im)
    options.NLM_MRP = False
    
    ### Use non-local relative difference prior (NLRD)
    options.NLRD = False
    
    ### Use non-local GGMRF (NLGGMRF)
    options.NLGGMRF = True
    
    
    ############################## RDP PROPERTIES #############################
    ### Edge weighting factor
    options.RDP_gamma = 1
     
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ###########################################################################
    
    ###########################################################################
    ########################## DEPTH OF INTERACTION ###########################
    ###########################################################################
    # Uncomment the below value to set a depth of interaction (mm)
    # NOTE: By default this is set to 0, i.e. it is assumed that all the
    # interactions occur at the surface of the detector crystals. What this
    # value changes is the depth of where the interactions are assumed to
    # occur, i.e. it only changes the detector coordinates such that the
    # transaxial coordinates are "deeper" in the crystal.
    # options.DOI = 0
    ###########################################################################
    
        
    
    # pz, the reconstructed image
    # fp, the forward projections (if selected, i.e. options.storeFP = True)
    pz, fp = reconstructions_main(options)
    
    """
    OMEGA CODE ENDS
    """
    
    # import matplotlib as plt
    
    # plt.pyplot.imshow(pz[:,:,39])
        
    from scipy.io import savemat

    mdic = {"pz": pz}

    savemat(sim.output_dir + "PETVereosRecon.mat", mdic)
        