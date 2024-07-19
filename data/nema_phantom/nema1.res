


              SIMIND Monte Carlo Simulation Program    V8.0  
------------------------------------------------------------------------------
 Phantom S : h2o       Crystal...: nai       InputFile.: nema              
 Phantom B : perspex   BackScatt.: lucite    OutputFile: nema1             
 Collimator: pb_sb2    SourceRout: smap      SourceImg.: nema              
 Cover.....: al        ScoreRout.: none      DensityImg: nema              
------------------------------------------------------------------------------
 PhotonEnergy.......: 140          Matrix    PhotonsPerProj....: 12579273       
 EnergyResolution...: 9            Spectra   Activity..........: 1              
 MaxScatterOrder....: 3            gi-lehr   DetectorLenght....: 25             
 DetectorWidth......: 0            SPECT     DetectorHeight....: 0.9525         
 UpperEneWindowTresh: 150.5        Random    Distance to det...: 19.478         
 LowerEneWindowTresh: 129.5        Phantom   ShiftSource X.....: 0              
 PixelSize  I.......: 0.4664       Resolut   ShiftSource Y.....: 0              
 PixelSize  J.......: 0.4664       Forced    ShiftSource Z.....: 0              
 HalfLength S.......: 11           Header    HalfLength P......: 11             
 HalfWidth  S.......: 0            SaveMap   HalfWidth  P......: 0              
 HalfHeight S.......: 0                      HalfHeight P......: 0              
 SourceType.........: 1ByteCodeMap           PhantomType.......: 1ByteCodeMap 
------------------------------------------------------------------------------
 GENERAL DATA
 keV/channel........: 0.5                    CutoffEnergy......: 0              
 Photons/Bq.........: 0.891                  StartingAngle.....: 0              
 CameraOffset X.....: 0                      CoverThickness....: 0              
 CameraOffset Y.....: 0                      BackscatterThickn.: 0              
 MatrixSize I.......: 128                    IntrinsicResolut..: 0.36           
 MatrixSize J.......: 128                    AcceptanceAngle...: 2.4432         
 Emission type......: 2                      Initial Weight....: 0.07083        
 NN ScalingFactor...: 1                      Energy Channels...: 512            
                                                                              
 SPECT DATA
 RotationMode.......: -360                   Nr of Projections.: 64             
 RotationAngle......: 5.625                  Projection.[start]: 1              
 Orbital fraction...: 0                      Projection...[end]: 64             
 Center of Rotation File: nema1.cor
                                                                              
 COLLIMATOR DATA FOR ROUTINE: Analytical          
 CollimatorCode.....: gi-lehr                CollimatorType....: Parallel 
 HoleSize X.........: 0.1212                 Distance X........: 0.012          
 HoleSize Y.........: 0.13995                Distance Y........: 0.08037        
 CenterShift X......: 0.0666                 X-Ray flag........: F              
 CenterShift Y......: 0.11535                CollimThickness...: 3.28           
 HoleShape..........: Hexagonal              Space Coll2Det....: 0              
 CollDepValue [57]..: 0                      CollDepValue [58].: 0              
 CollDepValue [59]..: 0                      CollDepValue [60].: 0              
                                                                              
 IMAGE-BASED PHANTOM DATA
 RotationCentre.....: 183,183                Bone definition...: 1100           
 CT-Pixel size......: 0.1                    Slice thickness...: 0.2            
 StartImage.........: 1                      No of CT-Images...: 110            
 MatrixSize I.......: 364                    CTmapOrientation..: 0              
 MatrixSize J.......: 364                    StepSize..........: 0.1            
 CenterPoint I......: 183                    ShiftPhantom X....: 0              
 CenterPoint J......: 183                    ShiftPhantom Y....: 0              
 CenterPoint K......: 56                     ShiftPhantom Z....: 0              
                                                                              
 INFO FOR TCT file
 MatrixSize I.......: 128                    MatrixSize J......: 128            
 MatrixSize K.......: 128                    Units.............: g/cm3*1000          
 Scout File.........: F
                                                                              
 PHANTOM DATA FROM FILE: phantom.zub section: 3

 Code Volume          Density   Voxels Volume(mL)        MBq     MBq/mL   Value
   1: sphere 1          1.000      844  0.169E+01  0.403E-02  0.238E+01  60.000
   2: sphere 2          1.000     1774  0.355E+01  0.705E-02  0.199E+01  50.000
   3: sphere 3          1.000     3662  0.732E+01  0.131E-01  0.179E+01  45.000
   4: sphere 4          1.000     7747  0.155E+02  0.246E-01  0.159E+01  40.000
   5: sphere 5          1.000    13397  0.268E+02  0.373E-01  0.139E+01  35.000
   6: sphere 6          1.000    51992  0.104E+03  0.124E+00  0.119E+01  30.000
   7: background        1.000  4545244  0.909E+04  0.723E+00  0.795E-01   2.000
   8: spine/lung        0.300   169224  0.338E+03  0.673E-01  0.199E+00   5.000
                                                                              
 Phantom volume....L: 9.59    
 Phantom mass.....Kg: 9.35    
                                                                              
 INTERACTIONS IN THE CRYSTAL
 MaxValue spectrum..: 40.15          
 MaxValue projection: 0.3590E-01     
 CountRate spectrum.: 41.77          
 CountRate E-Window.: 18.87          
                                                                              
 SCATTER IN ENERGY WINDOW
 Scatter/Primary....: 0.32081        
 Scatter/Total......: 0.24289        
 Scatter order 1....: 87.22 %        
 Scatter order 2....: 11.64 %        
 Scatter order 3....: 1.13 %         
                                                                              
 CALCULATED DETECTOR PARAMETERS
 Efficiency E-window: 0.4158         
 Efficiency spectrum: 0.9207         
 Sensitivity Cps/MBq: 18.8669        
 Sensitivity Cpm/uCi: 41.8846        
                                                                              
 Simulation started.: 2024:07:10 14:29:59
 Simulation stopped.: 2024:07:10 15:49:06
 Elapsed time.......: 1 h, 19 m and 7 s
 DetectorHits.......: 31521520       
 DetectorHits/CPUsec: 6647           
                                                                              
 OTHER INFORMATION
 Compiled 2024:05:03 with INTEL Mac   
 Random number generator: Intel RAN
 Comment:EMISSION
 Energy resolution as function of 1/sqrt(E)
 Header file: nema1.h00
 Linear angle sampling within acceptance angle
 Inifile: simind.ini
 Command: nema nema1/fz:phantom/45:3/tr:11/tr:15
