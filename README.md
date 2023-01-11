The code manages batch calculations with Gaussian quantum chemistry package.
This includes:
  
  - **Starting** calculations using
  
    -provided geometries and sample input files
    
    OR
    
    -given tree of existing calculation folders
    
    OR
    
    -pkl-file with calculation database
  
  - **Fixing** existing calculations using predetermined correction strategies applied to nested structure of input (the parameters usually have sub-parameters, e.g. _scf=(yqc, maxcycle=100))_, expressed via nested instructions on how to
  
    - Add/remove (sub)parameter
    
    - Increase/decrease numerical values of parameters
    
  - **Postprocessing** of results
    
