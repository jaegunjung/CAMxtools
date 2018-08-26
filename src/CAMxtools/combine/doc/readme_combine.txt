
Version update history

  - After formatted for CAMxtools package
    - Formatted for github release and drop version name from the file. Updates will be recorded by date.
    - 12/27/2017, jjung - Moved writing to file functionality to wrt_ioapi.py and removed the blank_ioapi function. Closed all input files.
    - 12/22/2017, jjung - Cleaned unnecessary lists to understand the script easily. Added more consistent comments. blank_ioapi is moved to other file.
  - After being shared with modelers by placing it on /models folder
    - v1.1 introduced main function to call the body of program as a function, "combine" so that other scripts can use "combine" function and formatted for github style with changes of variable name for easy understanding.
    - v1 cleaned the update logs in the header of the script. Species definition file can have blank line with white columns. "!" is okay to use for the commenting out lines (10/05/2017,jjung).
  - Used for projects but before being located on /models direcotry
    - a15 can handle uam input file (a14 generated an issue when uam files are fed by introducing xarray). Also, when it handles multiple input files, start was not an integer. Thus, I forced the variable type to integer. This is a new issue, which was not seen before. I guess adding xarray library may cause some issue... I don't know. (10/3/2017, jjung).
    - a14 improved the processing speed by introducing list comprehension and xarray. Also switched from python2 to python3. Do not create file too large (>= 2GB). Use original FORTRAN based combine in the case. For xarray, http://xarray.pydata.org/en/stable/.                                                                  
    - a13 the customized nc file does not have TSTEP. Fixed this issue.
    - a12 will not handle variables starting with number.
    - a11 fix a new issue from a8-a10 when underscore(_) is used in species name. Now it is resolved. But it will be a problem if species name include other species name such as HONO_DD has NO_DD. This causes a problem when adding "v" in front of species name from species definition file.                                     
    - a10 version fix a new bug from a6-a9. NLAYS was no. of vars instead of no. of layers. This bug is fixed.
    - a9 version fix a new bug in a8. The bug was repeating v as many as the same species name appear over the multiple input files.
    - a8 version can handle when variable names start with number such as DDM output variable. Additional regular expression is added, and internal variable names begin with "v" (6/23/2017,jjung).                                                  
    - a7 version updates main loop and fixes a bug when var[0] was used in spec_def. The bug was using always 0 to 24 of outputs instead of using the current time stamp.                                                                             
    - a6 version allows to specify only surface layer and fixed potential bugs if input files are 3D.
    - a5 version fixs TFLAG when the dates are spread from the last day of the previous year to the first day of the new year.
    - a4 version handles customized nc file as an input file (3/19/2017, jjung).
    - a3 version fixs when the date inlcudes more than two years (2/25/2017, jjung).
    - a2 version allow overriding values from input rather than starting from the end of existing output file (2/23/2017, jjung).
    - a1 version fixes file attributes so that m3tools can be applied (2/22/2017, jjung).
    - a0 version allows appending of output files (2/17/2017, jjung).
    - created the first version (1/27/2017, jjung).
