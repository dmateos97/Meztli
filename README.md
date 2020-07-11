# Meztli - Minimum length bell nozzle design
Method of characteristics nozzle optimization code:

This code uses the Prandtl-Meyer expansion function to find the maximum
angle by which a given exhaust flow can turn at the throat for a minimum
length nozzle design. The method of characteristics is then used to
parametrize the expansion fan and find the axisymmetric flow deflections
and new slopes in order to find the optimized nozzle wall contour. This
method delivers axial exit flow with uniform Mach number for given Pc/Pe
and specific heat ratio inputs. This code consists of 5 sections which
perform the following functionalities:
 
Section 1 - Input section:
    Creates GUI for modular user input parameters
 
Section 2 - Preliminary parameter calculations:
    Uses nozzle thermodynamic relations to find the exit Mach number for
    given combustion chamber conditions, and finds the maximum vacuum
    thrust coefficient Cfv for an ideal and reversible vacuum expansion
    process. Given user input thrust a vacuum optimized throat radius is
    found (Used if rt = 0 in section 1)
   
Section 3 - Expansion fan:
    The Prandtl-Meyer function is used to find and plot the expansion fan
    
Section 4 - Flow deflections:
    The MOC is used to find and plot the intersections of the left
    running and right running characteristics and the flow deflections at
    each intersection of the crossing expansion fans.
 
Section 5 - Wall intersections:
    The wall intersections are found given the final slopes computed in
    section 4 & the wall slopes starting with theta max found as a
    function of the Prandtl-Meyer function. The wall slopes iterate n
    times until they reach 0 at the nozzle exit to ensure an evenly
    distributed flow. The results are plotted to provide the nozzle
    contour.
