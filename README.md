# ABAQUS-US
## Author
- [Juan David Gómez Cataño](http://www.eafit.edu.co/docentes-investigadores/Paginas/juan-gomez.aspx), Professor at Universidad EAFIT.
This REPO contains a variety of ABAQUS user element (UEL) and user material (UMAT) subroutines. A list of input files and the related subroutine is defined in the file versheet.dat.In the future I will add more details on how to use the subroutines. The subroutines inlcude classical and Cosserat elemnts with various constitutive models in the form of UMATs.

## Installation
Download the subroutines and run the Abaqus with the input file (`.inp`) that accompanies the subroutines.

    abaqus job=filename user=user_routine

### Example
This is a uniaxial test with classical plasticity

    abaqus job=UNIUSER_CLA_KIN user=UMAT_PCLK

## License
It is licensed under the [MIT license](http://en.wikipedia.org/wiki/MIT_License).
