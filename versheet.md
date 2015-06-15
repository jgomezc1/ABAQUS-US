| Problem            | File               | Observations                                     | User Subroutine                |
|:-------------------|:------------------:|:-------------------------------------------------|:------------------------------:|
| Uniaxial extension | `UNIUSER_CLA`      | Uniaxial extension                               | `UEL8_PCLI_R.for`              |
|                    |                    | Classical plasticity theory                      |                                |
|                    |                    |                                                  |                                |
| Uniaxial extension | `UNIUSER_COS`      | Uniaxial extension                               | `UEL8_PCOI.for UEL8_PCOR.for`  |
|                    |                    | Cosserat plasticity theory                       |                                |
|                    |                    | Radial return and return mapping algorithms      |                                |
|                    |                    |                                                  |                                |
| Uniaxial extension | `UNIUSER9CS`       | Uniaxial extension                               | `UEL9_PCOI.for UEL9_PCOR.for`  |
|                    |                    | Cosserat plasticity theory                       |                                |
|                    |                    | Radial return and return mapping algorithms      |                                |
|                    |                    |                                                  |                                |
| Uniaxial extension | `UNIUSER_CLA_ELA`  | Uniaxial extension                               | `UEL8_ECL.for`                 |
|                    |                    | Classical elasticity                             |                                |
|                    |                    |                                                  |                                |
| Uniaxial extension | `UNIUSER_CLA_CAE`  | Uniaxial extension                               | `Abaqus model`                 |
|                    |                    | Classical elasticity                             |                                |
|                    |                    |                                                  |                                |
| Uniaxial extension | `UNIUSER_COS_ELA`  | Uniaxial extension                               | `UEL8_COS.for`                 |
|                    |                    | Cosserat elasticity                              |                                |
|                    |                    |                                                  |                                |
| Beam bending       | `STOLKEN`          | Beam with solid elements                         | `UEL9_PCOI.for UEL9_PCOR.for`  |
|                    |                    |                                                  |                                |
| Uniaxial extension | `UNIUSER_CLA_UEL9` | Uniaxial extension                               | `UEL9_PCLK.for`                |
|                    |                    | Classical plasticity theory                      |                                |
|                    |                    | Combined hardening                               |                                |
|                    |                    |                                                  |                                |
| Uniaxial extension | `UNIUSER_CLA_UEL8` | Uniaxial extension                               | `UEL8_PCLK.for`                |
|                    |                    | Classical plasticity theory                      |                                |
|                    |                    | Combined hardening                               |                                |
|                    |                    |                                                  |                                |
| Uniaxial extension | `UNIUSER_CLA_KIN`  | Uniaxial extension                               | `UMAT_PCLK.for`                |
|                    |                    | Classical plasticity theory                      |                                |
|                    |                    | Combined hardening                               |                                |
|                    |                    |                                                  |                                |
| Uniaxial extension | `UNIUSER_UEL8`     | Uniaxial extension                               | `UEL8_PCLK.for`                |
|                    |                    | Classical plasticity theory                      |                                |
|                    |                    | Combined hardening different boundary conditions |                                |
|                    |                    |                                                  |                                |
| Uniaxial extension | `UNIUSER_COS_KIN`  | Uniaxial extension                               | `UEL8_PCLK_KIN.for`            |
|                    |                    | Cosserat plasticity theory                       |                                |
|                    |                    | Combined hardening                               |                                |
|                    |                    |                                                  |                                |
| Beam bending       | `STOLKEN_CLA`      | Beam bending                                     | `UEL9_PCLK.for`                |
|                    |                    | Classical plasticity theory                      |                                |
|                    |                    | Combined hardening                               |                                |
|                    |                    |                                                  |                                |
| Beam bending       | `STOLKEN_KIN`      | Uniaxial extension                               | `UEL9_PCOR_KIN.for`            |
|                    |                    | Cosserat plasticity theory                       |                                |
|                    |                    | Combined hardening                               |                                |
|                    |                    |                                                  |                                |
| Uniaxial extension | `UNIAXIAL_AXI_UEL` | Uniaxial extension                               | `UEL8_ECL_AXY.for`             |
|                    |                    | Classical elasticity                             |                                |
|                    |                    | Axisymmetric formulation                         |                                |