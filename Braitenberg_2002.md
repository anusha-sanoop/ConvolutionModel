# Inverse Modelling of Elastic Thickness by Convolution Method – The Eastern Alps as a Case Example

**Authors:**  
Carla Braitenberg¹, Jörg Ebbing², Hans-Jürgen Götze²

**Affiliations:**  
¹ Department of Earth Sciences, University of Trieste, Via Weiss 1, 34100 Trieste, Italy  
² Institut für Geologische Wissenschaften, Geophysik, Freie Universität Berlin, Malteserstrasse 74-100, D-12249 Berlin, Germany

**Publication:**  
Earth and Planetary Science Letters, 202(3-4), 387-404 (2002)

---

## Summary

This study presents a method for determining spatial variations of the effective elastic thickness (Tₑ) of the lithosphere in the Eastern Alps using a convolution-based inverse modeling approach. The research integrates gravity data, topographic information, and detailed crustal density models to estimate the mechanical properties of the lithosphere.

**Key Findings:**
- Tₑ increases towards the Molasse Basin in the north and the Po Basin in the south
- The main Alpine range exhibits lower Tₑ values
- These variations reflect differences in lithospheric strength and rigidity, correlating with geological and tectonic features

---

## Methodology

The convolution method employed in this study involves the following steps:

### 1. Data Collection and Processing
- **Gravity Data:** Collection of terrain-corrected Bouguer gravity data over the study area
- **Data Processing:** 
  - Stripping the Bouguer field of contributions from the upper 10 km of the crust
  - Application of a low-pass filter with a cut-off period of 75 km
  - Extraction of crust-mantle interface (CMI) undulations

### 2. Forward Modeling
- Creation of a detailed 3D crustal density model
- Simulation of CMI undulations based on the density structure

### 3. Convolution Method
The core of the methodology involves:
- **Preparation of Unit Impulse Flexure Response Functions:** These functions represent the lithospheric response to a unit load and are essential for modeling the flexural behavior of the lithosphere
- **Convolution with Topographic Loads:** The unit impulse response functions are convolved with the actual topographic loads to compute the expected flexural response of the lithosphere
- **Comparison with Observed Data:** The modeled flexural response is compared with observed gravity data and crustal structure information obtained from independent geophysical investigations
- **Inverse Modeling:** Adjustment of flexure parameters, including Tₑ, to achieve the best fit between the modeled and observed data

### 4. Inverse Modeling Approach
The inverse modeling process:
- Fits the crustal structure (derived from independent geophysical investigations) with the thin plate flexure model
- Estimates spatial variations in Tₑ by analyzing the lithosphere's response to loading
- Utilizes Green's functions to relate surface loads to lithospheric deflections

---

## Key Equations

### 1. Convolution Integral for Lithospheric Deflection

The lithospheric deflection is modeled using a convolution integral:

\[ w(x, y) = \iint R(x - x', y - y') \cdot q(x', y') \, dx' \, dy' \]

Where:
- \( w(x, y) \) = lithospheric deflection at coordinates \( (x, y) \)
- \( R(x - x', y - y') \) = flexural response function (unit impulse response)
- \( q(x', y') \) = applied load at coordinates \( (x', y') \)
- The integral represents the convolution operation

### 2. Flexural Rigidity

The flexural rigidity \( D \) is related to the effective elastic thickness \( T_e \) by:

\[ D = \frac{E T_e^3}{12(1 - \nu^2)} \]

Where:
- \( D \) = flexural rigidity
- \( E \) = Young's modulus
- \( T_e \) = effective elastic thickness
- \( \nu \) = Poisson's ratio

### 3. Simplified Flexural Response

In simplified form, the flexural response can be expressed as:

\[ w(x) = \frac{q(x)}{D} * R(x) \]

Where:
- \( w(x) \) = deflection of the lithosphere
- \( q(x) \) = applied load
- \( D \) = flexural rigidity
- \( R(x) \) = unit impulse response function
- \( * \) = convolution operation

### 4. Relationship Between Load and Deflection

The relationship between the applied load and the resulting deflection:

\[ w(x, y) = \frac{q(x, y)}{D} \]

This equation illustrates that the lithospheric deflection is directly proportional to the applied load and inversely proportional to the flexural rigidity.

---

## Technical Details

### Data Processing Parameters
- **Crustal stripping depth:** Upper 10 km of the crust
- **Low-pass filter cut-off:** 75 km period
- **Target interface:** Crust-mantle interface (CMI)

### Modeling Approach
- The method uses Green's functions to relate surface loads to lithospheric deflections
- The convolution approach allows for efficient computation of flexural responses
- The inverse problem is solved by iteratively adjusting Tₑ to minimize the misfit between observed and modeled data

---

## Results and Interpretation

### Spatial Variations of Tₑ
- **Northern region (Molasse Basin):** Higher Tₑ values, indicating greater lithospheric strength
- **Southern region (Po Basin):** Higher Tₑ values, suggesting increased rigidity
- **Central Alpine range:** Lower Tₑ values, reflecting reduced lithospheric strength

### Geological Implications
- Variations in Tₑ correlate with tectonic features and crustal heterogeneities
- Lower Tₑ in the main Alpine range may reflect:
  - Active tectonic processes
  - Thermal weakening
  - Crustal thickening
- Higher Tₑ in the foreland basins suggests:
  - More rigid lithosphere
  - Less active deformation
  - Different thermal structure

---

## Conclusions

The convolution method proves effective in estimating spatial variations of Tₑ in complex geological settings like the Eastern Alps. The approach successfully:

1. Integrates multiple geophysical datasets (gravity, topography, crustal structure)
2. Provides spatially continuous estimates of Tₑ
3. Reveals correlations between lithospheric strength and tectonic features
4. Offers insights into the mechanical behavior of the lithosphere

The findings contribute to understanding:
- The mechanical properties of the lithosphere
- The response of the lithosphere to tectonic forces
- The relationship between lithospheric structure and surface geology

---

## References

Braitenberg, C., Ebbing, J., & Götze, H.-J. (2002). Inverse modelling of elastic thickness by convolution method – the eastern Alps as a case example. *Earth and Planetary Science Letters*, 202(3-4), 387-404.

---

## Notes

- The convolution method allows for efficient forward modeling of lithospheric flexure
- The inverse approach enables estimation of Tₑ without requiring a priori assumptions about its spatial distribution
- The method is particularly useful in regions with complex tectonic histories and varying lithospheric properties
- Integration of gravity data with crustal structure models provides strong constraints on the flexural response

