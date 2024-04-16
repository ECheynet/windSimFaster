# windSimFaster: WindSimFast, but faster

WindSimFaster is a variant of the WindSimFast repository, focusing on simulating spatially correlated wind fields more efficiently. Unlike WindSimFast, which simultaneously simulates all three velocity components (u, v, and w), WindSimFaster assumes these components are independent. Therefore, one component can be simulated at a time.  Both WindSImFast and WindSimFaster should be able to generate non-Gaussian and non-homogeneous wind fields. Both tools rely on the assumption of stationary turbulence, with WindSimFaster introducing spatial correlation through the coherence of turbulence after simulating first uncorrelated wind histories.

## Summary

WindSimFaster generates spatially correlated wind fields, focusing on the efficiency of simulations. It provides also a framework for using hindcast data to derive a more realistic turbulent wind field (See Documentation2.mlx). 

## Content

This repository includes:

- **Functions Folder**: functions for wind field simulation:
  - `KaimalModel.m`: Generates spectral densities according to the Kaimal model.
  - `getSamplingPara.m`: Calculates time and frequency vectors.
  - `windSimFaster.m`: Core function for generating spatially correlated wind fields.
  - `coherence.m`: Estimates coherence for validation purposes.
  - `randomProcess.m`: Creates time-domain one-point random processes from given power spectral densities.

- **Tutorials**:
  - `Documentation1.mlx`: A basic usage guide for WindSimFaster
  - `Documentation2.mlx`: Utilizes NORA3 hindcast data to simulate more realistic turbulent wind fields at Sørlige Nordsjø II

In both tutorials, only the along-wind velocity component is simulated. A similar approach can be sued for the other two components.

Any comment, suggestion or question is welcomed. This is the first version of the repository. Som bugs may still be present.

## References

WindSImFast is available at https://github.com/ECheynet/windSimFast

