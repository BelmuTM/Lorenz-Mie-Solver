# Lorenz-Mie Solver

<a href="https://patreon.com/Belmu"><img src="https://i.imgur.com/XNsj1Qx.png" title="Patreon Page" width="15%" height="15%"></a>
<a href="https://discord.gg/jjRrhpkH9e"><img src="https://i.imgur.com/egvRtMf.png" title="Discord Server" width="15%" height="15%"></a>

<p align="left">
  <a>
  <img
    src="https://i.imgur.com/glQ08GO.png"
    width="40%"
    height="40%"
  ></a>
  <a>
  <img
    src="https://i.imgur.com/0NStZzn.png"
    width="55%"
    height="55%"
  ></a>
</p>

Implementation of the generalized Lorenz-Mie solution to Maxwell's equations. It outputs tabulated data for extinction/scattering cross sections, asymmetry factors
and phase functions (unpolarized, perpendicular, parallel) from a set of parameters:
1. Indices of refraction of particle(s) and host medium
2. Volume fraction of particle(s) in the host medium
3. Range of particle radii (e.g.: 1'000nm - 20'000nm)
4. Range of wavelengths (e.g.: 390nm - 830nm)
5. Integration resolutions for particle radii, wavelengths and angles

# References
   - [Bohren, C. F., & Huffman, D. R. (1983). Absorption and Scattering of Light by Small Particles.](https://onlinelibrary.wiley.com/doi/book/10.1002/9783527618156)
   - [Frisvad et al. (2007). Computing the Scattering Properties of Participating Media Using Lorenz-Mie Theory.](https://cseweb.ucsd.edu//~henrik/papers/lorenz_mie_theory/)
   - [Wikipedia. (2026). Legendre polynomials.](https://en.wikipedia.org/wiki/Legendre_polynomials)
   - [Wikipedia. (2026). Mie scattering.](https://en.wikipedia.org/wiki/Mie_scattering)

# Compiling

> [!IMPORTANT]
> C++ 20 (or later) and a 64-bit operating system are required.

### Windows System Requirements

- [CMake 3.22 (or later)](https://github.com/Kitware/CMake/releases/)
    * [In the installer] Check the box indicating: "Add CMake to the system `PATH`".

### Linux Package Dependencies

Ubuntu / Debian Package Dependencies

```
sudo apt-get update
sudo apt-get install -y cmake git
```

Arch Linux dependencies

```bash
sudo pacman -Sy cmake git
```

# License

> [!NOTE]  
> Lorenz-Mie Solver is licensed under the [GNU General Public License V3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).
> Consider reading the terms before modifying or redistributing this project.

<div align="center">
  Made with ❤️
</div>
