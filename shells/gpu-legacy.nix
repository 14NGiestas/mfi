{ pkgs ? import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/refs/tags/24.11.tar.gz") { config.allowUnfree = true; } }:

let
  # Use CUDA 11.8 which is compatible with driver 470.x (CUDA 11.4)
  cudaPkgs = pkgs.cudaPackages_11_8;
  
  # Libraries needed for both build and runtime
  cudaLibs = [
    cudaPkgs.libcublas
    cudaPkgs.libcublas.dev
    cudaPkgs.libcusolver
    cudaPkgs.libcusolver.dev
    cudaPkgs.cuda_cudart
    cudaPkgs.cuda_cudart.dev
    cudaPkgs.cuda_nvcc
  ];

  # All libraries to be included in build and runtime paths
  allLibs = [
    pkgs.hdf5
    pkgs.hdf5-fortran
    pkgs.blas
    pkgs.lapack
  ] ++ cudaLibs;
in
pkgs.mkShell {
  nativeBuildInputs = [
    pkgs.pkg-config 
    pkgs.gfortran
  ];

  buildInputs = allLibs;

  shellHook = ''
    # Build-time and Run-time: Help gfortran and the dynamic linker find all libraries
    export LIBRARY_PATH="${pkgs.lib.makeLibraryPath allLibs}:$LIBRARY_PATH"
    export LD_LIBRARY_PATH="${pkgs.lib.makeLibraryPath allLibs}:$LD_LIBRARY_PATH"
    
    # Optional: Automatically enter zsh if not already in it
    if [[ $- == *i* && $SHELL != *"zsh"* ]]; then
      exec zsh
    fi
  '';
}
