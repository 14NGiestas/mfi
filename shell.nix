{ pkgs ? import <nixpkgs> { config.allowUnfree = true; } }:

let
  # Specific CUDA package set for consistency
  cudaPkgs = pkgs.cudaPackages;
  
  # Libraries needed for both build and runtime
  cudaLibs = [
    cudaPkgs.libcublas
    cudaPkgs.libcusolver
    cudaPkgs.cuda_cudart
    cudaPkgs.cuda_nvcc
  ];
in
pkgs.mkShell {
  nativeBuildInputs = [
    pkgs.pkg-config 
    pkgs.fortran-fpm
  ];

  buildInputs = [
    pkgs.hdf5
    pkgs.hdf5-fortran
    pkgs.blas
    pkgs.lapack
  ] ++ cudaLibs;

  shellHook = ''
    # Set aliases manually since programs.zsh doesn't work here
    alias fpm="fortran-fpm"
    
    # Build-time: Help gfortran find -lcublas
    export LIBRARY_PATH="${pkgs.lib.makeLibraryPath cudaLibs}:$LIBRARY_PATH"
    
    # Run-time: Help the executable find libcublas.so
    export LD_LIBRARY_PATH="${pkgs.lib.makeLibraryPath cudaLibs}:$LD_LIBRARY_PATH"
    
    # Optional: Automatically enter zsh if not already in it
    if [[ $- == *i* && $SHELL != *"zsh"* ]]; then
      exec zsh
    fi
  '';
}
