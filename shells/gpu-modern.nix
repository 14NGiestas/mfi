{ pkgs ? import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/refs/tags/24.11.tar.gz") { config.allowUnfree = true; } }:

let
  cudaPkgs = pkgs.cudaPackages_12_8;
  cudaLibs = [
    cudaPkgs.libcublas
    cudaPkgs.libcublas.dev
    cudaPkgs.cuda_cudart
    cudaPkgs.cuda_cudart.dev
    cudaPkgs.cuda_nvcc
  ];
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
    export LIBRARY_PATH="${pkgs.lib.makeLibraryPath allLibs}:$LIBRARY_PATH"
    export LD_LIBRARY_PATH="${pkgs.lib.makeLibraryPath allLibs}:$LD_LIBRARY_PATH"
  '';
}
