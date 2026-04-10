{ pkgs ? import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/refs/tags/24.11.tar.gz") { config.allowUnfree = true; } }:

let
  # Libraries needed for both build and runtime
  allLibs = [
    pkgs.hdf5
    pkgs.hdf5-fortran
    pkgs.blas
    pkgs.lapack
  ];
in
pkgs.mkShell {
  nativeBuildInputs = [
    pkgs.pkg-config 
    pkgs.gfortran
    pkgs.fortran-fpm
  ];

  buildInputs = allLibs;

  shellHook = ''
    alias fpm="fortran-fpm"
    
    # Build-time and Run-time: Help gfortran and the dynamic linker find all libraries
    export LIBRARY_PATH="${pkgs.lib.makeLibraryPath allLibs}:$LIBRARY_PATH"
    export LD_LIBRARY_PATH="${pkgs.lib.makeLibraryPath allLibs}:$LD_LIBRARY_PATH"
    
    # Optional: Automatically enter zsh if not already in it
    if [[ $- == *i* && $SHELL != *"zsh"* ]]; then
      exec zsh
    fi
  '';
}
