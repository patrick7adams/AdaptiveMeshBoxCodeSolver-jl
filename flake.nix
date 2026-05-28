# Pin nixpkgs for reproducibility. You can swap to a channel like nixos-24.05.
{
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.05";
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };

        # Common tools you want in every class shell
        commonTools = with pkgs; [
          git
          gnumake
          gnuplot
          pkg-config
          # add more global tools here (ripgrep, fd, just, etc.)
        ];

        # A Python you can share across classes; tune packages as you like
        commonPython =
          pkgs.python312.withPackages (ps: with ps; [
            pip
            numpy
            scipy
            matplotlib
            jupyterlab
            sympy
            # add pandas, networkx, etc. if you want them everywhere
          ]);

        # Helper to build a class shell with extra packages layered on top
        mkClassShell = {extraPkgs ? [], extraHook ? ""}: 
          pkgs.mkShell {
          packages = commonTools ++ [ commonPython ] ++ extraPkgs;

          # Put small niceties here that you want in *every* shell
          shellHook = ''
            export PYTHONNOUSERSITE=1
            export PIP_DISABLE_PIP_VERSION_CHECK=1
            echo "Entered $(pwd) dev shell... goon"
            ${extraHook}
          '';
          };
      in {
        devShells = {
          # fallback if you run `nix develop` with no attribute
          default = mkClassShell {};

          # === Per-class shells (edit these to your needs) ===

          P4est = mkClassShell {
            extraPkgs = [ 
              pkgs.gcc 
              pkgs.cmake 
              pkgs.gdb 
              pkgs.libtool 
              pkgs.automake 
              pkgs.autoconf
              pkgs.mpi
              pkgs.tmpi
              pkgs.paraview
              pkgs.gmsh
            ];
            extraHook = "echo 'welcome to research city!67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n67\n'";
          };

        };
      });
}

