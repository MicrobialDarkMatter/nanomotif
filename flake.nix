{
  description = "A flake for nanomotif dev in nix";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = [
            pkgs.cargo
            pkgs.rustc
            pkgs.linuxKernel.packages.linux_zen.perf
            pkgs.cargo-cross
            pkgs.cargo-release

            # For building
            pkgs.openssl
            pkgs.clang
            pkgs.pkg-config
            pkgs.libGL
            pkgs.stdenv.cc.cc.lib
            pkgs.glib

            # For python package
            pkgs.python3
            pkgs.python3Packages.python-lsp-server
            pkgs.python3Packages.flake8
            pkgs.python3Packages.pylint
            pkgs.python3Packages.pyflakes

            pkgs.pixi
          ];
          shellHook = ''
            export LD_LIBRARY_PATH="${pkgs.libGL}/lib/:${pkgs.stdenv.cc.cc.lib}/lib/:${pkgs.glib.out}/lib/:$LD_LIBRARY_PATH"
          '';
        };
        packages.default = pkgs.buildRustPackage {
          src = ./.;
          cargoLock = {
            lockFile = ./Cargo.lock;
          };
          buildInputs = [
            pkgs.openssl
            pkgs.pkg-config
          ];
          nativeBuildInputs = [
            pkgs.openssl_3_3
            pkgs.pkg-config
          ];
        };
      }
    );
}
