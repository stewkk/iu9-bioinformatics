{
  inputs.nixpkgs = {
    url = "github:NixOS/nixpkgs";
  };
  inputs.flake-utils.url = "github:numtide/flake-utils";
  inputs.flake-compat = {
    url = "github:edolstra/flake-compat";
    flake = false;
  };

  outputs = { self, nixpkgs, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };
        pythonEnv = pkgs.python313.withPackages (ps: [
          ps.pip
          ps.virtualenv
          ps.pytest
        ]);
      in {
        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            pythonEnv
            pyright
          ];

          shellHook = ''
    export VENV_DIR="$PWD/.venv"
    if [ ! -d "$VENV_DIR" ]; then
      ${pythonEnv}/bin/python -m venv $VENV_DIR
      source $VENV_DIR/bin/activate
      pip install pip setuptools wheel
    else
      source $VENV_DIR/bin/activate
    fi

    export PYTHONPATH="${pythonEnv}/lib/python3.13/site-packages:$PYTHONPATH"

    if [ -f requirements.txt ]; then
      pip install -r requirements.txt
    fi

    python --version
  '';
        };
      });
}

