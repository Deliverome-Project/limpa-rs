# SPDX-License-Identifier: GPL-3.0-or-later
import argparse
import json

from .runtime import runtime_info, setup_r, validate


def main():
    parser = argparse.ArgumentParser(
        description="Install and check the pinned LIMPA runtime"
    )
    parser.add_argument("command", choices=["setup-r", "doctor", "validate"])
    args = parser.parse_args()
    try:
        if args.command == "setup-r":
            print(setup_r())
        elif args.command == "doctor":
            print(json.dumps(runtime_info(), indent=2))
        else:
            validate()
    except (RuntimeError, OSError) as exc:
        parser.exit(1, f"{exc}\n")


if __name__ == "__main__":
    main()
