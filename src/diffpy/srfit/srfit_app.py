import argparse

from diffpy.srfit.version import __version__


def main():
    parser = argparse.ArgumentParser(
        prog="diffpy.srfit",
        description=(
            "Generalized code base for modeling problems.\n\n"
            "For more information, visit: "
            "https://github.com/diffpy/diffpy.srfit/"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--version",
        action="store_true",
        help="Show the program's version number and exit",
    )

    args = parser.parse_args()

    if args.version:
        print(f"diffpy.srfit {__version__}")
    else:
        # Default behavior when no arguments are given
        parser.print_help()


if __name__ == "__main__":
    main()
