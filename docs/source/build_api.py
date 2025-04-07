""" Script to generate the API for the peakqc package."""

import os
import glob


##############################################################################
# ---------------------------- Helper functions ---------------------------- #
##############################################################################

def hline():
    line = "-" * 30 + "\n\n"
    return line


def header(text, level):

    if level == 0:  # doc title
        s = "=" * len(text) + "\n"
        s += f"{text}\n"
        s += "=" * len(text) + "\n\n"

    else:
        level_symbol = {1: "=", 2: "-", 3: "~"}

        s = f"{text}\n"
        s += f"{level_symbol[level] * len(text)}\n\n"

    return s


def automodule(name):

    s = f"""
.. automodule:: {name}
    :members:
    :undoc-members:
    :show-inheritance:
    """
    s += "\n"

    return s


def get_modules(path):

    modules = [os.path.basename(f).replace(".py", "") for f in glob.glob(path + "/*")]
    modules = [f for f in modules if not f.startswith("_") and f != "data"]

    return modules


##############################################################################
# ---------------------------- Configure API ------------------------------- #
##############################################################################

package_dir = os.path.dirname(os.path.abspath(__file__)) + "/../../peakqc/"

##############################################################################
# ------------- Main function for generating API documentation ------------- #
##############################################################################

def main():

    # Generate the API documentation per module
    modules = get_modules(package_dir)

    # Create one page per module
    for module in modules:
        print(module)
        with open("API/" + module + ".rst", 'w') as fp:

            fp.write(header(module.capitalize(), 1))

            fp.write(automodule(f"peakqc.{module}"))

            fp.write(hline())


if __name__ == "__main__":
    main()
