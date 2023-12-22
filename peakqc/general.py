import beartype
from beartype.typing import Optional

import os
import matplotlib.pyplot as plt

@beartype
def _save_figure(prefix: str,
                 path: Optional[str],
                 dpi: int = 600) -> None:
    """Save the current figure to a file.

    Parameters
    ----------
    prefix : str
        Prefix to be used for the file name.
    path : Optional[str]
        Path to the file to be saved.
    dpi : int, default 600
        Dots per inch. Higher value increases resolution.
    """

    # 'path' can be None if _save_figure was used within a plotting function, and the internal 'save' was "None".
    # This moves the checking to the _save_figure function rather than each plotting function.
    if path is None:
        output_path = os.path.join(os.getcwd(), prefix + ".png")
        plt.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    else:
        output_path = os.path.join(path, prefix + ".png")
        plt.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")