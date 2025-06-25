import pathlib
import subprocess
from tempfile import TemporaryDirectory
from typing import Optional, Tuple

from predict.config import DEFAULT_IQTREE_EXE
from predict.custom_errors import IQTreeError

def run_iqtree_command(cmd: list[str]) -> None:
    """Helper method to run an IQ-TREE command.

    Args:
        cmd (list): List of strings representing the IQ-TREE command to run.

    Raises:
        IQTreeError: If the IQ-TREE command fails with a CalledProcessError.
        RuntimeError: If the IQ-TREE command fails with any other error.
    """
    try:
        subprocess.check_output(cmd, encoding="utf-8")
    except subprocess.CalledProcessError as e:
        raise IQTreeError(subprocess_exception=e)
    except Exception as e:
        raise RuntimeError("Running IQ-TREE command failed.") from e


# lấy số cuối của dòng
def _get_value_from_line(line: str, search_string: str) -> float:
    line = line.strip()
    if search_string in line:
        _, value = line.rsplit(" ", 1)
        return float(value)
    raise ValueError(
        f"The given line '{line}' does not contain the search string '{search_string}'."
    )

def get_iqtree_rfdist_results(rfdist_file: str) -> Tuple[int, float]:
    """
    Parse IQ-TREE RF distance matrix file to compute:
    - Number of unique topologies (số cây)
    - Average normalized RF distance

    Args:
        rfdist_file (str): Path to .rfdist file

    Returns:
        Tuple[int, float]: (number of trees, average RF distance)
    """
    with open(rfdist_file) as f:
        lines = f.readlines()

    n = int(lines[0].split()[0])
    if n < 2:
        raise ValueError("Không đủ cây để tính RF")
    
    # Đọc ma trận
    matrix = [
        list(map(float, line.strip().split()[1:]))  # Bỏ tên Tree0, Tree1...
        for line in lines[1:]
    ]

    # Tính trung bình RF trên tam giác trên
    total = 0.0
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            total += matrix[i][j]
            count += 1

    avg_rf = total / count if count > 0 else 0.0
    return n, avg_rf


class IQTree:
    """Class to interact with the IQ-TREE executable.

    Args:
        exe_path (Optional[pathlib.Path]): Path to the IQ-TREE executable.

    Raises:
        FileNotFoundError: If the executable does not exist.
        RuntimeError: If the executable is not a working IQ-TREE binary.
    """

    def __init__(self, exe_path: Optional[pathlib.Path] = DEFAULT_IQTREE_EXE):
        if exe_path is None or not exe_path.exists():
            raise FileNotFoundError("IQ-TREE executable not found.")

        try:
            out = subprocess.check_output([str(exe_path.absolute())], encoding="utf-8")
        except Exception as e:
            raise RuntimeError(
                f"Running `{exe_path}` failed: {e}"
            ) from e

        if "IQ-TREE" not in out:
            raise RuntimeError(
                f"The executable `{exe_path}` is not a valid IQ-TREE binary."
            )

        self.exe_path = exe_path

    def _base_cmd(self, msa_file: pathlib.Path, model: str, prefix: pathlib.Path, **kwargs) -> list[str]:
        additional_settings = []
        for key, value in kwargs.items():
            if value is None:
                additional_settings.append(f"--{key}")
            else:
                additional_settings.extend([f"--{key}", str(value)])

        return [
            str(self.exe_path.absolute()),
            "-s", str(msa_file.absolute()),
            "-m", model,
            "-pre", str(prefix.absolute()),
            *additional_settings,
        ]

    def _run_rfdist(self, trees_file: pathlib.Path, prefix: pathlib.Path, **kwargs) -> None:
        additional_settings = []
        for key, value in kwargs.items():
            if value is None:
                additional_settings.append(f"--{key}")
            else:
                additional_settings.extend([f"--{key}", str(value)])

        cmd = [
            str(self.exe_path.absolute()),
            "--rf_all",
            "-", str(trees_file.absolute()),
            "-pre", str(prefix.absolute()),
            *additional_settings,
        ]
        run_iqtree_command(cmd)

    def infer_parsimony_trees(
        self,
        msa_file: pathlib.Path,
        model: str,
        prefix: pathlib.Path,
        n_trees: int = 24,
        **kwargs,
    ) -> pathlib.Path:
        """
        Infer a set of trees using maximum parsimony.

        Args:
            msa_file (pathlib.Path): Path to the MSA file.
            model (str): Substitution model (e.g., "GTR+G", "LG+G").
            prefix (pathlib.Path): Output prefix for the run.
            n_trees (int): Number of trees to infer.
            **kwargs: Additional flags for IQ-TREE.

        Returns:
            pathlib.Path: Path to the resulting startTree file.
        """
        cmd = self._base_cmd(
            msa_file, model, prefix, start=None, tree=f"pars{{{n_trees}}}", **kwargs
        )
        run_iqtree_command(cmd)
        return pathlib.Path(f"{prefix}.iqtree.startTree")

    def get_rfdistance_results(
        self,
        trees_file: pathlib.Path,
        prefix: Optional[pathlib.Path] = None,
        **kwargs,
    ) -> tuple[int, float]:
        """
        Compute the number of unique topologies and relative RF-Distance.

        Args:
            trees_file (pathlib.Path): Path to the file with trees.
            prefix (Optional[pathlib.Path]): Prefix for IQ-TREE output.
            **kwargs: Additional flags for IQ-TREE.

        Returns:
            Tuple[int, float]: (num_topos, rel_rfdist)
        """
        with TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            if not prefix:
                prefix = tmpdir / "rfdist"
            self._run_rfdist(trees_file, prefix, **kwargs)
            log_file = pathlib.Path(f"{prefix}.iqtree.log")
            return get_iqtree_rfdist_results(log_file)
