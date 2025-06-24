import subprocess


class PyPythiaException(Exception):
    pass


class RAxMLNGError(Exception):
    """
    Custom RAxML-NG Exception used when running RAxML-NG commands.
    In case of a subprocess.CalledProcessError, the output of this Exception is either the entire RAxML-NG output,
    or only the lines containing the cause for the RAxML-NG error if the RAxML-NG output contains
    lines starting with "ERROR"
    """

    def __init__(self, subprocess_exception: subprocess.CalledProcessError):
        error_information = []

        for line in subprocess_exception.output.split("\n"):
            if line.strip().startswith("ERROR"):
                error_information.append(line.strip())

        if len(error_information) > 0:
            error_details = "RAxML-NG exited with the following error:\n"
            error_details += "\t" + "\n\t".join(error_information)
        else:
            error_details = "check the following RAxML-NG log for further information on the error(s):\n"
            error_details += subprocess_exception.output

        cmd = " ".join(subprocess_exception.cmd)
        self.message = f"Running RAxML-NG command failed: {cmd}\n" + error_details
        super().__init__(self.message)

class IQTreeError(Exception):
    """
    Custom IQ-TREE Exception used when running IQ-TREE commands.
    In case of a subprocess.CalledProcessError, the output of this Exception is either the entire IQ-TREE output,
    or only the lines containing the cause for the IQ-TREE error if the output contains lines starting with "ERROR".
    """

    def __init__(self, subprocess_exception: subprocess.CalledProcessError):
        error_information = []

        # Tách và kiểm tra từng dòng để tìm dòng báo lỗi
        for line in subprocess_exception.output.split("\n"):
            if line.strip().startswith("ERROR"):
                error_information.append(line.strip())

        if len(error_information) > 0:
            error_details = "IQ-TREE exited with the following error(s):\n"
            error_details += "\t" + "\n\t".join(error_information)
        else:
            error_details = "Check the following IQ-TREE log for further information on the error(s):\n"
            error_details += subprocess_exception.output

        cmd = " ".join(subprocess_exception.cmd)
        self.message = f"Running IQ-TREE command failed: {cmd}\n" + error_details
        super().__init__(self.message)