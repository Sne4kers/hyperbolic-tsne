import os

import numpy as np


def find_last_embedding(log_path):
    for subdir, dirs, files in reversed(list(os.walk(log_path))):
        for fi, file in enumerate(reversed(sorted(files, key=lambda x: int(x.split(", ")[0])))):
            root, ext = os.path.splitext(file)
            if ext == ".csv":
                total_file = subdir.replace("\\", "/") + "/" + file

                return np.genfromtxt(total_file, delimiter=',')