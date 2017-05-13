import numpy as np
import re
import sys

log_file = sys.argv[1]

with open(log_file, "r") as f:
    lines = f.readlines()

    nums = [[], [], [], []]
    i = 0
    while i < len(lines):
        if lines[i][:6] == "FGMRES":
            n_iter = float(re.findall("\d+", lines[i])[-3])
            print(n_iter)

            for j in range(1, 5):
                line = lines[i+j]
                nums[j-1].append((float(re.findall("\d+\.\d+", line)[0]))
                                 / n_iter)

            i += 5
        else:
            i += 1

    # print("Krylov =", nums[0])
    print("Krylov mean =", np.mean(nums[0]), "std =", np.std(nums[0]))
    print("LLS mean =", np.mean(nums[1]), "std =", np.std(nums[1]))
    print("Calculate Residual mean =", np.mean(nums[2]), "std =", np.std(nums[2]))
    print("Total mean =", np.mean(nums[3]), "std =", np.std(nums[3]))
