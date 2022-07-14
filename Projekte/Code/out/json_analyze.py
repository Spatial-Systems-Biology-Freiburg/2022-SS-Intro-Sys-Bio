from pathlib import Path
import numpy as np
import json


if __name__ == "__main__":
    
    count_stability_true = 0
    solving_reasons = {}
    count_stability_none = 0
    count_solving_true = 0
    count_total = 0

    for i, filename in enumerate(sorted(Path().rglob("*.json"))):
        print("[{:010.0f} Analyzing ...".format(i), end="\r")
        file = open(filename, "r")
        madata = json.load(file)

        count_total += 1
        msg = madata["solving"]["message"] 
        if msg in solving_reasons.keys():
            solving_reasons[msg] += 1
        else:
            solving_reasons[msg] = 0
        # stability_reasons[]
        if madata["LSA"]["succ"] == True:
            count_stability_true += 1
        elif madata["LSA"]["succ"] == None:
            count_stability_none += 1
        if madata["solving"]["success"] == True:
            count_solving_true += 1
        file.close()
    print()
    print("Total count:  {}".format(count_total))
    print("LSA Success:  {}".format(count_stability_true))
    print("LSA None:     {}".format(count_stability_none))
    print("Solving Succ: {}".format(count_solving_true))
    print(solving_reasons)
