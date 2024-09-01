from shared.hpc.commands import map_1, map_2, thermo, dock, reduce, g16
import dill
import sys
from rdkit import RDLogger
import os

RDLogger.DisableLog("rdApp.*")

if len(sys.argv) < 3:
    print("Usage: python apply.py <data_file_index> <dir>")
    print("intended for use in a distributed computing environment")
    sys.exit(1)

command = sys.argv[1]
index_value = sys.argv[2]
dir = sys.argv[3]  # /nobackup1/bmahjour/distributed
print(command, dir, index_value)


if command == "all":
    for i in ["map_1", "map_2", "thermo", "dock", "reduce"]:
        if not os.path.exists(f"{dir}/{i}_output_data/{index_value}"):
            os.mkdir(f"{dir}/{i}_output_data/{index_value}")

        input_file_loc = f"{dir}/{i}_input_data/input_{index_value}.txt"

        if not os.path.exists(input_file_loc):
            print(f"File {input_file_loc} does not exist")
            sys.exit(0)

        with open(input_file_loc, "r") as f:
            data = [x.strip() for x in f.readlines()]

        if i == "map_1":
            map_1(dir, data, index_value)
        if i == "map_2":
            map_2(dir, data, index_value)
        if i == "thermo":
            thermo(dir, data, index_value)
        if i == "dock":
            dock(dir, data, index_value)

elif command == "g16":
    input_file_loc = f"{dir}/input_data/input_{index_value}.txt"
    if not os.path.exists(input_file_loc):
        print(f"File {input_file_loc} does not exist")
        sys.exit(0)
    with open(input_file_loc, "r") as f:
        data = [x.strip() for x in f.readlines()]

    g16(dir, data, index_value)


else:

    input_file_loc = f"{dir}/{command}_input_data/input_{index_value}.txt"

    if not os.path.exists(input_file_loc):
        print(f"File {input_file_loc} does not exist")
        sys.exit(0)

    with open(input_file_loc, "r") as f:
        data = [x.strip() for x in f.readlines()]

    if not os.path.exists(f"{dir}/{command}_output_data/{index_value}"):
        os.mkdir(f"{dir}/{command}_output_data/{index_value}")

    if command == "map_1":
        map_1(dir, data, index_value)

    if command == "map_2":
        map_2(dir, data, index_value)

    if command == "reduce":
        reduce(dir, data, index_value)

print("all done!")
