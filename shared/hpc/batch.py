import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import os
import itertools
import dill
import numpy as np


class DataLoader:
    def __init__(
        self,
        d=300,
        filter_seen=False,
        output_path="../distributed",
        data_path="../data",
    ):
        self.comm_data = []
        self.output_path = output_path
        self.data_path = data_path
        self.divisions = d
        self.filter_seen = filter_seen
        self.seen = []
        self.invald = []

    def load_data(self, file="sigma_compiled_sanitized_20230108.xlsx"):
        self.comm_data = pd.read_excel(f"{self.data_path}/{file}")

    def filter_smiles(self, smi):

        if smi in self.invald:
            return False

        if self.filter_seen:
            if smi in self.seen:
                return False

        if "." in smi:
            return False

        mol = Chem.MolFromSmiles(smi)

        if mol is None:
            return False

        mw = Chem.Descriptors.ExactMolWt(mol)
        if mw > 355:
            return False

        o = mol.GetAtomsMatchingQuery(Chem.AtomFromSmarts("[H]"))
        if len(o) > 0:
            return False

        return True

    # def _sample_1(self):
    #     n = self.comm_data.sample(1)
    #     sm = list(n.to_dict()["sanitized_smiles"].values())[0]
    #     o = self.filter_smiles(sm)
    #     while not o:
    #         if sm not in self.invald:
    #             self.invald.append(sm)
    #         n = self.comm_data.sample(1)
    #         sm = list(n.to_dict()["sanitized_smiles"].values())[0]
    #         o = self.filter_smiles(sm)
    #     self.seen.append(sm)
    #     return sm

    # def _sample_3(self):
    #     hits = []
    #     for i in range(3):
    #         n = self._sample_1()
    #         hits.append(n)
    #     return ".".join(hits)

    # def create_dataset(self):
    #     batch_size = self.dataset_size / self.divisions
    #     d_idx = 0
    #     n_idx = 0
    #     smiles_sets_collected = []
    #     for i in range(self.dataset_size):
    #         smi = self._sample_3()
    #         smiles_sets_collected.append(smi)
    #         if n_idx == batch_size:
    #             with open(
    #                 f"{self.output_path}/map_1_input_data/input_{d_idx}.txt", "w"
    #             ) as f:
    #                 f.write("\n".join(smiles_sets_collected))
    #             print(f"data shard {d_idx} of {self.divisions} created")
    #             smiles_sets_collected = []
    #             d_idx += 1
    #             n_idx = 0
    #         else:
    #             n_idx += 1

    def create_all_trios(self, percent=0.3):
        unique_smiles = []
        for idx, i in self.comm_data.iterrows():
            sm = Chem.CanonSmiles(i["smiles"])
            o = self.filter_smiles(sm)
            if not o:
                print(sm, "filtered")
                continue
            if sm not in unique_smiles:
                unique_smiles.append(sm)
        # print(len(unique_smiles))
        all_3mers = list(itertools.combinations(unique_smiles, 3))
        all_3mers = all_3mers[: int(len(all_3mers) * percent)]
        # print(len(all_3mers))
        self.dataset_size = len(all_3mers)
        batch_size = int(self.dataset_size / self.divisions)
        # print(batch_size)
        d_idx = 0
        n_idx = 0
        smiles_sets_collected = []
        for i in range(self.dataset_size):
            smi = ".".join(all_3mers[i])
            smiles_sets_collected.append(smi)
            if n_idx + 1 == batch_size:
                with open(
                    f"{self.output_path}/map_1_input_data/input_{d_idx}.txt", "w"
                ) as f:
                    f.write("\n".join(smiles_sets_collected))

                with open(
                    f"{self.output_path}/all_input_data/input_{d_idx}.txt", "w"
                ) as f:
                    f.write("\n".join(smiles_sets_collected))

                # print(f"data shard {d_idx} of {self.divisions} created")
                smiles_sets_collected = []
                d_idx += 1
                n_idx = 0
            else:
                n_idx += 1

        return self.dataset_size


def generate_sh(filename, index, dir, name, command):
    content = f"""#!/bin/bash
#SBATCH -n 8           # 1 core
#SBATCH -t 2-10:00:00   # d-hh:mm:ss
#SBATCH -J {command}_{index} # job name
#SBATCH --output=log_out/{name}/{command}_sample_{index}.log   # Standard output and error log. 
#SBATCH -p sched_mit_ccoley,sched_mit_ccoley_r8 # partition
#SBATCH --mem=20G # 10 gb

source /etc/profile.d/modules.sh

export OMP_NUM_THREADS=8
module load xtb
module load openbabel
module load orca
python -u ./shared/hpc/apply.py {command} {index} {dir} > log_out/{name}/{command}_sample_{index}.out

"""
    with open(filename, "w") as f:
        f.write(content)


class DistributedPipeline:
    def __init__(
        self,
        output_path="/nobackup1/bmahjour/distributed",
        data_path="../data",
        d=300,
        filter_seen=False,
    ):
        self.path = output_path
        self.data_path = data_path
        self.dl = DataLoader(
            d=d,
            filter_seen=filter_seen,
            output_path=self.path,
            data_path=self.data_path,
        )
        self.name = output_path.split("/")[-1]

    def _make_data_folders(self, stage="map_1"):
        if not os.path.exists(f"{self.path}/{stage}_input_data"):
            os.makedirs(f"{self.path}/{stage}_input_data")
        if not os.path.exists(f"{self.path}/{stage}_output_data"):
            os.makedirs(f"{self.path}/{stage}_output_data")
        if not os.path.exists(f"{self.path}/{stage}_slurm_scripts"):
            os.makedirs(f"{self.path}/{stage}_slurm_scripts")

    def _reset_data(self, stage="map_1"):
        if os.path.exists(f"{self.path}/{stage}_input_data"):
            os.system(f"rm -rf {self.path}/{stage}_input_data")
        if os.path.exists(f"{self.path}/{stage}_output_data"):
            os.system(f"rm -rf {self.path}/{stage}_output_data")
        if os.path.exists(f"{self.path}/{stage}_slurm_scripts"):
            os.system(f"rm -rf {self.path}/{stage}_slurm_scripts")

        self._make_data_folders(stage)

    def _create_sh_files(self, stage):
        d = self.dl.divisions

        if stage == "calculate_novelty_askcos":
            generate_sh(
                f"{self.path}/{stage}_slurm_scripts/output_script-{d}.sh",
                d,
                self.path,
                self.name,
            )
            return

        for i in range(0, d):
            generate_sh(
                f"{self.path}/{stage}_slurm_scripts/output_script-{i}.sh",
                i,
                self.path,
                self.name,
                stage,
            )

    def _submit(self, stage):
        if stage == "calculate_novelty_askcos":
            os.system(
                f"sbatch {self.path}/{stage}_slurm_scripts/output_script-{self.dl.divisions}.sh"
            )
            return

        for i in range(0, self.dl.divisions):
            os.system(f"sbatch {self.path}/{stage}_slurm_scripts/output_script-{i}.sh")
            # break

    # def create_dataset(self, file):
    #     self.dl.load_data(file)
    #     self.dl.create_dataset()
    #     print("made dataset")

    def load_data(self, file):
        self.dl.load_data(file)

    def create_dataset_all(self, percent=0.3):
        total_data_size = self.dl.create_all_trios(percent)
        print(
            self.dl.dataset_size,
            "experiments distributed over",
            self.dl.divisions,
            "divisions",
        )
        print("batch size:", int(self.dl.dataset_size / self.dl.divisions))

    def pipe(self, stage, clear_data=False, create_inputs=False, percent=0.3):
        if not os.path.exists(f"example/log/{self.name}"):
            os.makedirs(f"example/log/{self.name}")
        else:
            os.system(f"rm -rf example/log/{self.name}")
            os.makedirs(f"example/log/{self.name}")

        if stage == "all":
            for i in ["all", "map_1", "map_2", "thermo", "dock", "reduce"]:
                if clear_data:
                    self._reset_data(i)
                else:
                    self._make_data_folders(i)
        else:
            if clear_data:
                self._reset_data(stage)
            else:
                self._make_data_folders(stage)

        if create_inputs:
            self.create_dataset_all(percent)

        if stage == "reduce":
            with open(
                f"{self.data_path}/hits_from_psql_with_docking_20240807.pkl", "rb"
            ) as f:
                hits_w_docking_score = dill.load(f)

                n = len(hits_w_docking_score)
                n_per_d = n // self.dl.divisions

                print(n_per_d)

                for i in range(self.dl.divisions):
                    with open(
                        f"{self.path}/reduce_input_data/input_{str(i)}.txt", "w"
                    ) as f:
                        for j in range(i * n_per_d, (i + 1) * n_per_d):
                            f.write(str(hits_w_docking_score[j]) + "\n")

        self._create_sh_files(stage=stage)
        self._submit(stage)


class SingleNetworkDistributedPipeline:
    def __init__(
        self,
        path="/nobackup1/bmahjour/single_network_distributed",
        node_indices=[],
        d=300,
    ):
        self.path = path
        self.name = path.split("/")[-1]
        self.node_indices = node_indices
        self.d = d
        self.batches = []

        n_per_d = int(np.ceil(len(node_indices) / d))
        for i in range(d):
            l = node_indices[i * n_per_d : (i + 1) * n_per_d]
            if len(l) == 0:
                continue
            self.batches.append(l)

        print([len(b) for b in self.batches])
        print(len(self.batches), "batches")

    def _create_sh_files(self, stage):
        if not os.path.exists(f"{self.path}/{stage}_slurm_scripts"):
            os.makedirs(f"{self.path}/{stage}_slurm_scripts")

        for i in range(0, len(self.batches)):
            generate_sh(
                f"{self.path}/{stage}_slurm_scripts/output_script-{i}.sh",
                i,
                self.path,
                self.name,
                stage,
            )

    def _submit(self, stage):
        for i in range(0, len(self.batches)):
            os.system(f"sbatch {self.path}/{stage}_slurm_scripts/output_script-{i}.sh")

    def pipe(self, stage):
        if not os.path.exists(f"{self.path}"):
            os.makedirs(f"{self.path}")

        if not os.path.exists(f"{self.path}/input_data"):
            os.makedirs(f"{self.path}/input_data")
        else:
            os.system(f"rm -rf {self.path}/input_data")
            os.makedirs(f"{self.path}/input_data")

        for d_idx, batch in enumerate(self.batches):
            with open(f"{self.path}/input_data/input_{d_idx}.txt", "w") as f:
                for b in batch:
                    f.write(str(b) + "\n")

        if not os.path.exists(f"log_out/{self.name}"):
            os.makedirs(f"log_out/{self.name}")
        else:
            os.system(f"rm -rf log_out/{self.name}")
            os.makedirs(f"log_out/{self.name}")

        self._create_sh_files(stage=stage)
        self._submit(stage)
