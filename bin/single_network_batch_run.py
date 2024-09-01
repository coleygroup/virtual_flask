from shared.hpc.batch import SingleNetworkDistributedPipeline
from rdkit import RDLogger
import sys
from shared.hpc.reset_tables import cancel_all_transactions, reset_psql
RDLogger.DisableLog("rdApp.*")

x = SingleNetworkDistributedPipeline("t", [1,4,5,6], 4)

x.pipe("g16")
