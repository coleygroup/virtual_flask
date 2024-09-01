from shared.hpc.batch import DistributedPipeline
from rdkit import RDLogger
import sys
from shared.hpc.reset_tables import cancel_all_transactions, reset_psql

# import psycopg2

RDLogger.DisableLog("rdApp.*")


server = sys.argv[1]
percent = float(sys.argv[2])
reset_tables = int(sys.argv[3])

cancel_all_transactions()

if reset_tables == 1:
    reset_psql()


x = DistributedPipeline(
    output_path="./example/distributed_test_1",
    data_path="../../shared/data",
    d=300,
    filter_seen=False,
)

x.load_data("inventory.xlsx")

x.pipe("all", clear_data=True, create_inputs=True, percent=percent)
