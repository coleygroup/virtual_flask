# import docker
import psycopg2
import time


def connect_to_local_psql():
    # conn = psycopg2.connect(host="3.138.199.28", port=5432, database="postgres", user="username", password="password")
    conn = psycopg2.connect(host="localhost", port=5432, database="postgres", user="bmahjour", password="postgres")
    cur = conn.cursor()
    return conn, cur

def cancel_all_transactions():
    try:
        # Connect to the PostgreSQL database
        conn, cur = connect_to_local_psql()
        conn.autocommit = True  # Required to run admin queries without transactions

        # Execute the query to terminate all active transactions except the current session
        cur.execute(
            """
            SELECT pg_terminate_backend(pid)
            FROM pg_stat_activity
            WHERE state = 'active' 
            AND pid <> pg_backend_pid()  
            AND usename != 'postgres';   
        """
        )

        print("All active transactions have been cancelled.")

    except Exception as e:
        print(f"An error occurred: {e}")

    finally:
        # Clean up the connection and cursor
        if cur:
            cur.close()
        if conn:
            conn.close()


def reset_psql(conn, name):

    cur = conn.cursor()

    table_name = f"{name}_networks"
    cur.execute(f"DROP TABLE IF EXISTS {name}_networks CASCADE;")
    print(3)

    create_table_command = f"""
    CREATE TABLE {table_name} (
        network_id SERIAL PRIMARY KEY,
        mapped_input_smiles_list TEXT[],
        reagents TEXT[],
        input_mapped_smiles TEXT,
        input_unmapped_smiles TEXT,
        starting_material_atom_maps JSONB,
        runtime JSONB,
        cum_runtime JSONB,
        mech_count_map JSONB,
        end_state TEXT,
        other_data JSONB
    );
    """

    cur.execute(create_table_command)
    print(4)

    table_name = f"{name}_nodes"
    cur.execute(f"DROP TABLE IF EXISTS {table_name} CASCADE;")
    print(5)

    create_table_command = f"""
    CREATE TABLE {table_name} (
        node_id SERIAL PRIMARY KEY,
        network_id int REFERENCES {name}_networks(network_id),
        mapped_smiles TEXT,
        unmapped_smiles TEXT,
        these_reacting_atoms_path JSONB,
        product_in_precalc BOOLEAN,
        other_data JSONB
    );
    """

    cur.execute(create_table_command)
    print(6)

    table_name = f"{name}_edges"
    cur.execute(f"DROP TABLE IF EXISTS {table_name} CASCADE;")
    print(7)

    create_table_command = f"""
    CREATE TABLE {table_name} (
        edge_id SERIAL PRIMARY KEY,
        network_id int REFERENCES {name}_networks(network_id),
        source_node_id int REFERENCES {name}_nodes(node_id),
        destination_node_id int REFERENCES {name}_nodes(node_id),
        other_data JSONB
    );
    """

    cur.execute(create_table_command)
    print(8)

    conn.commit()

    cur.close()
    conn.close()
