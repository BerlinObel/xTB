import sys
import os
from time import sleep, time
import sqlite3
from pickle import dumps, loads
import pandas as pd
from contextlib import contextmanager
import random

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from settings import DB_PATH

def connect_db():
    """Establish connection to the database."""
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA busy_timeout = 5000")  # Set a busy timeout of 5000 milliseconds (5 seconds)
    return conn

def close_db(conn):
    """Close connection to the database."""
    conn.close()

def serialize_object(obj):
    """Serialize an object to be stored in the database."""
    return dumps(obj)

def deserialize_object(blob):
    """Deserialize an object retrieved from the database."""
    return loads(blob)



@contextmanager
def locked_db_connection():
    conn = connect_db()
    try:
        with conn:
            yield conn
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        close_db(conn)
        
def create_table():
    """Create the table if it doesn't already exist."""
    conn = connect_db()
    cursor = conn.cursor()
    
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS MoleculeData (
            ID INTEGER PRIMARY KEY AUTOINCREMENT,
            HashedName TEXT,
            Smiles TEXT,
            Multiplicity REAL,
            Charge REAL,
            BatchID TEXT,
            CalculationStage TEXT,
            ProductObject BLOB,
            ReactantObject BLOB,
            TransitionStateObject BLOB,
            StorageEnergy REAL,
            BackReactionBarrier REAL,
            MaxAbsorptionProduct REAL,
            MaxOscillatorStrengthProduct REAL,
            MaxAbsorptionReactant REAL,
            MaxOscillatorStrengthReactant REAL,
            SolarConversionEfficiency REAL,
            GroundStateStats BLOB,
            TransistionStateStats BLOB,
            ExcitationStats BLOB
        )
    """)
    
    close_db(conn)
    
def create_orca_table():
    """Create the table if it doesn't already exist."""
    conn = connect_db()
    cursor = conn.cursor()
    
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS ORCAData (
            ID INTEGER PRIMARY KEY AUTOINCREMENT,
            HashedName TEXT,
            BatchID TEXT,
            CalculationStage TEXT,
            StorageEnergy REAL,
            BackReactionBarrier REAL,
            MaxAbsorptionProduct REAL,
            MaxOscillatorStrengthProduct REAL,
            MaxAbsorptionReactant REAL,
            MaxOscillatorStrengthReactant REAL,
            SolarConversionEfficiency REAL
        )
    """)
    
    close_db(conn)
    
def insert_data(data, db_table="MoleculeData"):
    """Insert data into the database."""
    conn = connect_db()
    cursor = conn.cursor()
    
    # Get the keys and values from the data
    keys = data.keys()
    values = [dumps(data[key]) if isinstance(data[key], (dict, list)) else data[key] for key in keys]
    
    # Construct the SQL query dynamically
    columns = ', '.join(keys)
    placeholders = ', '.join(['?' for _ in keys])
    
    cursor.execute(f"""
        INSERT INTO {db_table} ({columns})
        VALUES ({placeholders})
    """, values)
    
    conn.commit()
    close_db(conn)

def add_batch(input_file, batch_id):
    chunk_size = int(input("Batch size: "))
    data = pd.read_csv(input_file)
    chunks = [data[i:i + chunk_size] for i in range(0, data.shape[0], chunk_size)]

    for chunk in chunks:
        for _, row in chunk.iterrows():
            molecule_data = {
                'HashedName': row['comp_name'],
                'Smiles': row['smiles'],
                'Multiplicity': row['multiplicity'],
                'Charge': row['charge'],
                'BatchID': batch_id,
                'CalculationStage': 'storage'
            }
            insert_data(molecule_data)

def set_calculation_stage(wanted_stage):
    conn = connect_db()
    cursor = conn.cursor()
    cursor.execute(f"UPDATE MoleculeData SET CalculationStage = '{wanted_stage}'")
    conn.commit()
    close_db(conn)
    
def delete_table():
    """Delete a table from an SQLite database."""
    table_name = input("Delete Table: ")
    conn = None
    try:
        conn = connect_db()
        cursor = conn.cursor()
        cursor.execute(f"DROP TABLE {table_name}")
        conn.commit()
        print(f"Table {table_name} deleted successfully.")
    except sqlite3.Error as e:
        print(f"An error occurred: {e}")
    finally:
        if conn:
            close_db(conn)
  

def rename_table():
    """Rename a table in an SQLite database."""
    old_table_name = input("Old Table Name: ")
    new_table_name = input("New Table Name: ")
    conn = None
    try:
        conn = connect_db()
        cursor = conn.cursor()
        cursor.execute(f"ALTER TABLE {old_table_name} RENAME TO {new_table_name}")
        conn.commit()
        print(f"Table {old_table_name} renamed to {new_table_name} successfully.")
    except sqlite3.Error as e:
        print(f"An error occurred: {e}")
    finally:
        if conn:
            close_db(conn)
            
def backup_table(original_table_name, backup_table_name):
    """Create a backup of a table by copying it to a new table with a different name."""
    conn = None
    try:
        # Connect to the database
        conn = connect_db()

        # Create a cursor object to execute SQL queries
        cursor = conn.cursor()

        # Check if the backup table already exists
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (backup_table_name,))
        if cursor.fetchone():
            print(f"Backup table {backup_table_name} already exists.")
            return

        # Execute the query to create a copy of the table
        cursor.execute(f"CREATE TABLE {backup_table_name} AS SELECT * FROM {original_table_name}")

        # Commit the changes
        conn.commit()

        print(f"Table {original_table_name} backed up as {backup_table_name} successfully.")

    except sqlite3.Error as e:
        print(f"An error occurred: {e}")

    finally:
        # Close the database connection
        if conn:
            close_db(conn)

def retrieve_data(filters, db_table='MoleculeData'):
    """
    Retrieve data from the database based on the filters.

    Parameters:
    filters (dict): A dictionary where keys are column names and values are filter criteria.
    db_table (str): The database table to query. Defaults to 'MoleculeData'.

    Returns:
    pandas.DataFrame: A DataFrame containing the retrieved data, or None if an error occurs.
    """
    
    conn = None
    try:
        conn = connect_db()
        cursor = conn.cursor()
        query_parts = []
        params = []

        for key, value in filters.items():
            if isinstance(value, str) and "%" in value:
                query_parts.append(f"{key} LIKE ?")
                params.append(value)
            else:
                query_parts.append(f"{key} = ?")
                params.append(value)

        query = f"SELECT * FROM {db_table} WHERE " + " AND ".join(query_parts)
        
        # Debugging prints (you can remove these in production)
        print("Query:", query)
        print("Parameters:", params)
        print("Database Path:", DB_PATH)
        print("Current Working Directory:", os.getcwd())

        cursor.execute(query, params)
        data = cursor.fetchall()
        columns = [desc[0] for desc in cursor.description]
        data_df = pd.DataFrame(data, columns=columns)
        for col in ['ProductObject', 'ReactantObject', 'TransitionStateObject', 'GroundStateStats', 'TransistionStateStats', 'ExcitationStats']:
            if col in data_df.columns:
                data_df[col] = data_df[col].apply(lambda x: deserialize_object(x) if x is not None else x)
        return data_df
    except sqlite3.Error as e:
        print(f"An error occurred: {e}")
        return None
    finally:
        if conn:
            close_db(conn)


def update_data(data_df, MoleculeData='MoleculeData'):
    """Update data in the database with a DataFrame based on batch ID"""
    
    # Serialize complex object columns to binary format
    for col in ['ProductObject', 'ReactantObject', 'TransitionStateObject', 'GroundStateStats', 'TransistionStateStats', 'ExcitationStats']:
        if col in data_df.columns:
            data_df[col] = data_df[col].apply(lambda x: serialize_object(x) if x is not None else x)

    for _ in range(5):  # Retry up to 5 times
        try:
            with locked_db_connection() as conn:
                random_integer = random.randint(1, 1000)
                unique_table_name = f"NewTable_{int(time())}_{random_integer}"

                # Write DataFrame to SQLite
                data_df.to_sql(unique_table_name, conn, if_exists='replace', index=False)

                # Construct the SQL update query dynamically with conditional logic for 'CalculationStage'
                set_statements = ", ".join(
                    f"{col} = CASE WHEN '{col}' = 'CalculationStage' THEN (CASE WHEN CalculationStage = 'storage' THEN (SELECT {col} FROM {unique_table_name} WHERE {unique_table_name}.HashedName = {MoleculeData}.HashedName) WHEN instr(CalculationStage, (SELECT {col} FROM {unique_table_name} WHERE {unique_table_name}.HashedName = {MoleculeData}.HashedName)) > 0 THEN CalculationStage ELSE CalculationStage || ', ' || (SELECT {col} FROM {unique_table_name} WHERE {unique_table_name}.HashedName = {MoleculeData}.HashedName) END) ELSE (SELECT {col} FROM {unique_table_name} WHERE {unique_table_name}.HashedName = {MoleculeData}.HashedName) END"
                    for col in data_df.columns if col not in ['HashedName']
                )

                query = f'''
                UPDATE {MoleculeData}
                SET 
                    {set_statements}
                WHERE 
                    EXISTS (SELECT 1 FROM {unique_table_name} WHERE {unique_table_name}.HashedName = {MoleculeData}.HashedName);
                '''
                
                # Execute the query and drop the temporary table within the context manager
                conn.execute(query)
                conn.execute(f"DROP TABLE {unique_table_name}")
                break  # Break out of the loop if the operation succeeds
        except sqlite3.OperationalError as e:
            if "database is locked" in str(e):
                print("Database is locked, retrying...")
                sleep(100)  # Wait for 100 seconds before retrying
            else:
                raise  # Re-raise the exception if it's not a "database is locked" error

    conn.close()




