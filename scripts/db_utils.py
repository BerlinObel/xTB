import sys
import os
from time import sleep
import sqlite3
from pickle import dumps, loads
import pandas as pd
from contextlib import contextmanager

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

def insert_data(data):
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
        INSERT INTO MoleculeData ({columns})
        VALUES ({placeholders})
    """, values)
    
    conn.commit()
    close_db(conn)
    
def set_calculation_stage(wanted_stage):
    conn = connect_db()
    cursor = conn.cursor()
    cursor.execute(f"UPDATE MoleculeData SET CalculationStage = '{wanted_stage}'")
    conn.commit()
    close_db(conn)
    
def retrieve_data(filters):
    """Retrieve data from the database based on the filters."""
    conn = connect_db()
    cursor = conn.cursor()
    
    # Modify the query construction logic to use LIKE for 'CalculationStage'
    query = "SELECT * FROM MoleculeData WHERE "
    query += " AND ".join([f"{key} LIKE ?" if key == 'CalculationStage' else f"{key} = ?" for key in filters.keys()])
    
    print("Query:", query)
    print("Parameters:", tuple(filters.values()))
    print("Database Path:", DB_PATH)
    print("Current Working Directory:", os.getcwd())
    print(conn)
    
    # Adjust the parameters to use % for 'CalculationStage'
    params = [f"%{value}%" if key == 'CalculationStage' else value for key, value in filters.items()]
    
    try:
        cursor.execute(query, tuple(params))
    except Exception as e:
        print(f"An error occurred: {e}")
        close_db(conn)
        return None
    
    data = cursor.fetchall()
    
    # Get the column names from the cursor description
    columns = [desc[0] for desc in cursor.description]
    
    # Create a pandas DataFrame from the data
    data_df = pd.DataFrame(data, columns=columns)
    
    # Deserialize complex object columns from binary format
    for col in ['ProductObject', 'ReactantObject', 'TransitionStateObject', 'GroundStateStats', 'TransistionStateStats', 'ExcitationStats']:
        if col in data_df.columns:
            data_df[col] = data_df[col].apply(lambda x: deserialize_object(x) if x is not None else x)
    
    close_db(conn)
    
    return data_df    


def update_data(data_df):
    """Update data in the database with a DataFrame based on batch ID"""
    
    # Serialize complex object columns to binary format
    for col in ['ProductObject', 'ReactantObject', 'TransitionStateObject', 'GroundStateStats', 'TransistionStateStats', 'ExcitationStats']:
        if col in data_df.columns:
            data_df[col] = data_df[col].apply(lambda x: serialize_object(x) if x is not None else x)

    for _ in range(5):  # Retry up to 5 times
        try:
            with locked_db_connection() as conn:
                # Write DataFrame to SQLite
                data_df.to_sql('NewTable', conn, if_exists='replace', index=False)

                # Construct the SQL update query dynamically with conditional logic for 'CalculationStage'
                set_statements = ", ".join(
                    f"{col} = CASE WHEN '{col}' = 'CalculationStage' THEN (CASE WHEN CalculationStage = 'storage' THEN (SELECT {col} FROM NewTable WHERE NewTable.HashedName = MoleculeData.HashedName) WHEN instr(CalculationStage, (SELECT {col} FROM NewTable WHERE NewTable.HashedName = MoleculeData.HashedName)) > 0 THEN CalculationStage ELSE CalculationStage || ', ' || (SELECT {col} FROM NewTable WHERE NewTable.HashedName = MoleculeData.HashedName) END) ELSE (SELECT {col} FROM NewTable WHERE NewTable.HashedName = MoleculeData.HashedName) END"
                    for col in data_df.columns if col not in ['HashedName']
                )

                query = f'''
                UPDATE MoleculeData
                SET 
                    {set_statements}
                WHERE 
                    EXISTS (SELECT 1 FROM NewTable WHERE NewTable.HashedName = MoleculeData.HashedName);
                '''
                
                # Execute the query and drop the temporary table within the context manager
                conn.execute(query)
                conn.execute("DROP TABLE NewTable")
                break  # Break out of the loop if the operation succeeds
        except sqlite3.OperationalError as e:
            if "database is locked" in str(e):
                print("Database is locked, retrying...")
                sleep(1)  # Wait for 1 second before retrying
            else:
                raise  # Re-raise the exception if it's not a "database is locked" error

    # Close the connection
    conn.close()




