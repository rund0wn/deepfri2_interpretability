import pickle

def check_uniprot_id_exists(file_path, target_uniprot_id):
    """
    Reads a pickle file containing a list of dictionaries and checks 
    if a specific uniprot_id exists in any of the records.
    """
    try:
        with open(file_path, 'rb') as f:
            records = pickle.load(f)
            
        # Iterate through the list of dictionaries
        for record in records:
            # Use .get() to avoid KeyError if 'uniprot_id' is missing in some rows
            if record.get('uniprot_id') == target_uniprot_id:
                print(f"UniProt ID '{target_uniprot_id}' is in the dataset.")
                return True
                
        # If the loop finishes without returning, the ID wasn't found
        print(f"UniProt ID '{target_uniprot_id}' is missing from the dataset.")
        return False

    except FileNotFoundError:
        print(f"Error: Could not find the file at '{file_path}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Specify your file path and the UniProt ID you want to search for
    FILE_PATH = 'GPSite_processed/gpsite_binding_records_lt1020.pkl'
    
    # Replace 'P12345' with the actual UniProt ID you are looking for
    TARGET_ID = 'P13569' 
    
    check_uniprot_id_exists(FILE_PATH, TARGET_ID)