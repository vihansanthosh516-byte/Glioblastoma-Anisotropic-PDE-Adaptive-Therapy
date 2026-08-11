from idc_index import IDCClient

# Initialize the IDC client
client = IDCClient()

# 1. Print all available collections to find the exact string match for MU-Glioma-Post
collections = client.get_collections()
mu_collections = [c for c in collections if "glioma" in c.lower() or "mu" in c.lower()]
print("Matching collections found in IDC:", mu_collections)

# Use the exact string name found from the list above (e.g., 'mu_glioma_post' or similar)
if mu_collections:
    target_collection = mu_collections[0]
    print(f"\nQuerying collection: {target_collection}")
    
    # Query patients using the verified collection name
    query = f"""
        SELECT PatientID, COUNT(DISTINCT StudyInstanceUID) as study_count
        FROM index 
        WHERE collection_id = '{target_collection}'
        GROUP BY PatientID
        HAVING study_count >= 2
        LIMIT 3
    """
    results_df = client.sql_query(query)
    print("Longitudinal patients found:")
    print(results_df)

    if not results_df.empty:
        target_patient = results_df.iloc[0]["PatientID"]
        print(f"\nDownloading imaging series for patient: {target_patient}...")
        client.download_from_selection(patientId=target_patient, downloadDir="./data/tcia")
        print(f"Successfully downloaded data for {target_patient} into ./data/tcia/")
else:
    print("Collection filter keyword not matched. Check available collection names.")