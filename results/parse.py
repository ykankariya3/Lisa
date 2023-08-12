import os
import pandas as pd

def parse_out_file(file_path, df=None, prefix=""):
    with open(file_path, "r") as f:
        lines = f.readlines()

    current_file = None
    decomp_runtime = None
    decomposition = None
    runtime = None
    min_states = None

    for line in lines:
        line = line.strip()
        if line.startswith("AAAI2020-benchmarks-cap/"):
            if current_file:
                # Extract the folder name from the file path
                folder_name = os.path.basename(os.path.dirname(current_file))
                file_data = {
                    "FolderName": folder_name,
                    "FileName": os.path.basename(current_file),
                    f"{prefix}Decomp Runtime": float(decomp_runtime) if decomp_runtime is not None else None,
                    f"{prefix}Decomposition": int(decomposition) if decomposition is not None else None,
                    f"{prefix}Runtime": float(runtime) if runtime is not None else None,
                    f"{prefix}Min States": int(min_states) if min_states is not None else None,
                }

                if df is None:
                    df = pd.DataFrame(columns=["FolderName", "FileName"])

                # If the DataFrame already contains a row with the same FolderName and FileName, update it with the new data
                existing_row = df[(df["FolderName"] == folder_name) & (df["FileName"] == file_data["FileName"])]
                if not existing_row.empty:
                    idx = existing_row.index[0]
                    df.loc[idx, file_data.keys()] = file_data.values()
                else:
                    df = df.append(file_data, ignore_index=True)

            current_file = line
            decomp_runtime = None
            decomposition = None
            runtime = None
            min_states = None

        elif line.startswith("Decomp Runtime:"):
            decomp_runtime_str = line.split(":")[1].strip()
            decomp_runtime = float(decomp_runtime_str.split("ms")[0])
        elif line.startswith("Breakdown:"):
            decomposition = int(line.split(":")[1])
        elif line.startswith("Runtime:"):
            runtime_str = line.split(":")[1].strip()
            runtime = float(runtime_str.split("ms")[0])
        elif line.startswith("Min States:"):
            min_states = int(line.split(":")[1])

    # Store the data for the last benchmark
    if current_file:
        folder_name = os.path.basename(os.path.dirname(current_file))
        file_data = {
            "FolderName": folder_name,
            "FileName": os.path.basename(current_file),
            f"{prefix}Decomp Runtime": float(decomp_runtime) if decomp_runtime is not None else None,
            f"{prefix}Decomposition": int(decomposition) if decomposition is not None else None,
            f"{prefix}Runtime": float(runtime) if runtime is not None else None,
            f"{prefix}Min States": int(min_states) if min_states is not None else None,
        }

        if df is None:
            df = pd.DataFrame(columns=["FolderName", "FileName"])

        # If the DataFrame already contains a row with the same FolderName and FileName, update it with the new data
        existing_row = df[(df["FolderName"] == folder_name) & (df["FileName"] == file_data["FileName"])]
        if not existing_row.empty:
            idx = existing_row.index[0]
            df.loc[idx, file_data.keys()] = file_data.values()
        else:
            df = df.append(file_data, ignore_index=True)

    return df


if __name__ == "__main__":
    output_file = "caseBenchmark.xlsx"
    input_file = "Rport-2746820.out"
    sheet_name = "countResults"

    try:
        existing_data = pd.read_excel(output_file, sheet_name=sheet_name)
    except FileNotFoundError:
        existing_data = None

    new_data = parse_out_file(input_file, existing_data, "Lisa1Exp")
    new_data.to_excel(output_file, sheet_name=sheet_name, index=False)
