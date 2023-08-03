import os
import csv
import subprocess
import re
import pandas as pd
import shutil

# Paths to lisa1 and lisa2
lisa1_path = "src/lisa1/lisa"
lisa2_path = "src/lisa2/lisa"

# Command templates
lisa1_command = "{path} -exp -ltlf {file}"
lisa2_command = "{path} -ltlf {file}"

def extract_data(output):
    stdout_str = output.stdout
    runtime_match = re.search(r"Runtime: (\d+\.\d+)", stdout_str)
    states_match = re.search(r"Min States: (\d+)", stdout_str)
    breakdown_match = re.search(r"Breakdown: (\d+)", stdout_str)

    runtime = runtime_match.group(1) if runtime_match else None
    states = states_match.group(1) if states_match else None
    breakdown = breakdown_match.group(1) if breakdown_match else None
    return breakdown, runtime, states

# Path to the benchmarks folder
benchmark_folder = "quick"

# Initialize variables for total runtime
total_runtime_lisa1 = 0.0
total_runtime_lisa2 = 0.0

# Iterate over each folder in the benchmark directory
for folder in os.listdir(benchmark_folder):
    folder_path = os.path.join(benchmark_folder, folder)
    if os.path.isdir(folder_path):
        # CSV file path for the current folder
        print(folder)
        csv_file = os.path.join(folder_path, "benchmark_data.csv")

        # Collect the benchmark files
        benchmark_files = []
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                if not (file.startswith('._')) and file.endswith(".ltlf"):
                    benchmark_files.append(os.path.join(root, file))
                elif not (file.startswith('._')) and file.endswith(".tlsf"):
                    tlsf_file = os.path.join(root, file)
                    ltlf_file = os.path.splitext(file)[0] + ".ltlf"
                    ltlf_file_path = os.path.join(root, ltlf_file)

                    # Create new file with the same name but .ltlf extension
                    if not os.path.exists(ltlf_file_path):
                        # Use syfco to populate the new ltlf file
                        syfco_command = f"syfco {tlsf_file} --format ltlxba-fin"
                        syfco_output = subprocess.run(
                            syfco_command, shell=True, capture_output=True, text=True
                        )
                        ltlf_content = syfco_output.stdout
                        with open(ltlf_file_path, 'w') as ltlf_file:
                            ltlf_file.write(ltlf_content)
                        benchmark_files.append(ltlf_file_path)

        # Execute commands and collect data
        data = []
        for file in benchmark_files:
            file_name = os.path.basename(file)
            row = [file_name]

            print(file)

            # Run command in lisa1
            lisa1_cmd = lisa1_command.format(path=lisa1_path, file=file)
            lisa1_output = subprocess.run(
                lisa1_cmd, shell=True, capture_output=True, text=True
            )

            # Extract runtime, number of states, and number of subformulas from lisa1 output
            # Modify the code below based on the actual output format of lisa1
            lisa1_subformulas, lisa1_runtime, lisa1_states = extract_data(lisa1_output)  # Replace with actual number of subformulas

            row.extend([lisa1_runtime, lisa1_states, lisa1_subformulas])

            #total_runtime_lisa1 += float(lisa1_runtime)

            # Run command in lisa2
            lisa2_cmd = lisa2_command.format(path=lisa2_path, file=file)
            lisa2_output = subprocess.run(
                lisa2_cmd, shell=True, capture_output=True, text=True
            )

            # Extract runtime, number of states, and number of subformulas from lisa2 output
            # Modify the code below based on the actual output format of lisa2
            lisa2_subformulas, lisa2_runtime, lisa2_states = extract_data(lisa2_output)

            row.extend([lisa2_runtime, lisa2_states, lisa2_subformulas])

            #total_runtime_lisa2 += float(lisa2_runtime)

            breakdown_diff = int(lisa2_subformulas) - int(lisa1_subformulas)
            states_equal = lisa1_states == lisa2_states

            row.extend([breakdown_diff, states_equal])

            data.append(row)

        #total_runtime_row = ["Total Runtime", total_runtime_lisa1, "", "", total_runtime_lisa2, "", ""]
        #data.append(total_runtime_row)

        # Write data to CSV file
        with open(csv_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Benchmark", "Lisa1 Runtime", "Lisa1 States", "Lisa1 Subformulas",
                             "Lisa2 Runtime", "Lisa2 States", "Lisa2 Subformulas", 
                             "Breakdown Difference",  "Misstates Equal"])
            writer.writerows(data)

        print("Data collection completed for folder", folder)

        # Convert the CSV file to Excel
        excel_file = os.path.join(folder_path, f"{folder}.xlsx")
        df = pd.read_csv(csv_file)
        df.to_excel(excel_file, index=False)

        print("Excel file created for folder", folder)

print("All data collection and conversion completed.")
