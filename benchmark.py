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
lisa1_command = "timeout 600 {path} -exp -ltlf {file}"
lisa2_command = "timeout 600 {path} -ltlf {file}"
lydia_command = "timeout 600 lydia -l ltlf -f {file}"
lisa2_docker_command = "timeout 600 docker run -v /Volumes/STORAGE/LISA/lisa:/app/test lisa-app ./lisa -exp -ltlf /app/test/{file}"

def extract_data(output):
    stdout_str = output.stdout
    print(stdout_str)
    runtime_match = re.search(r"Runtime: (\d+)", stdout_str)
    states_match = re.search(r"Min States: (\d+)", stdout_str)
    breakdown_match = re.search(r"Breakdown: (\d+)", stdout_str)
    breakdown_runtime = re.search(r"Decomp: (\d+)", stdout_str)
    unique_match = re.search(r"Unique Constructs: (\d+)", stdout_str)
    repetitive_match = re.search(r"Repeated Constructs: (\d+)", stdout_str)

    runtime = runtime_match.group(1) if runtime_match else None
    states = states_match.group(1) if states_match else None
    breakdown = breakdown_match.group(1) if breakdown_match else None
    breakdowntime = breakdown_runtime.group(1) if breakdown_runtime else None
    unique = unique_match.group(1) if unique_match else None
    repetitive = repetitive_match.group(1) if repetitive_match else None
    return breakdown, runtime, states, breakdowntime, unique, repetitive

# Path to the benchmarks folder
benchmark_folder = "large"

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
            # lisa1_cmd = lisa1_command.format(path=lisa1_path, file=file)
            # lisa1_output = subprocess.run(
            #     lisa1_cmd, shell=True, capture_output=True, text=True
            # )

            # # Extract runtime, number of states, and number of subformulas from lisa1 output
            # # Modify the code below based on the actual output format of lisa1
            # lisa1_subformulas, lisa1_runtime, lisa1_states, lisa1_decomp = extract_data(lisa1_output)  # Replace with actual number of subformulas

            # row.extend([lisa1_runtime, lisa1_states, lisa1_subformulas, lisa1_decomp])


            #total_runtime_lisa1 += float(lisa1_runtime)

            # Run command in lisa2
            # lisa2_cmd = lisa2_command.format(path=lisa2_path, file=file)
            # lisa2_output = subprocess.run(
            #     lisa2_cmd, shell=True, capture_output=True, text=True
            # )

            # # Extract runtime, number of states, and number of subformulas from lisa2 output
            # # Modify the code below based on the actual output format of lisa2
            # lisa2_subformulas, lisa2_runtime, lisa2_states, decomp, unique, repetitive = extract_data(lisa2_output)

            # row.extend([lisa2_runtime, lisa2_states, lisa2_subformulas, decomp, unique, repetitive])

            #total_runtime_lisa2 += float(lisa2_runtime)
            lydia_runtime = 0
            count = 0

            for i in range(1):

                lydia_cmd = lydia_command.format(file=file)
                lydia_output = subprocess.run(
                    lydia_cmd, shell=True, capture_output=True, text=True
                )

                temp = re.search(r"Overall time elapsed: (\d+\.\d+)", lydia_output.stdout)
                if temp:
                    lydia_runtime = temp.group(1) 
                    count = count + 1

            if (count == 0):
                count = 1
            print(float(lydia_runtime) / count)
            print()

            row.extend([float(lydia_runtime) / count])

            # #Run command in lisa2
            # lisa2_docker_cmd = lisa2_docker_command.format(file=file)
            # lisa2_output = subprocess.run(
            #     lisa2_docker_cmd, shell=True, capture_output=True, text=True
            # )

            # # Extract runtime, number of states, and number of subformulas from lisa2 output
            # # Modify the code below based on the actual output format of lisa2
            # lisa2_subformulas, lisa2_runtime, lisa2_states, decomp, unique, repetitive = extract_data(lisa2_output)

            # row.extend([lisa2_runtime, lisa2_states, lisa2_subformulas, decomp, unique, repetitive])

            data.append(row)

        #total_runtime_row = ["Total Runtime", total_runtime_lisa1, "", "", total_runtime_lisa2, "", ""]
        #data.append(total_runtime_row)

        # Write data to CSV file
        with open(csv_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Benchmark", "Runtime", "States", "Subformulas", "Decomp Time",
                             "Unique Constructions", "Repetitive Constructions"])
            writer.writerows(data)

        print("Data collection completed for folder", folder)

        # Convert the CSV file to Excel
        excel_file = os.path.join(folder_path, f"{folder}.xlsx")
        df = pd.read_csv(csv_file)
        df.to_excel(excel_file, index=False)

        print("Excel file created for folder", folder)

print("All data collection and conversion completed.")
