import pandas as pd
import matplotlib.pyplot as plt

# Function to calculate cumulative times
def calculate_cumulative_times(run_times):
    return [sum(run_times[:i + 1]) for i in range(len(run_times))]

# Read data from the Excel files for solution1 and solution2
df_solution1 = pd.read_excel('countBenchmark.xlsx')
df_solution2 = pd.read_excel('nimBenchmark.xlsx')
df_solution3 = pd.read_excel('caseBenchmark.xlsx')

# Sort the data based on runtimes
df_solution1 = pd.concat([df_solution1.sort_values(by='Lisa1Runtime'), df_solution2.sort_values(by='Lisa1Runtime'), df_solution3.sort_values(by='Lisa1Runtime')])
df_solution2 = pd.concat([df_solution1.sort_values(by='Lisa2Runtime'), df_solution2.sort_values(by='Lisa2Runtime'), df_solution3.sort_values(by='Lisa2Runtime')])
df_solution3 = pd.concat([df_solution1.sort_values(by='Lisa1ExpRuntime'), df_solution2.sort_values(by='Lisa1expRuntime'), df_solution3.sort_values(by='Lisa1ExpRuntime')])

df_solution1 = df_solution1.sort_values(by='Lisa1Runtime')
df_solution2 = df_solution1.sort_values(by='Lisa2Runtime')
df_solution3 = df_solution1.sort_values(by='Lisa1ExpRuntime')


# Extract the runtimes column from the DataFrames
run_times_solution1 = df_solution1['Lisa1Runtime'].tolist()
run_times_solution2 = df_solution2['Lisa2Runtime'].tolist()
run_times_solution3 = df_solution3['Lisa1ExpRuntime'].tolist()

threshold = 600000

run_times_solution1 = [time for time in run_times_solution1 if time <= threshold]
run_times_solution2 = [time for time in run_times_solution2 if time <= threshold]
run_times_solution3 = [time for time in run_times_solution3 if time <= threshold]

# Calculate the cumulative times for each solution
cumulative_times_solution1 = calculate_cumulative_times(run_times_solution1)
cumulative_times_solution2 = calculate_cumulative_times(run_times_solution2)
cumulative_times_solution3 = calculate_cumulative_times(run_times_solution3)

# Create the plot
num_benchmarks_solution1 = range(1, len(run_times_solution1) + 1)
num_benchmarks_solution2 = range(1, len(run_times_solution2) + 1)
num_benchmarks_solution3 = range(1, len(run_times_solution3) + 1)

plt.plot(num_benchmarks_solution1, cumulative_times_solution1, marker='o', linestyle='-', color='b', label='Lisa1')
plt.plot(num_benchmarks_solution2, cumulative_times_solution2, marker='o', linestyle='-', color='r', label='Lisa2')
plt.plot(num_benchmarks_solution3, cumulative_times_solution3, marker='o', linestyle='-', color='g', label='Lisa1 Sym')

plt.xlabel('Number of Benchmarks Solved')
plt.ylabel('Cumulative Time')
plt.title('Benchmark')
plt.legend()
plt.grid(True)

# Display the plot
plt.show()
