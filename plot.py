import pandas as pd
import matplotlib.pyplot as plt

# Function to calculate cumulative times
def calculate_cumulative_times(run_times):
    return [sum(run_times[:i + 1]) for i in range(len(run_times))]

# Read data from the Excel files for solution1 and solution2
df_solution1 = pd.read_excel('CounterDockerResults.xlsx')
# df_solution2 = pd.read_excel('nim1.xlsx')
# df_solution3 = pd.read_excel('counterLydia.xlsx')
# df_solution4 = pd.read_excel('dcounterLydia.xlsx')

# Sort the data based on runtimes
# df_solution1 = pd.concat([df_solution1.sort_values(by='Benchmark'), df_solution2.sort_values(by='Benchmark')])
# df_solution2 = pd.concat([df_solution1.sort_values(by='Lisa1 Decomp TimeLisa2 Runtime'), df_solution2.sort_values(by='Lisa1 Decomp TimeLisa2 Runtime')])
# df_solution3 = pd.concat([df_solution3.sort_values(by='Lisa1 Runtime'), df_solution4.sort_values(by='Lisa1 Runtime')])

df_solution1 = df_solution1.sort_values(by='Lydia')
df_solution2 = df_solution1.sort_values(by='Lisa2')
# df_solution3 = df_solution1.sort_values(by='Lydia Runtime')
# df_solution4 = df_solution1.sort_values(by='Lisa2 inc/ Map Runtime')


# Extract the runtimes column from the DataFrames
run_times_solution1 = df_solution1['Lydia'].tolist()
run_times_solution2 = df_solution2['Lisa2'].tolist()
# run_times_solution3 = df_solution3['Lydia Runtime'].tolist()
# run_times_solution4 = df_solution4['Lisa2 inc/ Map Runtime'].tolist()

threshold = 600000

run_times_solution1 = [time for time in run_times_solution1 if time <= threshold and time > 0]
run_times_solution2 = [time for time in run_times_solution2 if time <= threshold and time > 0]
# run_times_solution3 = [time for time in run_times_solution3 if time <= threshold and time > 0]
# run_times_solution4 = [time for time in run_times_solution4 if time <= threshold and time > 0]

# Create the plot
num_benchmarks_solution1 = range(1, len(run_times_solution1) + 1)
num_benchmarks_solution2 = range(1, len(run_times_solution2) + 1)
# num_benchmarks_solution3 = range(1, len(run_times_solution3) + 1)
# num_benchmarks_solution4= range(1, len(run_times_solution4) + 1)

plt.plot(num_benchmarks_solution1, run_times_solution1, marker='o', linestyle='-', color='b', label='Lydia')
plt.plot(num_benchmarks_solution2, run_times_solution2, marker='o', linestyle='-', color='r', label='Lisa2')
# plt.plot(num_benchmarks_solution3, run_times_solution3, marker='o', linestyle='-', color='g', label='Lydia')
# plt.plot(num_benchmarks_solution4, run_times_solution4, marker='o', linestyle='-', color='y', label='Lisa2 Mapping Runtime')


plt.xlabel('Number of Benchmarks Solved')
plt.ylabel('Time')
plt.title('CounterDocker')
plt.legend()
plt.grid(True)

# Display the plot
plt.show()
