import sys

# Read input from stdin
data = [line.strip().split("\t") for line in sys.stdin]

# Convert second column to integers
for row in data:
    row[1] = int(row[1])

# Find the smallest value in column 2
min_value = min(row[1] for row in data)

# Compute scaling factors and print results
for row in data:
    scaling_factor = min_value / row[1]
    print(f"{row[0]}\t{row[1]}\t{scaling_factor:.6f}")
