import sys
import pandas as pd
import matplotlib.pyplot as plt

if len(sys.argv) < 3:
    print("Usage: python plot_anneal.py input.csv output.png")
    sys.exit(1)

csvfile = sys.argv[1]
pngfile = sys.argv[2]

# read data
data = pd.read_csv(csvfile)

# plot
plt.figure(figsize=(8,6))
plt.plot(data["T"], data["best_km"], label="Best distance")
plt.plot(data["T"], data["current_km"], alpha=0.4, label="Current distance")

plt.xscale("log")
plt.xlabel("Temperature")
plt.ylabel("Total Distance (km)")
plt.title("Simulated Annealing Schedule")
plt.legend()
plt.grid(True)

plt.savefig(pngfile, dpi=200)
plt.close()

