import matplotlib.pyplot as plt

enzymes = {
    "EcoRI": ("GAATTC", 1),
    "BamHI": ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
    "TaqI": ("TCGA", 1),
    "HaeIII": ("GGCC", 2),
}

seq = "".join(line.strip() for line in open("sequence.fasta") if not line.startswith(">")).upper()
print("Sequence length:", len(seq), "\n")

lanes = {"Ladder": [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500, 250, 100]}

for name, (site, cut) in enzymes.items():
    cuts = [i + cut for i in range(len(seq) - len(site) + 1) if seq[i:i + len(site)] == site]
    bounds = [0] + cuts + [len(seq)]
    frags = [bounds[i + 1] - bounds[i] for i in range(len(bounds) - 1)]
    lanes[name] = frags
    print(name, "site:", site, "cut after base", cut)
    print("Number of cleavages:", len(cuts))
    print("Cleavage positions:", cuts)
    print("Fragment lengths:", frags, "\n")

plt.figure(figsize=(8, 6), facecolor="black")
plt.gca().set_facecolor("black")
for x, frags in enumerate(lanes.values()):
    for f in frags:
        plt.hlines(f, x - 0.3, x + 0.3, color="white", linewidth=2)
plt.yscale("log")
plt.ylim(50, 15000)
plt.xticks(range(len(lanes)), lanes.keys(), color="white")
plt.yticks(lanes["Ladder"], lanes["Ladder"], color="white")
plt.title("Electrophoresis gel simulation", color="white")
plt.show()
