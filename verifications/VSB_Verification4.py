# ===== VSB Verifications: Velocity Profiles (Quiver plotting) =====
import matplotlib.pyplot as plt
import numpy as np
print("=== BEGINNING Verifications4.py ===")

# === Imports ===
print("   IMPORTING PACKAGES...")

# === Helper function ===


def toFloat_list(list):
    res = []
    for item in list:
        res.append(float(item))
    return res


# === Read in data ===
print("   READING IN FILES...")
filenames = ["X_mesh.txt", "Y_mesh.txt", "UX_mesh.txt", "UY_mesh.txt"]

file1 = open(filenames[0], 'r')
file2 = open(filenames[1], 'r')
file3 = open(filenames[2], 'r')
file4 = open(filenames[3], 'r')

X_mesh = []
Y_mesh = []
UX = []
UY = []

for line in file1:
    values = line.split(',')[:-1]
    X_mesh.append(toFloat_list(values))

for line in file2:
    values = line.split(',')[:-1]
    Y_mesh.append(toFloat_list(values))

for line in file3:
    values = line.split(',')[:-1]
    UX.append(toFloat_list(values))

for line in file4:
    values = line.split(',')[:-1]
    UY.append(toFloat_list(values))

# === Extract velocity profile points ===
X_left = []
X_right = []
for i in np.arange(1, 25, 1):
    X_left.append([X_mesh[12][12], Y_mesh[i][i]])
    X_right.append([X_mesh[12][12] + UY[i][i], Y_mesh[i][i]])

# === Plot data with quiver (b/c Julia can't)... ===
print("   PLOTTING...")
fig1, ax1 = plt.subplots()
ax1.set_title('Velocity Profile')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_xlim(0, 10)
ax1.set_ylim(-0.5, 5.5)
ax1.grid(True)
for i in np.arange(1, 25, 2):
    ax1.scatter(X_mesh[12][12], Y_mesh[i][i], marker='o', color='b')
    ax1.scatter(X_mesh[12][12] + UY[i][i], Y_mesh[i][i], marker='o', color='r')
Q = ax1.quiver(X_mesh[::4], Y_mesh[::4], UX[::4], UY[::4])
# QK = ax1.quiverkey(Q, )
plt.savefig('VSB_Velocity_Profile.pdf')

print("=== END OF Verifications4.py ===")
