import os
import pandas as pd
import matplotlib.pyplot as plt

i = 77
run_id = f"{i:03d}"
run_path = os.path.join("runs", run_id)
trap_force_path = os.path.join(run_path, "force.csv")
img_path = os.path.join(run_path, "force.png")
df_trap_force = pd.read_csv(trap_force_path)


plt.figure(figsize=(6, 4))
plt.plot(df_trap_force["step"]/1e5, df_trap_force["force_pN"], lw=1)

plt.xlabel("time (ns)")
plt.ylabel("force (pN)")
plt.title(f"{run_id}: Force vs. Time")
plt.tight_layout()
plt.savefig(img_path)
plt.close()