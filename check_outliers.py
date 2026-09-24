import json

with open("output/adaptive_geometry_metrics.json") as f:
    m = json.load(f)

print("=" * 60)
print("VR ratio < 1.0 (adaptive wins on final mass)")
print("=" * 60)
count = 0
for p in m:
    if p["vr_ratio"] < 1.0:
        count += 1
        print(f"{p['patient_id']}: VR={p['vr_ratio']:.3f}, "
              f"rho={p['rho_per_day']:.5f}, "
              f"drug%={p['drug_percent_of_mtd']:.1f}, "
              f"holidays={p['num_holidays']}, "
              f"rf_mtd={p['resistant_fraction_mtd']:.3f}, "
              f"rf_ad={p['resistant_fraction_adaptive']:.3f}")
if count == 0:
    print("(none)")

print()
print("=" * 60)
print("Adaptive resistance fraction > 0.95")
print("=" * 60)
count = 0
for p in m:
    if p["resistant_fraction_adaptive"] > 0.95:
        count += 1
        print(f"{p['patient_id']}: rf_ad={p['resistant_fraction_adaptive']:.3f}, "
              f"rf_mtd={p['resistant_fraction_mtd']:.3f}, "
              f"rho={p['rho_per_day']:.5f}, "
              f"holidays={p['num_holidays']}, "
              f"drug%={p['drug_percent_of_mtd']:.1f}")
if count == 0:
    print("(none)")

print()
print("=" * 60)
print("Adaptive progressed earlier than MTD (ttp_adaptive < ttp_mtd)")
print("=" * 60)
count = 0
for p in m:
    if p["ttp_adaptive"] < p["ttp_mtd"]:
        count += 1
        print(f"{p['patient_id']}: ttp_mtd={p['ttp_mtd']}, "
              f"ttp_ad={p['ttp_adaptive']}, "
              f"rho={p['rho_per_day']:.5f}, "
              f"vr={p['vr_ratio']:.2f}")
if count == 0:
    print("(none)")