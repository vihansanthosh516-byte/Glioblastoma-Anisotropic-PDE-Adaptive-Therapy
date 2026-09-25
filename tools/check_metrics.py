import json

with open("output/adaptive_geometry_metrics.json") as f:
    a = json.load(f)

with open("output/adaptive_geometry_metrics_80.json") as f:
    b = json.load(f)

print(f"Identical: {a == b}")
print()
print("metrics.json means:")
print(f"  drug %:   {sum(p['drug_percent_of_mtd'] for p in a)/len(a):.2f}")
print(f"  VR:       {sum(p['vr_ratio'] for p in a)/len(a):.3f}")
print(f"  holidays: {sum(p['num_holidays'] for p in a)/len(a):.2f}")
print()
print("metrics_80.json means:")
print(f"  drug %:   {sum(p['drug_percent_of_mtd'] for p in b)/len(b):.2f}")
print(f"  VR:       {sum(p['vr_ratio'] for p in b)/len(b):.3f}")
print(f"  holidays: {sum(p['num_holidays'] for p in b)/len(b):.2f}")