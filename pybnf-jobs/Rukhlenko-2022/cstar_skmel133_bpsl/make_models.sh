#!/bin/sh
# Generate the three dose-baked model variants from the baseline model.
# A BPSL constraint experiment can't apply a non-free-parameter condition, so (following
# the Miller2025 example's one-model-per-condition pattern) each inhibitor dose is baked
# into a copy of the model via its I_<x>_conc parameters. Doses are the published values.
base=cstar_skmel133_bpsl.bngl
sed -e 's/^I_AKT_conc\t= 0.0/I_AKT_conc\t= 3.0/'  "$base" > cstar_skmel133_bpsl_akt.bngl
sed -e 's/^I_IRS_conc\t= 0.0/I_IRS_conc\t= 4.0/'  "$base" > cstar_skmel133_bpsl_irs.bngl
sed -e 's/^I_IRS_conc\t= 0.0/I_IRS_conc\t= 2.0/' -e 's/^I_AKT_conc\t= 0.0/I_AKT_conc\t= 1.5/' "$base" > cstar_skmel133_bpsl_combo.bngl
echo "generated akt (I_AKT=3), irs (I_IRS=4), combo (I_IRS=2,I_AKT=1.5) from $base"
