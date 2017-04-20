"""
This script will make csv files that can be used for further analysis in JASP
(or another statistics software)

@author: mje
@email: mads [] cnru.dk
"""
import pandas as pd

# Load data
anova_sim = pd.read_csv("./data/anova_dt.csv")
mxm_sim = pd.read_csv("./data/mxm_dt.csv")

# Make dataframe for anova with only bonferroni corrected p-values
aov_data = pd.DataFrame()
aov_data["sig"] = anova_sim.Bonfsig
aov_data["test"] = "aov"

# Make dataframe for mixed models with only bonferroni corrected p-values
mxm_data = pd.DataFrame()
mxm_data["sig"] = mxm_sim.Bonfsig
mxm_data["test"] = "mxm"

# Save data
both_data = aov_data.append(mxm_data, ignore_index=True)
both_data.to_csv("./data/df_both_long_corrected.csv", index=False)

# Make dataframe for anova with only any significant p-values
aov_data = pd.DataFrame()
aov_data["sig"] = anova_sim.anysig
aov_data["test"] = "aov"

# Make dataframe for anova with only any significant p-values
mxm_data = pd.DataFrame()
mxm_data["sig"] = mxm_sim.anysig
mxm_data["test"] = "mxm"

# Save data
both_data = aov_data.append(mxm_data, ignore_index=True)
both_data.to_csv("./data/df_both_long_uncorrected.csv", index=False)

# Make dataframe for anova with only any significant p-values
aov_data = pd.DataFrame()
aov_data["anysig"] = anova_sim.anysig
aov_data["bonf"] = anova_sim.Bonfsig
aov_data["test"] = "aov"

# Make dataframe for anova with only any significant p-values
mxm_data = pd.DataFrame()
mxm_data["anysig"] = mxm_sim.anysig
mxm_data["bonf"] = mxm_sim.Bonfsig
mxm_data["test"] = "mxm"

# Save data
both_data = aov_data.append(mxm_data, ignore_index=True)
both_data.to_csv("./data/df_both_long.csv", index=False)
