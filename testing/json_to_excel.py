import pandas as pd
import json 

def jsonToPd(filename):
    with open(filename, "r") as f:
        data = json.load(f)

    all_data = []
    for entry in data:
        for kinase_id,inhibitors in entry.items():
            for inhibitor,details in inhibitors.items():
                row = {
                    "Kinase" : kinase_id,
                    "Inhibitor" : inhibitor,
                    **details
                    }
                all_data.append(row)

    df = pd.DataFrame(all_data)
    #df["FC"] = pd.to_numeric(df["FC"])
    return df


df = jsonToPd("KSEA_MTOR.json")
df1 = jsonToPd("KSEA_MTOR_sid.json")
df2 = jsonToPd("KSEA_MTOR_high.json")
df3 = jsonToPd("KSEA_MTOR_sid_high.json")

with pd.ExcelWriter("KSEA1.xlsx",engine="openpyxl",mode="a",if_sheet_exists="new") as writer:
    df.to_excel(writer, sheet_name="MTOR", index=False)
    df1.to_excel(writer, sheet_name="MTOR_sid", index=False)
    df2.to_excel(writer, sheet_name="MTOR_high", index=False)
    df3.to_excel(writer, sheet_name="MTOR_sid_high", index=False)