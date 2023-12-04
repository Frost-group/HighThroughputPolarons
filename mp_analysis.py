#%%
from mp_api.client import MPRester
import pandas as pd
import matplotlib.pyplot as plt
#%%
df = pd.read_csv('Data/merged_files/mobility_1_material.tsv',sep='\t')
# %%

with MPRester("tyls1Sy1dIhvUoxQn6NGyg4UXf9qLa1q") as mpr:
    docs = mpr.summary.search(material_ids=list(df["Name"]))
#%%
full_data = []
for i in range(len(docs)):
    example_doc = docs[i]
    mpid = str(example_doc.material_id)
    crystal = str(example_doc.symmetry.crystal_system)
    formula = list(example_doc.elements)
    str_formula = []
    for i in formula:
        str_formula.append("Element" +" " + str(i))
    full_data.append([mpid, crystal, str_formula])


#%%
print(full_data[1][2][1])
# %%
element_groups = {
    'Group13': ['Element B', 'Element Al', 'Element Ga', 'Element In', 'Element Tl'],
    'Group14': ['Element C', 'Element Si', 'Element Ge', 'Element Sn', 'Element Pb'],
    'Group15': ['Element N', 'Element P', 'Element As', 'Element Sb', 'Element Bi'],
    'Group16': ['Element O', 'Element S', 'Element Se', 'Element Te', 'Element Po'],
    'Group17': ['Element F', 'Element Cl', 'Element Br', 'Element I', 'Element At'],
}
def categorize_materials(materials):
    categorized_materials = []

    for material in materials:
        material_name, crystal_structure, material_elements = material
        group = "other"

        for group_name, elements_in_group in element_groups.items():
            if any(element in elements_in_group for element in material_elements):
                group = group_name
                break

        categorized_materials.append([material_name, crystal_structure, group])

    return categorized_materials


# Categorize materials
result = categorize_materials(full_data)

# Display the result
for material_info in result:
    print(f"Material: {material_info[0]}, Crystal Structure: {material_info[1]}, Group: {material_info[2]}")
# %%
df2 = pd.DataFrame(result, columns=["Name", "structure", "group"])

# Merge DataFrames on 'common_column'
merged_df = pd.merge(df, df2, on='Name')
merged_df.to_csv('Data/merged_files/mobility_1_material_crystal.tsv', sep='\t', index=False)
# %%
df = pd.read_csv('Data/merged_files/mobility_1_material_crystal.tsv',sep='\t')
#%%
unique_groups = df['group'].unique()

# Set up a color map for each group
colors = plt.cm.get_cmap('tab10', len(unique_groups))

# Plot each group with a different color
fig, ax = plt.subplots()

for i, group in enumerate(unique_groups[1:5]):
    print(i, group)
    group_data = df[df['group'] == group]
    ax.set_xlabel('Alpha')
    ax.set_ylabel('Mobility')
    
    plt.xlim(0.001, 100)  
    plt.ylim(0.00001, 10000000)
    plt.xscale('log')
    plt.yscale('log')
    ax.scatter(group_data['Alpha'], group_data['Mobility'], label=group, color=colors(i), alpha=0.7, marker='o', s= 10)
    #plt.show()
ax.legend(unique_groups)
# Add labels and legend

#plt.scatter(mobility_data["Alpha"], mobility_data["Mobility"])

# %%


for i, group in enumerate(unique_groups[:2]):
    fig, ax = plt.subplots()
    print(i, group)
    group_data = df[df['group'] == group]
    ax.set_xlabel('Alpha')
    ax.set_ylabel('Mobility')
   
    plt.xlim(0.001, 100)  
    plt.ylim(0.00001, 10000000)
    plt.xscale('log')
    plt.yscale('log')
    ax.scatter(group_data['Alpha'], group_data['Mobility'], label=group, color=colors(i), alpha=0.7, marker='o', s= 10)
    ax.legend(labels = [group])
    #plt.show()

# %%
unique_groups = df['structure'].unique()

# Set up a color map for each group
colors = plt.cm.get_cmap('tab10', len(unique_groups))

# Plot each group with a different color
fig, ax = plt.subplots()

for i, group in enumerate(unique_groups):
    print(i, group)
    group_data = df[df['structure'] == group]
    ax.set_xlabel('Alpha')
    ax.set_ylabel('Mobility')
    
    plt.xlim(0.001, 100)  
    plt.ylim(0.00001, 10000000)
    plt.xscale('log')
    plt.yscale('log')
    ax.scatter(group_data['Alpha'], group_data['Mobility'], label=group, color=colors(i), alpha=0.7, marker='o', s= 10)
    #plt.show()
ax.legend(unique_groups)
# Add labels and legend

#plt.scatter(mobility_data["Alpha"], mobility_data["Mobility"])

# %%


for i, group in enumerate(unique_groups):
    fig, ax = plt.subplots()
    print(i, group)
    group_data = df[df['structure'] == group]
    ax.set_xlabel('Alpha')
    ax.set_ylabel('Mobility')

    plt.xlim(0.001, 100)  
    plt.ylim(0.00001, 10000000)
    plt.xscale('log')
    plt.yscale('log')
    ax.scatter(group_data['Alpha'], group_data['Mobility'], label=group, color=colors(i), alpha=0.7, marker='o', s= 10)
    ax.legend(labels = [group])
    plt.show()
# %%
