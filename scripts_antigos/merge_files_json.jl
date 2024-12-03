import Pkg; Pkg.activate("/home/m3g/Documentos/1_analysis/1_analises/p100_analysis")
using ComplexMixtures, PDBTools
const CM = ComplexMixtures

include("programa/find_json.jl")
include("programa/read_data.jl")
includet("programa/myplots.jl")

save_folder = "json_median"

files_glu_full_05 = [
    load(find_json("1.1_100_s1_05/", "gluS*full.json")),
    load(find_json("2.1_100_s2_05/", "gluS*full.json")),
    load(find_json("3.1_100_s3_05/", "gluS*full.json"))
]

files_water_full_05  = [
    load(find_json("1.1_100_s1_05", "water*full.json")),
    load(find_json("2.1_100_s2_05", "water*full.json")),
    load(find_json("3.1_100_s3_05", "water*full.json"))
]

files_glu_compact_05 = [
    load(find_json("1.1_100_s1_05/", "gluS*compact.json")),
    load(find_json("2.1_100_s2_05/", "gluS*compact.json")),
    load(find_json("3.1_100_s3_05/", "gluS*compact.json"))
]

files_glu_extended_05   = [
    load(find_json("1.1_100_s1_05","glu*extended.json")),
    load(find_json("2.1_100_s2_05","glu*extended.json")),
    load(find_json("3.1_100_s3_05","glu*extended.json"))
]

files_water_compact_05  = [
    load(find_json("1.1_100_s1_05", "water*compact.json")),
    load(find_json("2.1_100_s2_05", "water*compact.json")),
    load(find_json("3.1_100_s3_05", "water*compact.json"))
]

files_water_extended_05 = [
    load(find_json("1.1_100_s1_05", "water*extended.json")),
    load(find_json("2.1_100_s2_05", "water*extended.json")),
    load(find_json("3.1_100_s3_05", "water*extended.json"))
]

files_glu_full_75    = [
    load(find_json("1.2_100_s1_75", "glu*full.json")),
    load(find_json("2.2_100_s2_75", "glu*full.json")),
    load(find_json("3.2_100_s3_75", "glu*full.json"))
]

files_water_full_75 = [
    load(find_json("1.2_100_s1_75", "water*full.json")),
    load(find_json("2.2_100_s2_75", "water*full.json")),
    load(find_json("3.2_100_s3_75", "water*full.json"))
]

files_glu_compact_75    = [
    load(find_json("1.2_100_s1_75", "glu*compact.json")),
    load(find_json("2.2_100_s2_75", "glu*compact.json")),
    load(find_json("3.2_100_s3_75", "glu*compact.json"))
]

files_glu_extended_75   = [
    load(find_json("1.2_100_s1_75", "glu*extended.json")),
    load(find_json("2.2_100_s2_75", "glu*extended.json")),
    load(find_json("3.2_100_s3_75", "glu*extended.json"))
]

files_water_compact_75  = [
    load(find_json("1.2_100_s1_75", "water*compact.json")),
    load(find_json("2.2_100_s2_75", "water*compact.json")),
    load(find_json("3.2_100_s3_75", "water*compact.json"))
]

files_water_extended_75 = [
    load(find_json("1.2_100_s1_75", "water*extended.json")),
    load(find_json("2.2_100_s2_75", "water*extended.json")),
    load(find_json("3.2_100_s3_75", "water*extended.json"))
]

files_glu_full_1    = [
    load(find_json("1.3_100_s1_1", "glu*full.json")),
    load(find_json("2.3_100_s2_1", "glu*full.json")),
    load(find_json("3.3_100_s3_1", "glu*full.json"))
]

files_water_full_1 = [
    load(find_json("1.3_100_s1_1", "water*full.json")),
    load(find_json("2.3_100_s2_1", "water*full.json")),
    load(find_json("3.3_100_s3_1", "water*full.json"))
]

files_glu_compact_1    = [
    load(find_json("1.3_100_s1_1", "glu*compact.json")),
    load(find_json("2.3_100_s2_1", "glu*compact.json")),
    load(find_json("3.3_100_s3_1", "glu*compact.json"))
]

files_glu_extended_1   = [
    load(find_json("1.3_100_s1_1", "glu*extended.json")),
    load(find_json("2.3_100_s2_1", "glu*extended.json")),
    load(find_json("3.3_100_s3_1", "glu*extended.json"))
]

files_water_compact_1  = [
    load(find_json("1.3_100_s1_1", "water*compact.json")),
    load(find_json("2.3_100_s2_1", "water*compact.json")),
    load(find_json("3.3_100_s3_1", "water*compact.json"))
]

files_water_extended_1 = [
    load(find_json("1.3_100_s1_1", "water*extended.json")),
    load(find_json("2.3_100_s2_1", "water*extended.json")),
    load(find_json("3.3_100_s3_1", "water*extended.json"))
]

# Criando os vetores
vec_glu_compact_05, vec_glu_extended_05, vec_water_compact_05, vec_water_extended_05 = [],[],[],[]
vec_glu_compact_75, vec_glu_extended_75, vec_water_compact_75, vec_water_extended_75 = [],[],[],[]
vec_glu_compact_1, vec_glu_extended_1, vec_water_compact_1, vec_water_extended_1     = [],[],[],[]

# Mesclando arquivos para compactos e estendidos, para diferentes concentrações

#sistema glu full
vec_glu_full_05       = merge(files_glu_full_05)
vec_glu_full_75       = merge(files_glu_full_75)
vec_glu_full_1        = merge(files_glu_full_1)

#sistema water full
vec_water_full_05     = merge(files_water_full_05)
vec_water_full_75     = merge(files_water_full_75)
vec_water_full_1      = merge(files_water_full_1)

# Concentração 0.50
vec_glu_compact_05    = merge(files_glu_compact_05)
vec_glu_extended_05   = merge(files_glu_extended_05)
vec_water_compact_05  = merge(files_water_compact_05)
vec_water_extended_05 = merge(files_water_extended_05)

# Concentração 0.75
vec_glu_compact_75    = merge(files_glu_compact_75)
vec_glu_extended_75   = merge(files_glu_extended_75)
vec_water_compact_75  = merge(files_water_compact_75)
vec_water_extended_75 = merge(files_water_extended_75)

# Concentração 1.0
vec_glu_compact_1     = merge(files_glu_compact_1)
vec_glu_extended_1    = merge(files_glu_extended_1)
vec_water_compact_1   = merge(files_water_compact_1)
vec_water_extended_1  = merge(files_water_extended_1)

# Salvando os resultados mesclados em arquivos JSON
# Analise Full glu
save(vec_glu_full_05,      "$save_folder/glu_full_05_average.json")
save(vec_glu_full_75,      "$save_folder/glu_full_75_average.json")
save(vec_glu_full_1,      "$save_folder/glu_full_1_average.json")

# Analise Full water
save(vec_water_full_05,      "$save_folder/water_full_05_average.json")
save(vec_water_full_75,      "$save_folder/water_full_75_average.json")
save(vec_water_full_1,      "$save_folder/water_full_1_average.json")

# Concentração 0.05
save(vec_glu_compact_05,    "$save_folder/glu_compact_05_average.json")
save(vec_glu_extended_05,   "$save_folder/glu_extended_05_average.json")
save(vec_water_compact_05,  "$save_folder/water_compact_05_average.json")
save(vec_water_extended_05, "$save_folder/water_extended_05_average.json")

# Concentração 0.75
save(vec_glu_compact_75,    "$save_folder/glu_compact_75_average.json")
save(vec_glu_extended_75,   "$save_folder/glu_extended_75_average.json")
save(vec_water_compact_75,  "$save_folder/water_compact_75_average.json")
save(vec_water_extended_75, "$save_folder/water_extended_75_average.json")

# Concentração 1.0
save(vec_glu_compact_1,    "$save_folder/glu_compact_1_average.json")
save(vec_glu_extended_1,   "$save_folder/glu_extended_1_average.json")
save(vec_water_compact_1,  "$save_folder/water_compact_1_average.json")
save(vec_water_extended_1, "$save_folder/water_extended_1_average.json")

# Carregue os dados
glu_full_05, glu_full_75, glu_full_1, water_full_05, water_full_75, 
water_full_1, glu_compact_05, glu_extended_05, water_compact_05,
water_extended_05, glu_compact_75, glu_extended_75, water_compact_75,
water_extended_75, glu_compact_1, glu_extended_1, water_compact_1,
water_extended_1 = read_data()

# Chame a função de plotagem
myplots(glu_full_05, glu_full_75, glu_full_1, water_full_05, water_full_75, 
water_full_1, glu_compact_05, glu_extended_05, water_compact_05,
water_extended_05, glu_compact_75, glu_extended_75, water_compact_75,
water_extended_75, glu_compact_1, glu_extended_1, water_compact_1,
water_extended_1)



